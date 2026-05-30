#!/usr/bin/env python3
"""
Run ModelAngelo as a CryoSPARC External job.

The script takes a CryoSPARC job containing a volume output, optionally takes
protein/RNA/DNA FASTA files, runs relion_python_modelangelo, streams the
ModelAngelo stdout/stderr live into the CryoSPARC job log, and registers the
resulting CIF files and input sequence files as downloadable CryoSPARC output groups.
A tar.gz archive of the full ModelAngelo output directory is always written in the External job directory.

Examples:

    python3 cs_modelangelo_external.py P44 W2 J123 \
        --source-group volume \
        --sequence sequence.fasta \
        --device 0

    python3 cs_modelangelo_external.py P44 W2 J123 \
        --source-group volume \
        --device cpu

Additional ModelAngelo arguments can be appended after --:

    python3 cs_modelangelo_external.py P44 W2 J123 --sequence seq.fasta -- \
        --keep-intermediate-results
"""

from __future__ import annotations

import argparse
import inspect
import json
import os
import shlex
import shutil
import subprocess
import sys
import tarfile
import time
import warnings
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional, Sequence, Tuple

import numpy as np

DEFAULT_EXECUTABLE = "relion_python_modelangelo"
DEFAULT_SOURCE_GROUPS = (
    # Prefer sharpened outputs/groups when the source job exposes them.
    "map_sharp",
    "volume_sharp",
    "volume",
    "volume_class_0",
    "volume_class_1",
    "volume_class_2",
    "volume_class_3",
    "map",
    "volume_masked",
)
MAP_FIELD_CANDIDATES = (
    # Prefer sharpened maps when present inside a volume output.
    "map_sharp/path",
    "volume_sharp/path",
    "sharp_map/path",
    "map/path",
    "volume/path",
    "volume_blob/path",
)
MASK_FIELD_CANDIDATES = (
    "mask/path",
    "mask_fsc/path",
    "mask_refine/path",
    "volume_mask/path",
)
MODEL_OUTPUT_NAME = "modelangelo_cif"
ARCHIVE_NAME = "modelangelo_results.tar.gz"


def norm_text(x: object) -> str:
    if x is None:
        return ""
    if isinstance(x, bytes):
        return x.decode("utf-8", errors="replace")
    if isinstance(x, np.bytes_):
        return x.decode("utf-8", errors="replace")
    return str(x).strip()


def shell_join(cmd: Sequence[str]) -> str:
    return " ".join(shlex.quote(str(x)) for x in cmd)


def existing_field(ds, field: str) -> bool:
    try:
        ds[field]
        return True
    except Exception:
        return False


def field_names(ds) -> List[str]:
    try:
        return sorted(norm_text(f) for f in ds.fields())
    except Exception:
        return []


def dataset_slots(ds) -> List[str]:
    try:
        return sorted(p for p in ds.prefixes() if p != "uid")
    except Exception:
        out = set()
        for f in field_names(ds):
            if "/" in f:
                out.add(f.split("/", 1)[0])
        return sorted(out)


def safe_obj_dir(obj) -> Path:
    for attr_name in ("dir", "project_dir", "path"):
        attr = getattr(obj, attr_name, None)
        if attr is None:
            continue
        try:
            value = attr() if callable(attr) else attr
            if value:
                return Path(norm_text(value)).expanduser().resolve()
        except Exception:
            pass
    raise RuntimeError(f"Could not infer filesystem directory for {obj!r}")


def load_instance_info(path: Optional[str]) -> Tuple[Dict[str, Any], Path]:
    candidates: List[Path] = []
    if path:
        candidates.append(Path(path).expanduser())
    candidates.extend(
        [
            Path("/home/user/instance_info.json"),
            Path("/home/exx/instance_info.json"),
            Path.home() / "instance_info.json",
        ]
    )
    seen = set()
    for candidate in candidates:
        candidate = candidate.expanduser()
        key = str(candidate)
        if key in seen:
            continue
        seen.add(key)
        if candidate.is_file():
            try:
                return json.loads(candidate.read_text()), candidate
            except Exception as exc:
                raise RuntimeError(f"Failed to read instance_info.json from {candidate}") from exc
    raise FileNotFoundError(
        "Could not find instance_info.json. Pass --instance-info explicitly, "
        "or place it at ~/instance_info.json or /home/user/instance_info.json."
    )


def connect_cryosparc(instance_info_path: Optional[str]):
    with warnings.catch_warnings():
        warnings.simplefilter("default")
        from cryosparc.tools import CryoSPARC

        info, info_path = load_instance_info(instance_info_path)
        cs = CryoSPARC(**info)
    if hasattr(cs, "test_connection") and not cs.test_connection():
        raise RuntimeError(f"Could not connect to CryoSPARC using {info_path}")
    return cs, str(info_path)


def find_project_and_job(cs, project_uid: str, job_uid: str):
    project = cs.find_project(project_uid)
    source_job = project.find_job(job_uid)
    return project, source_job


def create_external_job(project, workspace_uid: str, title: str, description: str):
    kwargs = {"workspace_uid": workspace_uid, "title": title}
    try:
        params = set(inspect.signature(project.create_external_job).parameters)
    except Exception:
        params = set()
    if description:
        if "desc" in params:
            kwargs["desc"] = description
        elif "description" in params:
            kwargs["description"] = description
    try:
        return project.create_external_job(**kwargs)
    except TypeError:
        # Older cryosparc-tools versions may take workspace UID positionally.
        if description:
            try:
                return project.create_external_job(workspace_uid, title=title, desc=description)
            except TypeError:
                pass
        return project.create_external_job(workspace_uid, title=title)


def infer_project_dir(project, override: Optional[str]) -> Path:
    if override:
        p = Path(override).expanduser().resolve()
        if not p.exists():
            raise FileNotFoundError(f"--project-dir does not exist: {p}")
        return p
    return safe_obj_dir(project)


def project_path(project_dir: Path, path_text: object) -> Path:
    p = Path(norm_text(path_text))
    if p.is_absolute():
        return p
    return project_dir / p


def relative_to_project(project_dir: Path, path: Path) -> str:
    try:
        return str(path.resolve().relative_to(project_dir.resolve()))
    except Exception:
        return str(path.resolve())


def load_source_volume(source_job, requested_group: Optional[str]):
    groups: List[str] = []
    if requested_group:
        groups.append(requested_group)
    for g in DEFAULT_SOURCE_GROUPS:
        if g not in groups:
            groups.append(g)
    last_exc = None
    for g in groups:
        try:
            return g, source_job.load_output(g)
        except Exception as exc:
            last_exc = exc
    raise RuntimeError(f"Could not load a volume output from {source_job.uid}. Tried: {', '.join(groups)}") from last_exc


def detect_field(ds, explicit: Optional[str], candidates: Sequence[str], kind: str, required: bool) -> Optional[str]:
    fields = set(field_names(ds))
    if explicit:
        if explicit in fields or existing_field(ds, explicit):
            return explicit
        raise RuntimeError(f"Requested {kind} field is not present: {explicit}")
    for f in candidates:
        if f in fields or existing_field(ds, f):
            return f
    if not required:
        return None
    path_fields = sorted(f for f in fields if f.endswith("/path"))
    raise RuntimeError(
        f"Could not detect {kind} path field. Available path fields: "
        + (", ".join(path_fields) if path_fields else "none")
    )


def first_path_from_field(ds, field: str, project_dir: Path, row: int = 0) -> Path:
    arr = np.asarray(ds[field])
    if arr.ndim == 0:
        value = arr.item()
    else:
        if len(arr) <= row:
            raise RuntimeError(f"Field {field} has {len(arr)} rows, cannot use row {row}")
        value = arr[row]
    p = project_path(project_dir, value)
    if not p.exists():
        raise FileNotFoundError(f"File referenced by {field} does not exist: {p}")
    return p


def connect_input_volume(external_job, source_job_uid: str, source_group: str, slots: Sequence[str]) -> None:
    try:
        external_job.connect(
            target_input="volume",
            source_job_uid=source_job_uid,
            source_output=source_group,
            slots=list(slots),
        )
        return
    except TypeError:
        pass
    external_job.add_input(type="volume", name="volume", slots=list(slots))
    external_job.connect("volume", source_job_uid, source_group)


def check_modelangelo(executable: str) -> str:
    exe = shutil.which(executable) if os.sep not in executable else executable
    if not exe:
        raise RuntimeError(
            f"Could not find ModelAngelo executable {executable!r} on PATH. "
            "Pass --executable /path/to/relion_python_modelangelo if needed."
        )
    try:
        proc = subprocess.run(
            [exe, "--version"],
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
            timeout=30,
        )
        version_text = (proc.stdout or "").strip()
    except Exception:
        version_text = ""
    return exe if not version_text else f"{exe}\n{version_text}"


def build_modelangelo_command(
    executable_path: str,
    volume_path: Path,
    output_dir: Path,
    protein_fasta: Optional[Path],
    rna_fasta: Optional[Path],
    dna_fasta: Optional[Path],
    mask_path: Optional[Path],
    device: Optional[str],
    keep_intermediate_results: bool,
    extra_args: Sequence[str],
) -> List[str]:
    executable = executable_path.splitlines()[0]
    has_sequence = bool(protein_fasta or rna_fasta or dna_fasta)
    if protein_fasta is None and (rna_fasta or dna_fasta):
        # ModelAngelo's build command technically requires protein-fasta in this version.
        # Use no-sequence mode rather than constructing an invalid build command.
        has_sequence = False
    subcommand = "build" if has_sequence else "build_no_seq"
    cmd = [
        executable,
        subcommand,
        "--volume-path",
        str(volume_path),
        "--output-dir",
        str(output_dir),
    ]
    if has_sequence and protein_fasta is not None:
        cmd += ["--protein-fasta", str(protein_fasta)]
    if has_sequence and rna_fasta is not None:
        cmd += ["--rna-fasta", str(rna_fasta)]
    if has_sequence and dna_fasta is not None:
        cmd += ["--dna-fasta", str(dna_fasta)]
    if mask_path is not None:
        cmd += ["--mask-path", str(mask_path)]
    if device:
        cmd += ["--device", str(device)]
    if keep_intermediate_results:
        cmd += ["--keep-intermediate-results"]
    cmd += list(extra_args)
    return cmd


def stream_command_to_job_log(cmd: Sequence[str], job, cwd: Optional[Path] = None, env: Optional[Dict[str, str]] = None) -> int:
    proc = subprocess.Popen(
        list(cmd),
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        stdin=subprocess.DEVNULL,
        text=True,
        bufsize=1,
        universal_newlines=True,
        cwd=str(cwd) if cwd else None,
        env=env,
    )
    assert proc.stdout is not None
    for line in proc.stdout:
        text = line.rstrip("\n")
        if text:
            try:
                job.log(text)
            except Exception:
                print(text, flush=True)
    return int(proc.wait())


def copy_sequence_file(path_text: Optional[str], work_dir: Path) -> Optional[Path]:
    if not path_text:
        return None
    src = Path(path_text).expanduser().resolve()
    if not src.exists():
        raise FileNotFoundError(f"Sequence file does not exist: {src}")
    dest = work_dir / src.name
    if src != dest:
        shutil.copy2(src, dest)
    return dest


def collect_files(root: Path, suffixes: Sequence[str]) -> List[Path]:
    lowered = tuple(s.lower() for s in suffixes)
    out: List[Path] = []
    if root.exists():
        for p in root.rglob("*"):
            if p.is_file() and p.name.lower().endswith(lowered):
                out.append(p)
    return sorted(out)


def classify_modelangelo_file_source(src: Path) -> Optional[str]:
    """Return a short label for special ModelAngelo output subdirectories.

    ModelAngelo often writes CIFs with duplicate basenames in different output
    subdirectories. A duplicate copied as ``foo_1.cif`` is ambiguous in the
    CryoSPARC output panel, so label known special directories explicitly.
    """
    parts = [p.lower() for p in src.parts]
    if "entropy_score" in parts or "entropy-score" in parts:
        return "entropy_score"
    return None


def destination_name_for_modelangelo_output(src: Path, used_names: Dict[str, int]) -> str:
    label = classify_modelangelo_file_source(src)
    if label:
        name = f"{src.stem}_{label}{src.suffix}"
        if name not in used_names:
            used_names[name] = 1
            return name

    name = src.name
    count = used_names.get(name, 0)
    used_names[name] = count + 1
    if count:
        name = f"{src.stem}_{count}{src.suffix}"
    return name


def copy_outputs_to_job_dir(files: Sequence[Path], dest_dir: Path) -> List[Path]:
    dest_dir.mkdir(parents=True, exist_ok=True)
    copied: List[Path] = []
    used_names: Dict[str, int] = {}
    for src in files:
        name = destination_name_for_modelangelo_output(src, used_names)
        dest = dest_dir / name
        if src.resolve() != dest.resolve():
            shutil.copy2(src, dest)
        copied.append(dest)
    return copied


def make_archive(source_dir: Path, archive_path: Path) -> None:
    if archive_path.exists():
        archive_path.unlink()
    with tarfile.open(archive_path, "w:gz") as tar:
        tar.add(source_dir, arcname=source_dir.name)




def sanitize_group_name(text: str, fallback: str = "cif") -> str:
    """Return a valid CryoSPARC result-group name."""
    import re

    stem = Path(text).stem if text else fallback
    name = re.sub(r"[^0-9A-Za-z_]", "_", stem)
    name = re.sub(r"_+", "_", name).strip("_")
    if not name:
        name = fallback
    if not name[0].isalpha():
        name = f"{fallback}_{name}"
    return name[:48]


def make_downloadable_map_dataset(path: Path, project_dir: Path):
    """Create a one-row volume/map Dataset pointing to an arbitrary file.

    CryoSPARC does not provide a generic arbitrary-file output type in the
    External-job schema exposed on this installation. `annotation_model` only
    accepts `checkpoint`, which creates a checkpoint `.cs` file rather than a
    directly downloadable CIF. This wrapper therefore registers each CIF as a
    deliberately named one-row `volume` output with its `map/path` pointing to
    the CIF file. Do not connect these outputs to downstream map-processing jobs;
    they are file-download outputs only.
    """
    try:
        from cryosparc.dataset import Dataset  # type: ignore
    except Exception:
        from cryosparc.tools import Dataset  # type: ignore

    rel = relative_to_project(project_dir, path)
    if os.path.isabs(rel):
        raise RuntimeError(
            "ModelAngelo output is not inside the CryoSPARC project directory, "
            f"so it cannot be saved as a portable output group: {path}"
        )
    return Dataset([
        ("uid", np.asarray([1], dtype=np.uint64)),
        ("map/path", np.asarray([rel], dtype=object)),
        ("map/shape", np.asarray([[1, 1, 1]], dtype=np.uint32)),
        ("map/psize_A", np.asarray([1.0], dtype=np.float32)),
    ])


def register_one_downloadable_file(job, path: Path, project_dir: Path, output_name: str, title: str) -> None:
    ds = make_downloadable_map_dataset(path, project_dir)
    errors: List[str] = []
    attempts: List[Tuple[str, Any]] = [
        ("volume/map direct slots", lambda: _register_downloadable_volume_direct(job, output_name, ds, title)),
        ("volume/map explicit slot spec", lambda: _register_downloadable_volume_spec(job, output_name, ds, title)),
        ("volume/map alloc dataset", lambda: _register_downloadable_volume_alloc(job, output_name, ds, title)),
    ]
    for label, fn in attempts:
        try:
            fn()
            job.log(f"Registered downloadable output `{output_name}` for `{path.name}` using {label}.")
            return
        except Exception as exc:
            errors.append(f"{label}: {exc}")

    msg = [
        f"Failed to register `{path.name}` as a downloadable CryoSPARC output group.",
        "The file is still present in the External job directory.",
        "Registration errors:",
    ]
    msg.extend(errors)
    raise RuntimeError("\n".join(msg))


def _register_downloadable_volume_direct(job, output_name: str, ds, title: str) -> None:
    job.add_output(type="volume", name=output_name, slots=["map"], title=title)
    job.save_output(output_name, ds)


def _register_downloadable_volume_spec(job, output_name: str, ds, title: str) -> None:
    job.add_output(
        type="volume",
        name=output_name,
        slots=[{"name": "map", "dtype": "map"}],
        title=title,
    )
    job.save_output(output_name, ds)


def _register_downloadable_volume_alloc(job, output_name: str, ds, title: str) -> None:
    allocated = job.add_output(type="volume", name=output_name, slots=["map"], alloc=ds, title=title)
    job.save_output(output_name, allocated)


def register_model_output(job, model_paths: Sequence[Path], project_dir: Path, output_name: str) -> None:
    """Expose ModelAngelo CIF files as downloadable outputs."""
    if not model_paths:
        raise RuntimeError("No CIF files supplied for output registration")

    used: Dict[str, int] = {}
    for path in model_paths:
        base = sanitize_group_name(path.name, fallback="cif")
        group = f"{output_name}_{base}"
        if len(group) > 60:
            group = group[:60]
        count = used.get(group, 0)
        used[group] = count + 1
        if count:
            suffix = f"_{count + 1}"
            group = f"{group[:60-len(suffix)]}{suffix}"
        register_one_downloadable_file(job, path, project_dir, group, f"CIF: {path.name}")
    job.log(f"Registered {len(model_paths)} CIF file(s) as downloadable output group(s).")


def register_archive_output(job, archive_path: Path, project_dir: Path) -> None:
    register_one_downloadable_file(
        job,
        archive_path,
        project_dir,
        "modelangelo_results_archive",
        "ModelAngelo results archive",
    )


def register_sequence_outputs(job, fasta_paths: Sequence[Tuple[str, Optional[Path]]], project_dir: Path) -> None:
    for label, path in fasta_paths:
        if path is None:
            continue
        group_name = sanitize_group_name(f"modelangelo_input_{label}_{path.name}", fallback="fasta")
        register_one_downloadable_file(
            job,
            path,
            project_dir,
            group_name,
            f"ModelAngelo input {label} FASTA: {path.name}",
        )

def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description="Run ModelAngelo as a CryoSPARC External job from an input volume.",
    )
    parser.add_argument("project", help="CryoSPARC project UID, e.g. P44")
    parser.add_argument("workspace", help="CryoSPARC workspace UID, e.g. W2")
    parser.add_argument("job", help="Input CryoSPARC job UID containing a volume output, e.g. J123")
    parser.add_argument("--source-group", default=None, help="Input volume output group; inferred if omitted")
    parser.add_argument("--map-field", default=None, help="Volume dataset path field to use; default auto-detects sharpened map fields first, then map/path")
    parser.add_argument("--row", type=int, default=0, help="Volume dataset row to use")
    parser.add_argument("--sequence", "--protein-fasta", dest="protein_fasta", default=None, help="Protein FASTA file; if omitted, build_no_seq is used")
    parser.add_argument("--rna-fasta", default=None, help="Optional RNA FASTA file for ModelAngelo build mode")
    parser.add_argument("--dna-fasta", default=None, help="Optional DNA FASTA file for ModelAngelo build mode")
    parser.add_argument("--mask", dest="mask_path", default=None, help="Optional explicit mask MRC path")
    parser.add_argument("--mask-field", default=None, help="Mask dataset path field; default auto-detects mask/path when present")
    parser.add_argument("--no-mask", action="store_true", help="Do not auto-use a mask from the input volume output")
    parser.add_argument("--device", default=None, help="ModelAngelo device argument, e.g. 0, 1, cpu. If omitted, ModelAngelo chooses")
    parser.add_argument("--executable", default=DEFAULT_EXECUTABLE, help="ModelAngelo executable")
    parser.add_argument("--keep-intermediate-results", action="store_true", help="Pass --keep-intermediate-results to ModelAngelo")
    parser.add_argument("--output", default=MODEL_OUTPUT_NAME, help="CryoSPARC output group name for CIF files")
    parser.add_argument("--project-dir", default=None, help="Override CryoSPARC project directory")
    parser.add_argument("--instance-info", default=None, help="Path to instance_info.json; default searches /home/user/instance_info.json, /home/exx/instance_info.json, then ~/instance_info.json")
    return parser


def main(argv: Optional[Sequence[str]] = None) -> int:
    parser = build_parser()
    args, extra_args = parser.parse_known_args(argv)
    if extra_args and extra_args[0] == "--":
        extra_args = extra_args[1:]

    cs, instance_info_used = connect_cryosparc(args.instance_info)
    project_uid = norm_text(args.project).upper()
    workspace_uid = norm_text(args.workspace).upper()
    source_job_uid = norm_text(args.job).upper()
    project, source_job = find_project_and_job(cs, project_uid, source_job_uid)
    project_dir = infer_project_dir(project, args.project_dir)

    source_group, volume_ds = load_source_volume(source_job, args.source_group)
    map_field = detect_field(volume_ds, args.map_field, MAP_FIELD_CANDIDATES, "volume map", required=True)
    volume_path = first_path_from_field(volume_ds, map_field, project_dir, row=int(args.row))

    explicit_mask = Path(args.mask_path).expanduser().resolve() if args.mask_path else None
    mask_path = None
    mask_field = None
    if explicit_mask is not None:
        if not explicit_mask.exists():
            raise FileNotFoundError(f"Explicit mask file does not exist: {explicit_mask}")
        mask_path = explicit_mask
    elif not args.no_mask:
        mask_field = detect_field(volume_ds, args.mask_field, MASK_FIELD_CANDIDATES, "mask", required=False)
        if mask_field:
            try:
                mask_path = first_path_from_field(volume_ds, mask_field, project_dir, row=int(args.row))
            except Exception:
                mask_path = None

    exe_info = check_modelangelo(args.executable)
    exe_path = exe_info.splitlines()[0]

    mode = "build" if args.protein_fasta else "build_no_seq"
    fasta_label = Path(args.protein_fasta).name if args.protein_fasta else "no sequence"
    title = f"ModelAngelo {mode}: {source_group} from {source_job_uid}"
    desc = f"ModelAngelo {mode} on `{source_group}` from {source_job_uid}; {fasta_label}."
    external_job = create_external_job(project, workspace_uid, title=title, description=desc)

    # External jobs created via cryosparc-tools are started/stopped by the run()
    # context manager. Calling queue() is invalid on newer cryosparc-tools.
    with external_job.run():
        job_dir = safe_obj_dir(external_job)
        work_dir = job_dir / "modelangelo_work"
        output_dir = work_dir / "modelangelo_output"
        copied_dir = job_dir / "modelangelo_cif_files"
        work_dir.mkdir(parents=True, exist_ok=True)
        output_dir.mkdir(parents=True, exist_ok=True)

        external_job.log(f"Using CryoSPARC instance_info: {instance_info_used}")

        try:
            connect_input_volume(external_job, source_job_uid, source_group, dataset_slots(volume_ds))
        except Exception as exc:
            external_job.log(f"Warning: could not connect input volume for provenance: {exc}")

        protein_fasta = copy_sequence_file(args.protein_fasta, work_dir)
        rna_fasta = copy_sequence_file(args.rna_fasta, work_dir)
        dna_fasta = copy_sequence_file(args.dna_fasta, work_dir)

        cmd = build_modelangelo_command(
            executable_path=exe_path,
            volume_path=volume_path,
            output_dir=output_dir,
            protein_fasta=protein_fasta,
            rna_fasta=rna_fasta,
            dna_fasta=dna_fasta,
            mask_path=mask_path,
            device=args.device,
            keep_intermediate_results=bool(args.keep_intermediate_results),
            extra_args=extra_args,
        )

        external_job.log(f"ModelAngelo executable: {exe_path}")
        if len(exe_info.splitlines()) > 1:
            for line in exe_info.splitlines()[1:]:
                external_job.log(f"ModelAngelo version: {line}")
        if map_field in {"map_sharp/path", "volume_sharp/path", "sharp_map/path"} or "sharp" in source_group.lower():
            external_job.log("Using sharpened map for ModelAngelo input.")
        external_job.log(f"Input volume: {source_job_uid}:{source_group}:{map_field}[{int(args.row)}]")
        external_job.log(f"Volume path: {volume_path}")
        if mask_path is not None:
            mask_source = args.mask_path if args.mask_path else f"{mask_field}[{int(args.row)}]"
            external_job.log(f"Mask path: {mask_path} ({mask_source})")
        if protein_fasta is not None:
            external_job.log(f"Protein FASTA: {protein_fasta}")
        else:
            external_job.log("No protein FASTA supplied; running build_no_seq.")
        if rna_fasta is not None:
            external_job.log(f"RNA FASTA: {rna_fasta}")
        if dna_fasta is not None:
            external_job.log(f"DNA FASTA: {dna_fasta}")
        if extra_args:
            external_job.log(f"Additional ModelAngelo arguments: {shell_join(extra_args)}")
        external_job.log(f"Running command: {shell_join(cmd)}")

        t0 = time.time()
        ret = stream_command_to_job_log(cmd, external_job, cwd=work_dir)
        dt = time.time() - t0
        if ret != 0:
            raise RuntimeError(f"ModelAngelo failed with exit code {ret}")
        external_job.log(f"ModelAngelo finished successfully in {dt / 60.0:.1f} min.")

        cif_files = collect_files(output_dir, (".cif",))
        if not cif_files:
            raise RuntimeError(f"ModelAngelo completed but no CIF files were found under {output_dir}")
        copied_cifs = copy_outputs_to_job_dir(cif_files, copied_dir)
        for p in copied_cifs:
            external_job.log(f"CIF output file: {p.name}")

        archive_path = job_dir / ARCHIVE_NAME
        make_archive(output_dir, archive_path)
        external_job.log(f"Wrote archive: {archive_path.name}")

        register_model_output(external_job, copied_cifs, project_dir, args.output)
        register_archive_output(external_job, archive_path, project_dir)
        register_sequence_outputs(
            external_job,
            (("protein", protein_fasta), ("rna", rna_fasta), ("dna", dna_fasta)),
            project_dir,
        )

    print(f"Created External job: {project_uid}/{external_job.uid}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
