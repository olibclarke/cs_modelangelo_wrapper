# cs_modelangelo_wrapper

`cs_modelangelo_wrapper.py` runs [ModelAngelo](https://www.nature.com/articles/s41586-024-07215-4) (Jiamali et al., 2024) using a CryoSPARC job as input and records the run as a CryoSPARC External job.

## Requirements

- `cryosparc-tools` in the Python environment running the wrapper
- A readable CryoSPARC `instance_info.json`
- [ModelAngelo](https://github.com/3dem/model-angelo) available as:

```bash
relion_python_modelangelo
```

The wrapper checks that the executable is available before running.

## Basic usage

Run ModelAngelo with a protein FASTA:

```bash
python3 cs_modelangelo_wrapper.py P27 W3 J1486 --sequence protein.fasta
```

Run ModelAngelo without a sequence:

```bash
python3 cs_modelangelo_wrapper.py P27 W3 J1486
```

Specify a GPU:

```bash
python3 cs_modelangelo_wrapper_labeled_cifs_sequence_output.py P27 W3 J1486 --sequence protein.fasta --device 0
```

Pass extra ModelAngelo options after `--`:

```bash
python3 cs_modelangelo_wrapper_labeled_cifs_sequence_output.py P27 W3 J1486 --sequence protein.fasta --device 0 -- --keep-intermediate-results
```

## Inputs

The positional arguments are:

```text
project workspace job
```

where `job` is a CryoSPARC job containing a volume output.

By default, the script infers the input volume group and preferentially uses the sharpened map when available. You can override this with:

```bash
--source-group volume
--map-field map/path
```

Optional sequence inputs:

```bash
--sequence protein.fasta
--rna-fasta rna.fasta
--dna-fasta dna.fasta
```

If `--sequence` is provided, the wrapper runs:

```bash
relion_python_modelangelo build
```

If no sequence is provided, it runs:

```bash
relion_python_modelangelo build_no_seq
```

## Masks

If a mask is present in the input volume output, the wrapper tries to use it automatically. Override or disable this behavior with:

```bash
--mask mask.mrc
--mask-field mask/path
--no-mask
```

## Outputs

The External job directory contains:

```text
modelangelo_cif_files/
modelangelo_results.tar.gz
```

The wrapper registers downloadable outputs for:

- each generated `.cif` file (labeled as a volume, as a workaround)
- the complete `modelangelo_results.tar.gz` archive
- input FASTA files, when supplied

CIFs found inside an `entropy_score` folder are renamed with an `_entropy_score.cif` suffix.

## Instance info convention

If `--instance-info` is not supplied, the wrapper searches:

```text
~/instance_info.json
```

Use an explicit path if needed:

```bash
--instance-info /path/to/instance_info.json
```

## Common options

```text
--source-group GROUP       Input volume output group; inferred if omitted
--map-field FIELD          Map path field; sharpened maps are preferred by default
--row N                    Row of the volume dataset to use; default 0
--sequence FASTA           Protein FASTA file
--rna-fasta FASTA          Optional RNA FASTA file
--dna-fasta FASTA          Optional DNA FASTA file
--mask MRC                 Explicit mask path
--no-mask                  Do not use an input mask
--device DEVICE            ModelAngelo device, e.g. 0, 1, cpu
--executable PATH          ModelAngelo executable; default relion_python_modelangelo
--output NAME              Base CryoSPARC output name for CIF files
--project-dir PATH         Override CryoSPARC project directory
--instance-info PATH       Path to instance_info.json
```
