# CAAARal-Reefs

This repository contains the code, notebooks, data products, and figures for the CAAARal Reefs project.

For MTFC math modelling competition judges, the main point to know is that this project uses Git LFS for large files, but it also includes Google Drive download scripts in case Git LFS is unavailable or incomplete on your machine.

## Repository layout

- `notebooks/cleaning/`: data cleaning notebooks
- `notebooks/analysis/`: analysis notebooks
- `scripts/`: Google Drive download helpers
- `data/`: raw, processed, and analyzed datasets
- `outputs/`: generated maps and figures
- `MATLAB/`: MATLAB models and exported figures

## Option 1: Clone with Git LFS

Use this if `git lfs` is available on your system.

1. Install Git LFS if needed.
2. Clone the repository:

```bash
git clone <REPO_URL>
cd CAAARal-Reefs
```

3. Initialize Git LFS and download large files:

```bash
git lfs install
git lfs pull
```

If Git LFS completes successfully, the repository should contain the large CSV and ZIP assets required by the project.

## Option 2: Clone normally and use the Google Drive download scripts

Use this if:

- `git lfs` is not installed
- Git LFS files appear as small pointer text files instead of real datasets
- `git lfs pull` fails or is blocked

### Prerequisites

- Git
- Python 3.6 or newer for the download scripts
- `pip`

### Setup

1. Clone the repository normally:

```bash
git clone <REPO_URL>
cd CAAARal-Reefs
```

2. Install the download dependency:

```bash
python3 -m pip install gdown
```

3. Run the recommended data download script:

```bash
python3 scripts/download_necessary_data.py
```

This script downloads the required Google Drive folders into `data/`.

## Available download scripts

### `scripts/download_necessary_data.py`

Recommended for judges. This downloads the required shared Google Drive resources into the repository's `data/` directory.

Run:

```bash
python3 scripts/download_necessary_data.py
```

### `scripts/download_data.py`

This is a more general single-link downloader. In its current configuration it downloads a shared Google Drive folder into the repository root.

Run:

```bash
python3 scripts/download_data.py
```

Use this only if you specifically want the full shared bundle referenced in that script.

## Quick verification

After setup, confirm that the repository contains populated data directories such as:

- `data/raw/`
- `data/processed/`
- `data/analyzed/`

If those folders contain actual datasets rather than Git LFS pointer text, the project data has been retrieved correctly.

## Reproducing the project

After the data is present, judges can inspect the project through:

- Jupyter notebooks in `notebooks/cleaning/` and `notebooks/analysis/`
- MATLAB code in `MATLAB/`
- generated artifacts in `outputs/` and `MATLAB/Figures/`

## Judge workflow summary

For MTFC judges, the simplest workflow is:

```bash
git clone <REPO_URL>
cd CAAARal-Reefs
python3 -m pip install gdown
python3 scripts/download_necessary_data.py
```

If Git LFS is installed and working, you may use `git lfs pull` instead of the Google Drive step.
