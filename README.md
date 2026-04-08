# CAAARal-Reefs

This repository contains the code, notebooks, data files, and figures for the 2025-26 Modeling the Future Challenge project titled "A Risk Analysis of the Effects of Coral Reef Degradation in Tropical Climates" by Team #23914 (CAAARal Reefs).

## Project Description
This project presents a comprehensive quantitative risk analysis of coral reef degradation, highlighting its severe ecological, economic, and societal impacts. It shows that reefs, which support over half a billion people, are rapidly declining due to rising ocean temperatures, acidification, and increased storm activity, leading to losses in fisheries, tourism, coastal protection, and biodiversity. Using over 1.6 million data points and multiple integrated models, including stochastic fishery dynamics, SARIMA-based coral cover forecasts, economic loss models, and a pharmaceutical value framework, the study projects substantial future damages such as billions in annual economic losses and high probabilities of fishery collapse by 2100. It also identifies irreversible losses in pharmaceutical potential due to species extinction as a uniquely critical risk. The project concludes by proposing a multi-faceted policy framework that combines restoration, behavioral incentives, insurance mechanisms, and biodiversity preservation to mitigate these risks and protect both ecosystems and global economies before critical thresholds are reached.

# Information for Judges
For Modeling the Future Challenge (MTFC) judges, the most important thing to know is that this project utilizes Git LFS (large file storage), but it also includes Google Drive download scripts as a fallback, in case Git LFS is unavailable or incomplete on your machine. Python dependencies for the notebooks and download scripts are listed in `requirements.txt` and can be installed accordingly.

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

4. Install the Python dependencies used by the notebooks and download helpers:

```bash
python3 -m pip install -r requirements.txt
```

## Option 2: Clone normally and use the Google Drive download scripts

Use this if:

- `git lfs` is not installed
- Git LFS files appear as small pointer text files instead of real datasets
- `git lfs pull` fails or is blocked

### Prerequisites

- Git
- Python 3
- `pip`

### Setup

1. Clone the repository normally:

```bash
git clone <REPO_URL>
cd CAAARal-Reefs
```

2. Install the Python dependencies:

```bash
python3 -m pip install -r requirements.txt
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

After the data is present and `requirements.txt` has been installed, judges can inspect the project through:

- Jupyter notebooks in `notebooks/cleaning/` and `notebooks/analysis/`
- MATLAB code in `MATLAB/`
- generated artifacts in `outputs/` and `MATLAB/Figures/`

## MATLAB workflow (current)

Core MATLAB modeling scripts live in `MATLAB/CoralReefDegradation/`.

1. Run `benthicCoverARIMA.m`
2. Run `WaveEnergy_Tourism.m`
3. Run `AtharvSensitivity.m`

These scripts are wired to share CSV inputs and outputs directly from `MATLAB/CoralReefDegradation/`.

Key generated files include:

- `MATLAB/CoralReefDegradation/benthic_cover_forecast.csv`
- `MATLAB/CoralReefDegradation/benthic_cover_observed_annual.csv`
- `MATLAB/CoralReefDegradation/reef_economic_loss_forecast.csv`
- `MATLAB/CoralReefDegradation/reef_economic_loss_observed.csv`
- `MATLAB/CoralReefDegradation/reef_economic_mitigation.csv`
- `MATLAB/CoralReefDegradation/reef_economic_npv.csv`

Generated MATLAB figures are saved to `MATLAB/Figures/`.

## Judge workflow summary

For MTFC judges, the simplest workflow is:

```bash
git clone <REPO_URL>
cd CAAARal-Reefs
python3 -m pip install -r requirements.txt
python3 scripts/download_necessary_data.py
```

If Git LFS is installed and working, you may use `git lfs pull` instead of the Google Drive step.

Best,
CAAARal Reefs
