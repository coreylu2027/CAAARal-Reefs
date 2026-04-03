#!/usr/bin/env python3
"""
download_data.py

Downloads a Google Drive file or folder into a target directory in your repo.
Good for research projects that moved large data off Git LFS and into Google Drive.

Setup:
    pip install gdown

Run:
    python3 scripts/download_data.py
"""

import sys
from pathlib import Path

from drive_download_utils import download_url, target_has_files

# =========================
# EDIT THESE VALUES
# =========================

# Paste your Google Drive share link here.
# Examples:
# DRIVE_URL = "https://drive.google.com/drive/folders/1AbCdEfGhIjKlMnOpQrStUvWxYz"
# DRIVE_URL = "https://drive.google.com/file/d/1AbCdEfGhIjKlMnOpQrStUvWxYz/view?usp=sharing"
DRIVE_URL = "https://drive.google.com/drive/folders/1XeJsYMXnk2fmRX29OH0kY12XAZyOmZ0Z?usp=share_link"

# Folder inside your repo where the data should go.
# Example: "data", "datasets/raw", "../project_data"
TARGET_DIR = "../"

# Set to True if you want to skip downloading when the target directory
# already contains files.
SKIP_IF_TARGET_NOT_EMPTY = False


def main() -> int:
    if DRIVE_URL == "PASTE_YOUR_GOOGLE_DRIVE_LINK_HERE":
        print("Error: Please edit DRIVE_URL in the script first.", file=sys.stderr)
        return 1

    scripts_dir = Path(__file__).resolve().parent
    target_dir = (scripts_dir / TARGET_DIR).resolve()
    target_dir.mkdir(parents=True, exist_ok=True)

    if SKIP_IF_TARGET_NOT_EMPTY and target_has_files(target_dir):
        print(f"Skipping download: target directory already has files: {target_dir}")
        return 0

    try:
        download_url(DRIVE_URL, target_dir)
        print("Download complete.")
        return 0

    except Exception as exc:
        print(f"Download failed: {exc}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
