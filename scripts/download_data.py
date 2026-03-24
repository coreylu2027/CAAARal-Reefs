#!/usr/bin/env python3
"""
download_data.py

Downloads a Google Drive file or folder into a target directory in your repo.
Good for research projects that moved large data off Git LFS and into Google Drive.

Setup:
    pip install gdown

Run:
    python download_data.py
"""

from __future__ import annotations

import re
import sys
from pathlib import Path

import gdown

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


def extract_drive_id(url: str) -> str | None:
    patterns = [
        r"/folders/([a-zA-Z0-9_-]+)",
        r"/file/d/([a-zA-Z0-9_-]+)",
        r"[?&]id=([a-zA-Z0-9_-]+)",
    ]
    for pattern in patterns:
        match = re.search(pattern, url)
        if match:
            return match.group(1)
    return None


def is_folder_url(url: str) -> bool:
    return "/drive/folders/" in url or "folders/" in url


def is_file_url(url: str) -> bool:
    return "/file/d/" in url or "id=" in url


def target_has_files(path: Path) -> bool:
    return path.exists() and any(path.iterdir())


def download_drive_folder(url: str, output_dir: Path) -> None:
    folder_id = extract_drive_id(url)
    if not folder_id:
        raise ValueError("Could not extract folder ID from Google Drive URL.")

    print(f"Downloading Google Drive folder into: {output_dir}")
    gdown.download_folder(
        id=folder_id,
        output=str(output_dir),
        quiet=False,
        use_cookies=False,
        remaining_ok=True,
    )


def download_drive_file(url: str, output_dir: Path) -> None:
    file_id = extract_drive_id(url)
    if not file_id:
        raise ValueError("Could not extract file ID from Google Drive URL.")

    output_dir.mkdir(parents=True, exist_ok=True)
    print(f"Downloading Google Drive file into: {output_dir}")
    gdown.download(
        id=file_id,
        output=str(output_dir),
        quiet=False,
        fuzzy=True,
        use_cookies=False,
    )


def main() -> int:
    if DRIVE_URL == "PASTE_YOUR_GOOGLE_DRIVE_LINK_HERE":
        print("Error: Please edit DRIVE_URL in the script first.", file=sys.stderr)
        return 1

    repo_dir = Path(__file__).resolve().parent
    target_dir = (repo_dir / TARGET_DIR).resolve()
    target_dir.mkdir(parents=True, exist_ok=True)

    if SKIP_IF_TARGET_NOT_EMPTY and target_has_files(target_dir):
        print(f"Skipping download: target directory already has files: {target_dir}")
        return 0

    try:
        if is_folder_url(DRIVE_URL):
            download_drive_folder(DRIVE_URL, target_dir)
        elif is_file_url(DRIVE_URL):
            download_drive_file(DRIVE_URL, target_dir)
        else:
            raise ValueError(
                "Unsupported Google Drive URL. Use a shared file or folder link."
            )

        print("Download complete.")
        return 0

    except Exception as exc:
        print(f"Download failed: {exc}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())