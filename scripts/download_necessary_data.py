#!/usr/bin/env python3
"""
download_necessary_data.py

Download all required Google Drive resources listed in DRIVE_URLS.

Setup:
    pip install gdown

Run:
    python3 scripts/download_necessary_data.py
"""

import sys
from pathlib import Path
from typing import Dict, List, Tuple, Union

from drive_download_utils import build_download_jobs, download_url

# =========================
# EDIT THESE VALUES
# =========================

# Supports:
# 1) plain URL string (uses DEFAULT_TARGET_DIR), or
# 2) dict with "url" and optional "target_dir".
DRIVE_URLS: List[Union[str, Dict[str, str]]] = [
    "https://drive.google.com/drive/folders/1lXqczOkri_ISOMOnxuQT8vlHoBQUwSot?usp=sharing",
    "https://drive.google.com/drive/folders/1q_ZZ-S-bb2_M5zdjImzXOiC1mKsdfRRa?usp=sharing",
    # {"url": "https://drive.google.com/file/d/<file_id>/view?usp=sharing", "target_dir": "../data/raw"},
]

# Used when an item in DRIVE_URLS is a plain URL string.
DEFAULT_TARGET_DIR = "../data"

# If True, stop at first failure. If False, continue and report all failures.
FAIL_FAST = False


def main() -> int:
    if not DRIVE_URLS:
        print(
            "Error: DRIVE_URLS is empty. Add at least one Google Drive URL.",
            file=sys.stderr,
        )
        return 1

    scripts_dir = Path(__file__).resolve().parent
    failures: List[Tuple[str, str]] = []

    jobs = build_download_jobs(DRIVE_URLS, scripts_dir, DEFAULT_TARGET_DIR)
    for index, (url, target_dir) in enumerate(jobs, start=1):
        print(f"[{index}/{len(jobs)}] Processing: {url}")
        try:
            download_url(url, target_dir)
        except Exception as exc:
            failures.append((url, str(exc)))
            print(f"Failed: {exc}", file=sys.stderr)
            if FAIL_FAST:
                break

    if failures:
        print("\nSome downloads failed:", file=sys.stderr)
        for url, error in failures:
            print(f"- {url}\n  {error}", file=sys.stderr)
        return 1

    print("All downloads completed.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
