#!/usr/bin/env python3
"""
download_necessary_data.py

Download all required Google Drive resources listed in DRIVE_URLS.

Setup:
    pip install gdown

Run:
    python scripts/download_necessary_data.py
"""

import re
import sys
from pathlib import Path
from typing import Dict, Iterable, List, Tuple, Union

import gdown

# =========================
# EDIT THESE VALUES
# =========================

# Supports:
# 1) plain URL string (uses DEFAULT_TARGET_DIR), or
# 2) dict with "url" and optional "target_dir".
DRIVE_URLS: List[Union[str, Dict[str, str]]] = [
    "https://drive.google.com/drive/folders/1lXqczOkri_ISOMOnxuQT8vlHoBQUwSot?usp=sharing",
    "https://drive.google.com/drive/folders/1q_ZZ-S-bb2_M5zdjImzXOiC1mKsdfRRa?usp=sharing"
    # {"url": "https://drive.google.com/file/d/<file_id>/view?usp=sharing", "target_dir": "../data/raw"},
]

# Used when an item in DRIVE_URLS is a plain URL string.
DEFAULT_TARGET_DIR = "../data"

# If True, stop at first failure. If False, continue and report all failures.
FAIL_FAST = False


def extract_drive_id(url: str):
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
    return "/drive/folders/" in url or "/folders/" in url


def is_file_url(url: str) -> bool:
    return "/file/d/" in url or "id=" in url


def download_drive_folder(url: str, output_dir: Path) -> None:
    folder_id = extract_drive_id(url)
    if not folder_id:
        raise ValueError("Could not extract folder ID from Google Drive URL.")

    output_dir.mkdir(parents=True, exist_ok=True)
    print(f"Downloading folder to: {output_dir}")
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
    print(f"Downloading file to: {output_dir}")
    gdown.download(
        id=file_id,
        output=str(output_dir),
        quiet=False,
        fuzzy=True,
        use_cookies=False,
    )


def iter_download_jobs(
    drive_urls: Iterable[Union[str, Dict[str, str]]],
    scripts_dir: Path,
) -> Iterable[Tuple[str, Path]]:
    default_target = (scripts_dir / DEFAULT_TARGET_DIR).resolve()
    for item in drive_urls:
        if isinstance(item, str):
            yield item, default_target
            continue

        if "url" not in item:
            raise ValueError("Each dict item in DRIVE_URLS must include a 'url' key.")

        target = item.get("target_dir", DEFAULT_TARGET_DIR)
        yield item["url"], (scripts_dir / target).resolve()


def main() -> int:
    if not DRIVE_URLS:
        print("Error: DRIVE_URLS is empty. Add at least one Google Drive URL.", file=sys.stderr)
        return 1

    scripts_dir = Path(__file__).resolve().parent
    failures: List[Tuple[str, str]] = []

    jobs = list(iter_download_jobs(DRIVE_URLS, scripts_dir))
    for index, (url, target_dir) in enumerate(jobs, start=1):
        print(f"[{index}/{len(jobs)}] Processing: {url}")
        try:
            if is_folder_url(url):
                download_drive_folder(url, target_dir)
            elif is_file_url(url):
                download_drive_file(url, target_dir)
            else:
                raise ValueError("Unsupported Google Drive URL. Use shared file/folder links.")
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
