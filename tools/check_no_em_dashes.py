#!/usr/bin/env python3
"""
Pre-commit / verification script to ensure no em dashes or en dashes
exist in any comment, docstring, or documentation file in the repository.
"""

import os
import sys

EXCLUDED_DIRS = {'.git', '.zig-cache', 'dist-newstyle', '.cabal-sandbox'}
TEXT_EXTS = {'.hs', '.lhs', '.c', '.h', '.py', '.md', '.html', '.css', '.js', '.cabal', '.json', '.yaml', '.yml', '.txt', 'Makefile', 'LICENSE'}

def check_file(filepath):
    errors = []
    with open(filepath, 'r', encoding='utf-8', errors='replace') as f:
        for idx, line in enumerate(f, start=1):
            # Check for unicode em dash (U+2014) or en dash (U+2013)
            if '\u2014' in line:
                errors.append((idx, "Contains Unicode em dash (U+2014)"))
            if '\u2013' in line:
                errors.append((idx, "Contains Unicode en dash (U+2013)"))
            # In markdown/docs/prose, check for ' -- ' used as an em dash
            if filepath.endswith(('.md', '.txt', '.html')):
                if ' -- ' in line:
                    errors.append((idx, "Contains ' -- ' used as punctuation em dash"))
    return errors

def main():
    root = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
    total_errors = 0
    checked_files = 0

    for dirpath, dirnames, filenames in os.walk(root):
        dirnames[:] = [d for d in dirnames if d not in EXCLUDED_DIRS]
        for fname in filenames:
            ext = os.path.splitext(fname)[1]
            if ext in TEXT_EXTS or fname in ('Makefile', 'LICENSE'):
                fpath = os.path.join(dirpath, fname)
                checked_files += 1
                errs = check_file(fpath)
                if errs:
                    for line_no, msg in errs:
                        print(f"ERROR: {os.path.relpath(fpath, root)}:{line_no}: {msg}")
                        total_errors += 1

    if total_errors > 0:
        print(f"\nFAILED: Found {total_errors} dash violations across {checked_files} files.")
        sys.exit(1)
    else:
        print(f"PASSED: Verified {checked_files} files. Zero em dashes found.")
        sys.exit(0)

if __name__ == '__main__':
    main()
