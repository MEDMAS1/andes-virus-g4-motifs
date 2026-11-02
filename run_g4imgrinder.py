#!/usr/bin/env python
"""Python-2-safe launcher wrapper for run_g4imgrinder_main.py."""

import os
import sys
try:
    from shutil import which
except ImportError:
    from distutils.spawn import find_executable as which


def main():
    if sys.version_info[0] < 3:
        python3 = which("python3")
        if python3 is None:
            python3 = which("python3.11") or which("python3.10") or which("python3.9")
        if python3:
            os.execvp(python3, [python3] + sys.argv)
        raise SystemExit("This script requires Python 3. Run with python3 run_g4imgrinder.py")

    from run_g4imgrinder_main import main as run_main
    run_main()


if __name__ == "__main__":
    main()
