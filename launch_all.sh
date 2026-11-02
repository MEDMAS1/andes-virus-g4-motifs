#!/usr/bin/env bash
set -euo pipefail

cd "$(dirname "$0")"

python3 run_qgrs.py --score 0
python3 run_qparse.py
python3 run_g4predict.py
python3 run_g4imgrinder.py
python3 run_g4promfinder.py
python3 run_g4catchall.py
python3 run_g4hunter.py
python3 run_g4rna_screener.py
python3 run_pqsfinder.py
python3 run_g4mismatch.py
python3 run_rg4detector.py
#Optional: Quadron requires library/Quadron/Quadron.lib and Rscript
#python3 run_quadron.py
