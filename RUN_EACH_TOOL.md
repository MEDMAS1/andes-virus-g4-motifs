# Running Each Tool Separately

Put FASTA files in `data/`. Each launcher below reads one file or all FASTA files from `data/` and writes results under `results/<tool>/<input-name>/`.

Run from the project folder with Python 3:

```powershell
py -3 run_qgrs.py
py -3 run_qparse.py
py -3 run_g4predict.py
py -3 run_g4imgrinder.py
py -3 run_g4promfinder.py
py -3 run_g4catchall.py
py -3 run_g4rna_screener.py
py -3 run_g4hunter.py
py -3 run_g4mismatch.py
py -3 run_rg4detector.py
# Optional: Quadron requires library/Quadron/Quadron.lib
py -3 run_quadron.py
```

To run only one input file:

```powershell
py -3 run_qgrs.py --input "data/L reverse.fasta"
```

## Compare common predictions across tools

After running the tools, create combined reports for common sequences, common
start positions, overlapping predictions, and sequences included inside longer
predictions:

```powershell
py -3 compare_g4_outputs.py
```

The reports are written to `results/common/`.

The Python 2 tools use the bundled `library/g4rna_env/bin/python` through WSL on Windows. The TensorFlow tools use `python3` by default; set `PYTHON3` if your TensorFlow environment has a different command.

`G4mismatch` and `rG4detector` need DNA alphabet input, so their launchers create a normalized copy in the result folder with `U` changed to `T`. Your original file in `data/` is not changed.
