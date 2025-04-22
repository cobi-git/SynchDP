# SynchDP: A Correlation Based Sequence Alignment Algorithm for Synchronizing Longitudinal Clinical Data

Time-series data in the real world—such as biological signals or power consumption logs—often contain asynchronous patterns due to delays or sampling differences.  
**SynchDP** is a Python-based local alignment algorithm designed to synchronize such time-series data using **dynamic programming based on Pearson correlation**.

This algorithm was originally developed for aligning disease progression patterns in clinical data, but its flexible structure allows application to any domain involving time-series analysis.

## 📌 Features

- Local alignment of time-series data using Pearson correlation
- Handles both synchronous and asynchronous sequence comparison
- Alignment score matrix and alignment path output
- Visualization utilities for alignment result
- Network construction and clustering tools (e.g., Louvain clustering)
- Tutorial codes for full pipeline usage

![workflow](./images/SynchDP_Method.pdf)

---

## 🔧 Installation

We recommend using a virtual environment with the following packages:

```bash
pip install -r requirements.txt
```

Main dependencies:
```python
numpy
scipy
pandas
matplotlib
networkx
python-louvain
```

---


## 🧪 Tutorial Usage

You can follow the step-by-step usage examples in the `tutorial/` folder.

### 📁 Input Data

Example input files are provided in `tutorial/input/`:

- `omics_info.csv` — Metadata for sequences
- `query_set.pkl` — Query time-series data
- `reference_seq_set.pkl` — Reference time-series data
- `date_info.pkl` — Timestamps for each data point

Output files will be saved in `tutorial/output/`.

You can generate the input files with:

```bash
python tutorial/generate_example_inputs.py
```

---

### 🔹 Reference Mode

Align query sequences to reference sequences based on time alignment.

```bash
python src/synchdp.py \
  -m reference \
  -q tutorial/input/querry_set.pkl \
  -t tutorial/input/reference_seq_set.pkl \
  --date_info tutorial/input/date_info.pkl \
  --omics_info tutorial/input/omics_info.csv \
  -p 0.1 -c 0.65 -w 6 -e 3 \
  -o tutorial/out
```

**Output:**
- Aligned sequences with assigned reference labels
- Synchronization results (PDF)

---

### 🔹 De Novo Mode

Align all sequences to each other without a reference.

```bash
python src/synchdp.py \
  -m denovo \
  -q tutorial/input/querry_set.pkl \
  --date_info tutorial/input/date_info.pkl \
  --omics_info tutorial/input/omics_info.csv \
  -r 95 -p 0.01 -c 0.7 -w 6 -e 3.5 \
  -o tutorial/outdenovo
```

**Output:**
- Pairwise alignment matrix
- Network graph (PDF)
- Reference sequence by hub node
- Synchronization results (PDF)

---

### 🧩 Notes

- `--mode` can be either `reference` or `denovo`
- `-o` (or `--output_dir`) specifies the output directory
- All output files will be saved inside the directory specified by `-o`

---

### 🔧 Command Line Options

| Option           | Description                                                    |
|------------------|----------------------------------------------------------------|
| `-m`             | Mode of alignment: `reference` or `denovo`                     |
| `-q`             | Query sequence `.pkl` file                                     |
| `-t`             | Reference sequence `.pkl` file (required in reference mode)    |
| `--date_info`    | Date information `.pkl` file                                   |
| `--omics_info`   | Metadata `.csv` file containing sample labels                  |
| `-o`             | Output directory to save results                               |
| `-p`             | Penalty score used in alignment scoring (e.g., 0.1)            |
| `-c`             | Pearson correlation threshold for alignment (e.g., 0.65)       |
| `-w`             | Alignment window size (e.g., 6)                                |
| `-e`             | Euclidean distance threshold for sequence filtering (e.g., 3.0)|
| `-r`             | Retain top-r percentile edges for network (only for denovo)    |

---

You can also check available options by running:

```bash
python src/synchdp.py --help
```

---

## 📁 Project Structure

```
SynchDP/
│
├── src/                        # Core algorithm modules
│   ├── argparser.py            # Argument parser for CLI usage
│   ├── synchdp.py              # Script for executing alignment and synchronization
│   ├── synchdpmod.py           # Core SynchDP algorithm implementation
│   └── synchronize.py          # Integrated module combining de novo and reference-based modes
│
│
├── tutorial/                   # Tutorial and usage examples
│   ├── input/                  # Input data files
│   │   ├── date_info.pkl
│   │   ├── omics_info.csv
│   │   ├── querry_set.pkl
│   │   └── reference_seq_set.pkl
│   │
│   ├── output/                 # Output results
│   │   ├── SynchDP_Network(De novo).pdf
│   │   ├── SynchDP_results(De_novo).pdf
│   │   └── SynchDP_results(reference_guide).pdf
│   │
│   └── generate_example_inputs.py   # Script to generate example input files
│
└── README.md                  # Main project documentation
```


## 📬 Contact

For questions, issues, or contributions, please open an issue or reach out via GitHub.
