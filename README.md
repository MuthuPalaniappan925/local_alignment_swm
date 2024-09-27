# **Local Sequence Alignment using SWM Algo**
- Supports DNA sequence alignment (A, G, C, T) or any custom sequences.
- Flexible scoring system: 
  - Default scoring for match, mismatch, and gap penalties.
  - Custom scoring via a CSV file for advanced scoring.
- Outputs the alignment and matrices (result matrix, traceback matrix).
- Traceback visualization showing optimal alignment.
- Argument parsing for easy command-line interface (CLI) usage.

## **Installation**

### **1. Clone the Repository**
```bash
git clone https://github.com/MuthuPalaniappan925/local_alignment_swm.git
cd local_alignment_swm
```

### **2. Install Dependencies**
This project uses `numpy` and `pandas` for matrix calculations and custom CSV handling. Install the required Python packages via pip:
```bash
pip install numpy pandas
```

## **Usage**

You can run the alignment via the command line. Below are some options to customize the execution.

### **Command Line Arguments**

- **`seq_a_file`** (required): Path to the file containing sequence A.
- **`seq_b_file`** (required): Path to the file containing sequence B.
- **`output_file`** (required): Path to the output file where alignment results are saved.
- **`--match`** (optional): Score for a match (default: `2`).
- **`--mismatch`** (optional): Penalty for a mismatch (default: `-1`).
- **`--gap`** (optional): Penalty for a gap (default: `-2`).
- **`--load_csv`** (optional): Load a custom scoring matrix from a CSV file (`scores.csv`).

### **Example**
```
python3 cpu_swm_baseline.py seq_1.txt seq_2.txt alg_out.txt --match 2 --mismatch -1 --gap -2
```

## **Custom Scoring CSV**
To use custom scoring, you need a CSV file (`scores.csv`) with the following structure:

```
  A G C T
A 2 -1 -1 -1
G -1 2 -1 -1
C -1 -1 2 -1
T -1 -1 -1 2
```

The first row and first column represent the possible sequence characters, while the cells contain the respective match/mismatch scores.

## **Classes Overview**

### **1. `ScoringSystem`**
Responsible for determining the score between any two characters in the sequence, based on the match, mismatch, or gap penalties. It can load a custom scoring matrix from a CSV file.

- **Methods**:
  - `__init__(self, match: int=2, mismatch: int=-1, gap: int=-2)`: Initializes with default match, mismatch, and gap scores.
  - `load_csv(self, filename: str)`: Loads a custom scoring matrix from a CSV file.
  - `score(self, a: str, b: str)`: Returns the score between two characters based on the scoring matrix.
  
### **2. `SequencesAnalyzer`**
Performs local sequence alignment using the Smith-Waterman algorithm.

- **Methods**:
  - `__init__(self, seq_a_file: str, seq_b_file: str, match: int=2, mismatch: int=-1, gap: int=-2, load_csv: bool=False)`: Initializes the analyzer with sequences and scoring options.
  - `local_alignment(self, output_filename: str)`: Runs the alignment and saves the results in the output file.
  - `smith_waterman_algorithm(self)`: Implements the Smith-Waterman algorithm.
  - `_traceback(self, result_matrix, traceback_matrix, start_pos: Tuple[int, int])`: Traces back through the matrices to generate the final aligned sequences.

## **Output**
The output consists of:
- The input sequences.
- The alignment score.
- Aligned sequences with match/mismatch indicators (`|` for matches, space for mismatches).
- The result and traceback matrices.
- Execution time.

Example output (saved in the output file):
```
Input Sequence A: AGCT
Input Sequence B: GACT

Alignment Score: 3

Alignment:
AGCT
 | |
GACT

Result Matrix:
 [[0 0 0 0]
  [0 0 2 0]
  [0 2 1 1]
  [0 0 1 3]
  [0 0 0 0]]

Traceback Matrix:
 [['' '' '' '']
  ['' '↖' '←' '↑']
  ['' '↑' '↖' '←']
  ['' '↑' '↑' '↖']
  ['' '←' '←' '']]
```
