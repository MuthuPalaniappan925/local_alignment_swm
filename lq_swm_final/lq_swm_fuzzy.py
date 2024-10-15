##Importing Pacakges
import numpy as np
import os
import time
import argparse
from typing import Tuple, Dict, Any

## Define a Scoring System for the fuzzy-based Smith-Waterman alignment
class ScoringSystem:
    def __init__(self, match: int = 2, mismatch: int = -1, gap: int = -2, quality_weights: Dict[str, int] = None) -> None:
        """
        Initialize the scoring system with match, mismatch, gap, and quality weights.

        :param match: Score for a match.
        :param mismatch: Penalty for a mismatch.
        :param gap: Penalty for a gap.
        :param quality_weights: A dictionary mapping quality characters to their weights.
        """
        self.match = match
        self.mismatch = mismatch
        self.gap = gap
        quality_weights = quality_weights = {
    '!': 0.0, '"': 0.01, '#': 0.02, '$': 0.03, '%': 0.04,
    '&': 0.05, "'": 0.06, '(': 0.07, ')': 0.08, '*': 0.09,
    '+': 0.1, ',': 0.11, '-': 0.12, '.': 0.13, '/': 0.14,
    '0': 0.15, '1': 0.16, '2': 0.17, '3': 0.18, '4': 0.19,
    '5': 0.2, '6': 0.21, '7': 0.22, '8': 0.23, '9': 0.24,
    ':': 0.25, ';': 0.26, '<': 0.27, '=': 0.28, '>': 0.29,
    '?': 0.3, '@': 0.31, 'A': 0.32, 'B': 0.33, 'C': 0.34,
    'D': 0.35, 'E': 0.36, 'F': 0.37, 'G': 0.38, 'H': 0.39,
    'I': 0.4
}
        self.quality_weights = quality_weights

    def _fuzzy_similarity(self, a: str, b: str, quality_a: str) -> float:
        """
        Calculate a fuzzy similarity score based on quality weights.

        :param a: Character from sequence A.
        :param b: Character from sequence B.
        :param quality_a: Quality score character for sequence A.
        :return: Fuzzy similarity score.
        """
        weight_a = self.quality_weights.get(quality_a, 1)
        return weight_a

    def _default_scoring(self, a: str, b: str, quality_a: str) -> float:
        """
        Compute the default scoring based on match, mismatch, and gap penalties.

        :param a: Character from sequence A.
        :param b: Character from sequence B.
        :param quality_a: Quality score character for sequence A.
        :return: Computed score.
        """
        if a == '-' or b == '-':
            return self.gap
        fuzzy_score = self._fuzzy_similarity(a, b, quality_a)
        if a != b:
            score = self.mismatch * fuzzy_score
        else:
            score = self.match * fuzzy_score
        return score

    def score(self, a: str, b: str, quality_a: str) -> float:
        """
        Return the score between two characters considering quality weights.

        :param a: Character from sequence A.
        :param b: Character from sequence B.
        :param quality_a: Quality score character for sequence A.
        :return: Calculated score.
        """
        assert isinstance(a, str) and isinstance(b, str)
        assert len(a) == 1 and len(b) == 1
        return self._default_scoring(a, b, quality_a)

## Define the sequence alignment class
class SequencesAnalyzer:
    traceback_symbols = {
        0: '↖',  # Diagonal
        1: '↑',  # Up
        2: '←'   # Left
    }

    def __init__(self, seq_a: str, seq_b: str, quality_a: str, match: int = 2, mismatch: int = -1, gap: int = -2, quality_weights: Dict[str, int] = None) -> None:
        """
        Initialize the SequencesAnalyzer with sequences and scoring parameters.

        :param seq_a: First sequence.
        :param seq_b: Second sequence.
        :param quality_a: Quality string corresponding to seq_a.
        :param match: Match score.
        :param mismatch: Mismatch penalty.
        :param gap: Gap penalty.
        :param quality_weights: Dictionary of quality weights.
        """
        self.seq_a = seq_a
        self.seq_b = seq_b
        self.quality_a = quality_a
        self.scoring_sys = ScoringSystem(match, mismatch, gap, quality_weights)

    def smith_waterman_algorithm(self) -> Dict[str, Any]:
        """
        Execute the Smith-Waterman algorithm for local sequence alignment.

        :return: A dictionary containing result matrix, traceback matrix, score, and score position.
        """
        rows, cols = len(self.seq_a) + 1, len(self.seq_b) + 1
        H = np.zeros(shape=(rows, cols), dtype=float)
        traceback = np.zeros(shape=(rows, cols), dtype=np.dtype('U5'))

        max_score = 0
        max_pos = (0, 0)

        for row in range(1, rows):
            for col in range(1, cols):
                a = self.seq_a[row - 1]
                b = self.seq_b[col - 1]
                qa = self.quality_a[row - 1] if row - 1 < len(self.quality_a) else '-'

                score_diag = H[row - 1, col - 1] + self.scoring_sys.score(a, b, qa)
                score_up = H[row - 1, col] + self.scoring_sys.gap
                score_left = H[row, col - 1] + self.scoring_sys.gap

                H[row, col] = max(0, round(score_diag, 2), round(score_up, 2), round(score_left, 2))
                
                if H[row, col] == round(score_diag, 2):
                    traceback[row, col] = self.traceback_symbols[0]
                elif H[row, col] == round(score_up, 2):
                    traceback[row, col] = self.traceback_symbols[1]
                elif H[row, col] == round(score_left, 2):
                    traceback[row, col] = self.traceback_symbols[2]

                if H[row, col] > max_score:
                    max_score = H[row, col]
                    max_pos = (row, col)

        return {
            'result_matrix': H,
            'traceback_matrix': traceback,
            'score': max_score,
            'score_pos': max_pos
        }

    def _traceback(self, result_matrix, traceback_matrix, start_pos: Tuple[int, int]) -> Tuple[str, str, str]:
        """
        Traceback through the alignment matrices to build the final aligned sequences.

        :param result_matrix: Matrix of alignment scores.
        :param traceback_matrix: Matrix of traceback symbols.
        :param start_pos: Position to start the traceback.
        :return: A tuple of aligned sequence A, sequence B, and match symbols.
        """
        seq_a_aligned = ''
        seq_b_aligned = ''
        matches = ''

        row, col = start_pos
        while result_matrix[row, col] > 0:
            symbol = traceback_matrix[row, col]

            if symbol == '↖':
                seq_a_aligned += self.seq_a[row - 1]
                seq_b_aligned += self.seq_b[col - 1]
                matches += '|' if self.seq_a[row - 1] == self.seq_b[col - 1] else ' '
                row -= 1
                col -= 1
            elif symbol == '↑':
                seq_a_aligned += self.seq_a[row - 1]
                seq_b_aligned += '-'
                matches += ' '
                row -= 1
            elif symbol == '←':
                seq_a_aligned += '-'
                seq_b_aligned += self.seq_b[col - 1]
                matches += ' '
                col -= 1
            else:
                break

        return seq_a_aligned[::-1], seq_b_aligned[::-1], matches[::-1]

def generate_random_sequence(base_seq: str, n: int) -> str:
    """
    Generate a new sequence with the first `n` chars of `base_seq` same, rest randomized.
    
    :param base_seq: Base sequence to derive from.
    :param n: Number of initial bases to keep identical.
    :return: Generated sequence.
    """
    from random import choice
    bases = ['A', 'T', 'C', 'G']
    length = len(base_seq)
    return base_seq[:n] + ''.join(choice(bases) for _ in range(length - n))

def parse_fastq_file(filename: str) -> Tuple[list, list]:
    """
    Parse a FASTQ file and return lists of read sequences and their quality scores.
    
    :param filename: Path to the FASTQ file.
    :return: A tuple (sequences, quality_scores) where:
             sequences is a list of read sequences.
             quality_scores is a list of quality score strings.
    """
    sequences = []
    quality_scores = []
    with open(filename, 'r', encoding='utf-8') as file:
        lines = file.readlines()
        for i in range(0, len(lines), 4):
            if lines[i].startswith('@') and (i + 1 < len(lines)) and (i + 3 < len(lines)):
                read_sequence = lines[i + 1].strip()
                quality_value = lines[i + 3].strip()
                sequences.append(read_sequence)
                quality_scores.append(quality_value)
    return sequences, quality_scores

def run_alignment_from_fastq(fastq_filename: str, n: int, k: int, quality_weights: Dict[str, int], output_filename: str, iterations: int = 1) -> None:
    """
    Run alignment for sequences parsed from a FASTQ file and save the results to a file.
    
    :param fastq_filename: Path to the FASTQ file.
    :param n: Number of initial bases to retain while generating the randomized sequence.
    :param k: Length of the sequence to generate.
    :param quality_weights: Dictionary of quality weights.
    :param output_filename: Output file path.
    :param iterations: Number of iterations to run for each read.
    """
    reads, qualities = parse_fastq_file(fastq_filename)
    assert len(reads) > 0, f"ERROR: No read sequences found in {fastq_filename}"

    with open(output_filename, 'w', encoding='utf-8') as f:
        for idx, (seq_a, quality_a) in enumerate(zip(reads, qualities), start=1):
            for i in range(k):
                seq_b = generate_random_sequence(seq_a, n)
                analyzer = SequencesAnalyzer(seq_a, seq_b, quality_a, quality_weights=quality_weights)
                result = analyzer.smith_waterman_algorithm()

                alignment_a, alignment_b, matches = analyzer._traceback(
                    result['result_matrix'], 
                    result['traceback_matrix'], 
                    result['score_pos']
                )

                f.write(f"Read {idx} - Iteration {i+1}:\n")
                f.write(f"Input Sequence A: {seq_a}\n")
                f.write(f"Generated Sequence B: {seq_b}\n")
                f.write(f"Alignment Score: {result['score']}\n\n")
                f.write("Alignment:\n")
                f.write(f"{alignment_a}\n")
                f.write(f"{matches}\n")
                f.write(f"{alignment_b}\n\n")
                f.write(f"Result Matrix:\n{result['result_matrix']}\n\n")
                f.write(f"Traceback Matrix:\n{result['traceback_matrix']}\n\n")
                f.write(f"Score Position: {result['score_pos']}\n")
                f.write("-" * 80 + "\n")

                print(f"[INFO] Completed Read {idx}, Iteration {i+1}")

    print(f"[INFO] Results written to {output_filename}")
    

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Sequence Alignment Tool.')
    parser.add_argument('--fastq', type=str, required=True, help='Path to the FASTQ file')
    parser.add_argument('--k', type=int, default=3, help='Number of iterations (default: 3)')
    parser.add_argument('--n', type=int, default=10, help='Number of bases to retain in the generated sequence (default: 10)')
    parser.add_argument('--output', type=str, required=True, help='Output file path')
    parser.add_argument('--quality_weights', type=str, help='Optional custom quality weights (e.g., "A:3,B:2,C:1")')

    args = parser.parse_args()

    quality_weights = None
    if args.quality_weights:
        quality_weights = {kv.split(':')[0]: int(kv.split(':')[1]) for kv in args.quality_weights.split(',')}

    assert os.path.exists(args.fastq), f"ERROR: {args.fastq} file not found"

    start_time = time.time()
    run_alignment_from_fastq(args.fastq, args.n, args.k, quality_weights, args.output)
    elapsed = round(time.time() - start_time, 2)
    print(f"Program completed. Elapsed: {elapsed} sec")

##python lq_swm_fuzzy.py --fastq SRR26136035_2.fastq --k 5 --n 23 --output result_run_1.txt