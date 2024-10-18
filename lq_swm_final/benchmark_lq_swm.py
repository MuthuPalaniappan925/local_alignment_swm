import numpy as np
import random
import time
import argparse
from typing import Tuple, Dict, Any

class ScoringSystem:
    def __init__(self, match: int = 2, mismatch: int = -1, gap: int = -2, quality_weights: Dict[str, int] = None) -> None:
        self.match = match
        self.mismatch = mismatch
        self.gap = gap
        self.quality_weights = quality_weights or {
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

    def _fuzzy_similarity(self, a: str, b: str, quality_a: str) -> float:
        weight_a = self.quality_weights.get(quality_a, 1)
        return weight_a

    def _default_scoring(self, a: str, b: str, quality_a: str) -> float:
        if a == '-' or b == '-':
            return self.gap
        fuzzy_score = self._fuzzy_similarity(a, b, quality_a)
        return self.match * fuzzy_score if a == b else self.mismatch * fuzzy_score

    def score(self, a: str, b: str, quality_a: str) -> float:
        return self._default_scoring(a, b, quality_a)

class SequencesAnalyzer:
    traceback_symbols = {0: '↖', 1: '↑', 2: '←'}

    def __init__(self, seq_a: str, seq_b: str, quality_a: str, match: int = 2, mismatch: int = -1, gap: int = -2, quality_weights: Dict[str, int] = None) -> None:
        self.seq_a = seq_a
        self.seq_b = seq_b
        self.quality_a = quality_a
        self.scoring_sys = ScoringSystem(match, mismatch, gap, quality_weights)

    def smith_waterman_algorithm(self) -> Dict[str, Any]:
        rows, cols = len(self.seq_a) + 1, len(self.seq_b) + 1
        H = np.zeros((rows, cols), dtype=float)
        traceback = np.zeros((rows, cols), dtype=np.dtype('U5'))
        max_score, max_pos = 0, (0, 0)

        for row in range(1, rows):
            for col in range(1, cols):
                a, b = self.seq_a[row - 1], self.seq_b[col - 1]
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

        return {'result_matrix': H, 'traceback_matrix': traceback, 'score': max_score, 'score_pos': max_pos}

    def _traceback(self, result_matrix, traceback_matrix, start_pos: Tuple[int, int]) -> Tuple[str, str, str]:
        seq_a_aligned, seq_b_aligned, matches = '', '', ''
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

def generate_random_sequence(length: int) -> str:
    return ''.join(random.choice("ACGT") for _ in range(length))

def generate_random_quality(length: int) -> str:
    return ''.join(random.choice("!\"#$%&'()*+,-./0123456789:;<=>?@ABCDEFGHI") for _ in range(length))

def run_alignment_random_sequences(output_filename: str, iterations: int = 70, num_references: int = 5) -> None:
    #lengths = [32, 64, 128, 256]
    lengths = [32, 64]
    
    with open(output_filename, 'w', encoding='utf-8') as f:
        for length in lengths:
            print(f"[INFO] Starting Alignments for Length {length} with {iterations} Iterations and {num_references} Reference Sequences")

            total_alignment_time = 0
            total_traceback_time = 0
            reference_sequences = [generate_random_sequence(length) for _ in range(num_references)]

            for i in range(iterations):
                seq_a = generate_random_sequence(length)
                quality_a = generate_random_quality(length)

                for ref_seq in reference_sequences:
                    start_time_alignment = time.time()
                    analyzer = SequencesAnalyzer(seq_a, ref_seq, quality_a)
                    result = analyzer.smith_waterman_algorithm()
                    alignment_time = round(time.time() - start_time_alignment, 4)
                    start_time_traceback = time.time()
                    alignment_a, alignment_b, matches = analyzer._traceback(
                        result['result_matrix'], result['traceback_matrix'], result['score_pos']
                    )
                    traceback_time = round(time.time() - start_time_traceback, 4)
                    total_alignment_time += alignment_time
                    total_traceback_time += traceback_time

                    ##avoid excessive I/O
                    if (i + 1) % 10 == 0:
                        f.write(f"Iteration {i + 1}/{iterations} - Length {length} with Reference:\n")
                        f.write(f"Input Sequence A: {seq_a}\n")
                        f.write(f"Generated Reference Sequence: {ref_seq}\n")
                        f.write(f"Alignment Score: {result['score']}\n")
                        f.write(f"Avg Alignment Time: {total_alignment_time / (num_references * (i + 1))} sec\n")
                        f.write(f"Avg Traceback Time: {total_traceback_time / (num_references * (i + 1))} sec\n\n")
                        f.write("Alignment:\n")
                        f.write(f"{alignment_a}\n")
                        f.write(f"{matches}\n")
                        f.write(f"{alignment_b}\n\n")
                        f.write("-" * 80 + "\n")
                        print(f"[INFO] Completed {i + 1} Alignments for Length {length} with {num_references} References")

            ##final statistics for this length
            f.write(f"[STAT] For length {length}, Total Alignment Time: {total_alignment_time}, Avg Alignment Time: {total_alignment_time / (iterations * num_references)}\n")
            f.write(f"[STAT] For length {length}, Total Traceback Time: {total_traceback_time}, Avg Traceback Time: {total_traceback_time / (iterations * num_references)}\n")

def main():
    parser = argparse.ArgumentParser(description="Sequence Alignment Experiment Script")
    parser.add_argument("--output_filename", type=str, default="alignment_results.txt", help="Name of the output file")
    parser.add_argument("--iterations", type=int, default=70, help="Number of iterations")
    parser.add_argument("--num_references", type=int, default=5, help="Number of reference sequences")
    args = parser.parse_args()
    run_alignment_random_sequences(args.output_filename, args.iterations, args.num_references)

if __name__ == "__main__":
    main()
