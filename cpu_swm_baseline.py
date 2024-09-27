## Importing Packages
import numpy as np
from typing import Tuple, Dict, Any
import pandas as pd
import time
import argparse

##Creating a Class for the Scoring System
class ScoringSystem:
    '''This class is responsible for returning the proper edit cost for any AGCT letter combination'''
    def __init__(self, match: int=2, mismatch: int=-1, gap: int=-2) -> None:
        self.match = match
        self.mismatch = mismatch
        self.gap = gap
        self.custom_scoring = None

    def load_csv(self, filename: str) -> None:
        self.custom_scoring = pd.read_csv(filename, header=0, index_col=0, sep=' ')
    
    def _default_scoring(self, a: str, b: str) -> int:
        '''Check if there is a match/mismatch/gap between the letters AGCT in the input files'''
        if a == b:
            return self.match
        elif a == "-" or b == "-":
            return self.gap
        return self.mismatch

    def score(self, a: str, b: str) -> int:
        '''This method should be used by algorithms'''
        assert isinstance(a, str) and isinstance(b, str)
        assert len(a) == 1 and len(b) == 1

        if self.custom_scoring is not None:
            try:
                return self.custom_scoring[a][b]
            except KeyError:
                print(f'WARNING: Key not found. Using defaults: {self}')
                return self._default_scoring(a, b)
        else:
            return self._default_scoring(a, b)
        
    def __str__(self):
        return f'Match: {self.match}, Mismatch: {self.mismatch}, Gap: {self.gap}'

##Craeting a class for the sequence analyzer
class SequencesAnalyzer:
    traceback_symbols = {
        0: '↖',
        1: '↑',
        2: '←'
    }

    def __init__(self, seq_a_file: str, seq_b_file: str, match: int = 2, mismatch: int = -1, gap: int = -2, load_csv: bool = False) -> None:
        ##Load sequences from files
        self.seq_a = self._load_sequence_from_file(seq_a_file)
        self.seq_b = self._load_sequence_from_file(seq_b_file)
        self.scoring_sys = ScoringSystem(match, mismatch, gap)

        if load_csv:
            self.scoring_sys.load_csv('scores.csv')
            print('[Scoring system]\n', self.scoring_sys)

    def _load_sequence_from_file(self, filename: str) -> str:
        '''Loads a sequence from a given file'''
        with open(filename, 'r') as file:
            return file.read().strip()

    def local_alignment(self, output_filename: str) -> Tuple[str, str, str, np.ndarray]:
        start_time = time.time()
        result: Dict[str, Any] = self.smith_waterman_algorithm()
        alignment_a, alignment_b, matches = self._traceback(
            result_matrix=result['result_matrix'],
            traceback_matrix=result['traceback_matrix'],
            start_pos=result['score_pos'])

        temp = len(result['traceback_matrix'])
        for i in range(1, temp):
            result['traceback_matrix'][:, 0][i] = self.seq_a[i - 1]
        for i in range(1, temp):
            result['traceback_matrix'][0, :][i] = self.seq_b[i - 1]

        end_time = time.time()
        elapsed_time = end_time - start_time

        with open(output_filename, "w", encoding='utf-8') as f:
            f.write(f"Input Sequence A: {self.seq_a}\n")
            f.write(f"Input Sequence B: {self.seq_b}\n")
            f.write(f"Alignment Score: {result['score']}\n\n")
            f.write("Alignment:\n")
            f.write(f"{alignment_a}\n")
            f.write(f"{matches}\n")
            f.write(f"{alignment_b}\n\n")
            f.write(f"Result Matrix:\n {result['result_matrix']}\n\n")
            f.write(f"Traceback Matrix:\n {result['traceback_matrix']}\n\n")
            f.write(f"Elapsed Time: {elapsed_time:.4f} seconds\n")

        print(
            f"[Local Alignment] Score={result['score']}\n"
            f"Result Matrix:\n {result['result_matrix']}\n"
            f"Traceback Matrix:\n {result['traceback_matrix']}\n"
            f"Alignment:\n {alignment_a}\n {matches}\n {alignment_b}\n"
        )
        return alignment_a, alignment_b, matches, result['traceback_matrix']

    def smith_waterman_algorithm(self) -> Dict[str, Any]:
        rows, cols = len(self.seq_a) + 1, len(self.seq_b) + 1

        ##Initialize matrices
        H = np.zeros(shape=(rows, cols), dtype=int)
        traceback = np.zeros(shape=(rows, cols), dtype=np.dtype('U5'))

        max_score = 0
        max_pos = (0, 0)

        for row in range(1, rows):
            for col in range(1, cols):
                a = self.seq_a[row - 1]
                b = self.seq_b[col - 1]

                score_diag = H[row - 1, col - 1] + self.scoring_sys.score(a, b)
                score_up = H[row - 1, col] + self.scoring_sys.gap
                score_left = H[row, col - 1] + self.scoring_sys.gap

                H[row, col] = max(0, score_diag, score_up, score_left)
                if H[row, col] == score_diag:
                    traceback[row, col] = self.traceback_symbols[0]  # ↖
                elif H[row, col] == score_up:
                    traceback[row, col] = self.traceback_symbols[1]  # ↑
                elif H[row, col] == score_left:
                    traceback[row, col] = self.traceback_symbols[2]  # ←

                ##Update maximum score and position
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
        seq_a_aligned = ''
        seq_b_aligned = ''
        matches = ''  ##store the vertical lines -> '|' for the matches

        row, col = start_pos

        while result_matrix[row, col] > 0:
            symbol = traceback_matrix[row, col]
            if symbol == '↖':
                seq_a_aligned += self.seq_a[row - 1]
                seq_b_aligned += self.seq_b[col - 1]
                if self.seq_a[row - 1] == self.seq_b[col - 1]:
                    matches += '|'
                else:
                    matches += ' '
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

        return seq_a_aligned[::-1], seq_b_aligned[::-1], matches[::-1]

if __name__ == "__main__":
    ##Setup argument parser
    parser = argparse.ArgumentParser(description="Perform local sequence alignment using Smith-Waterman algorithm.")
    parser.add_argument('seq_a_file', type=str, help="Input file containing sequence A.")
    parser.add_argument('seq_b_file', type=str, help="Input file containing sequence B.")
    parser.add_argument('output_file', type=str, help="Output file to save the alignment result.")
    parser.add_argument('--match', type=int, default=2, help="Score for a match (default: 2).")
    parser.add_argument('--mismatch', type=int, default=-1, help="Score for a mismatch (default: -1).")
    parser.add_argument('--gap', type=int, default=-2, help="Penalty for a gap (default: -2).")
    parser.add_argument('--load_csv', action='store_true', help="Load custom scoring from a CSV file.")

    args = parser.parse_args()

    ##create an instance of SequencesAnalyzer with the provided arguments
    analyzer = SequencesAnalyzer(args.seq_a_file, args.seq_b_file, match=args.match, mismatch=args.mismatch, gap=args.gap, load_csv=args.load_csv)
    
    ##perform local alignment and output results to file
    alignment_a, alignment_b, matches, trac_mat = analyzer.local_alignment(output_filename=args.output_file)
