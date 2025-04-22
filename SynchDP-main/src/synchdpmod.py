import pandas as pd
import numpy as np
import copy
import math
from enum import IntEnum
from scipy import stats

class Trace(IntEnum):
    STOP = 0
    LEFT = 1
    UP = 2
    DIAGONAL = 3

class synchDP:
    def __init__(self, data1: pd.DataFrame, data2: pd.DataFrame, penalty: float, pcut: float, window: int, ecut: float = None, candidate: int = 10):
        """
        Initialize SynchDP with sequence data and parameters.
        
        Args:
            data1: Sequence data for alignment (DataFrame format).
            data2: Sequence data for alignment (DataFrame format).
            penalty: Penalty for time gap scoring.
            pcut: Pearson-correlation cut-off for valid alignments.
            window: Window size for local alignment.
            ecut: Euclidean distance cut-off (optional).
        """
        # Parameters
        self.penalty = penalty
        self.candidate = candidate
        self.pcut = pcut
        self.window = window
        self.ecut = ecut
        self.data1, self.data2 = data1, data2
        
        # Sequence Initialization
        self.seq1, self.seq1_tgs = self._calculate_tgs(self.data1)
        self.seq2, self.seq2_tgs = self._calculate_tgs(self.data2)
        self.seq1, self.seq1_tgs, self.seq2, self.seq2_tgs = self._order_sequences(self.seq1, self.seq1_tgs, self.seq2, self.seq2_tgs)

        # Matrices
        self.row, self.col = len(self.seq1) + 1, len(self.seq2) + 1
        self.cumulative_matrix, self.backtrace_matrix, self.length_matrix = self._compute_alignment_matrix()

        # Average Matrix
        self.average_matrix = np.divide(
            self.cumulative_matrix, self.length_matrix, 
            out=np.zeros_like(self.cumulative_matrix), 
            where=self.length_matrix != 0 
            )
        self.average_matrix[np.isnan(self.average_matrix)] = 0
        
        # Results
        self.results = self._perform_backtrace()
        
    def _calculate_tgs(self, data: pd.DataFrame) -> tuple:
        """
        Calculate Time Gap Scores (TGS) for a given sequence.
        """
        tgs = []
        current_time = int(data.columns[0])
        for time in data.columns[1:]:
            tgs.append((int(time) - current_time) * self.penalty)
            current_time = int(time)
        sequence = data.iloc[0].to_list()
        return sequence, tgs
    
    def _order_sequences(self, seq1, seq1_tgs, seq2, seq2_tgs) -> tuple:
        """
        Ensure that seq1 is shorter or equal in length to seq2.
        """
        if len(seq1) <= len(seq2):
            return seq1, seq1_tgs, seq2, seq2_tgs
        return seq2, seq2_tgs, seq1, seq1_tgs
    
    def _compute_alignment_matrix(self):
        """
        Compute alignment matrix, backtrace matrix, and length matrix.
        """
        cumulative_matrix = np.zeros((self.row, self.col), dtype=float)
        backtrace_matrix = np.zeros((self.row, self.col), dtype=int)
        length_matrix = np.zeros((self.row, self.col), dtype=float)

        for i in range(self.window, self.row):
            for j in range(self.window, self.col):
                x_window, y_window = self.seq1[i - self.window:i], self.seq2[j - self.window:j]

                # Pearson Correlation
                if np.std(x_window) == 0 and np.std(y_window) == 0:
                    pscore = 1
                elif np.std(x_window) == 0 or np.std(y_window) == 0:
                    pscore = 0
                else:
                    pscore = stats.pearsonr(x_window, y_window)[0]

                # Matching Score Calculation
                ip = abs(round(sum(self.seq2_tgs[j - self.window:j - 1]) - sum(self.seq1_tgs[i - self.window:i - 1]), 3))
                match_score = pscore - ip

                # Scores for alignment
                diagonal_score = cumulative_matrix[i - 1, j - 1] + match_score
                vertical_score = cumulative_matrix[i - 1, j] - self.seq1_tgs[i - 2]
                horizontal_score = cumulative_matrix[i, j - 1] - self.seq2_tgs[j - 2]

                # Update matrices
                cumulative_matrix[i, j] = max(0, diagonal_score, vertical_score, horizontal_score)
                if cumulative_matrix[i, j] == diagonal_score:
                    backtrace_matrix[i, j] = Trace.DIAGONAL
                    length_matrix[i, j] = length_matrix[i - 1, j - 1] + 1
                elif cumulative_matrix[i, j] == vertical_score:
                    backtrace_matrix[i, j] = Trace.UP
                    length_matrix[i, j] = length_matrix[i - 1, j]
                elif cumulative_matrix[i, j] == horizontal_score:
                    backtrace_matrix[i, j] = Trace.LEFT
                    length_matrix[i, j] = length_matrix[i, j - 1]
                else:
                    backtrace_matrix[i, j] = Trace.STOP
                    length_matrix[i, j] = 0
        return cumulative_matrix, backtrace_matrix, length_matrix

    def _perform_backtrace(self):
        """
        Perform backtrace to extract alignment results, considering `ecut` if provided.
        """
        results = []
        temp_average_matrix = copy.deepcopy(self.average_matrix)

        while np.any(temp_average_matrix > 0):
            # Find the maximum average score in the matrix
            max_index = np.argmax(temp_average_matrix)
            i, j = divmod(max_index, temp_average_matrix.shape[1])
            avg_score = temp_average_matrix[i, j]

            path = []
            while self.backtrace_matrix[i, j] != Trace.STOP:
                temp_average_matrix[i, j] = 0  # Mark the cell as visited
                path.append((i, j))

                if self.cumulative_matrix[i, j] < (self.pcut * self.length_matrix[i, j]):
                    path.pop()  # Discard the current path if it doesn't meet the cut-off
                    break

                # Move based on the backtrace direction
                if self.backtrace_matrix[i, j] == Trace.DIAGONAL:
                    i, j = i - 1, j - 1
                elif self.backtrace_matrix[i, j] == Trace.LEFT:
                    j -= 1
                elif self.backtrace_matrix[i, j] == Trace.UP:
                    i -= 1
            path = path[::-1]
            
            if path:
                # Add initial window to the path
                x, y = path[0]
                for _ in range(self.window - 1):  # Add previous path by window size
                    x -= 1
                    y -= 1
                    wind_path = (x, y)
                    path = [wind_path] + path

                align_score = avg_score * math.log2(len(path))
                path_length = len(path)
                e_dist_sum = 0

                # Calculate Euclidean distance (if needed)
                for k, (x_idx, y_idx) in enumerate(path):
                    x_val = self.seq1[x_idx - 1] if x_idx > 0 else 0
                    y_val = self.seq2[y_idx - 1] if y_idx > 0 else 0
                    e_dist_sum += x_val - y_val

                self.e_score = abs(e_dist_sum) / path_length

                if len(self.data1.columns) > len(self.data2.columns):
                    path = [(y, x) for x, y in path]
                    
                # Evaluate ecut if provided
                if self.ecut is None or self.e_score < self.ecut:
                    results.append([align_score, avg_score, path])
                
        # results.sort(key=lambda x: -x[0])
        # results = results[:self.candidate]
        results.sort(key= lambda x:-x[1])
        return results
