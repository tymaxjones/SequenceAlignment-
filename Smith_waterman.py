"""
This is a python implementation of the Smith-Waterman algorithm.

The Smith-Waterman algorithm is a dynamic programming algorithm for
finding regions of similarity between two strings. It is used in
bioinformatics to find regions of local similarity between sequences.
"""

import numpy as np
from typing import Union, Tuple


def smith_waterman(
    seq1: str,
    seq2: str,
    match_score: int,
    mismatch_score: int,
    gap_penalty: int,
) -> Tuple[int, int, str, str, str]:
    """Smith Waterman algorithm to find the local alignment of two sequences.

    Smith-waterman algorithm is a dynamic programming algorithm to find the \
        best local alignment of two sequences.
    The length of sequence a is m.
    The length of sequence b is n.
    It has a time complexity of O(m*n) for Traceback using a nexted 'for' loop.
    and space complexity of O(m*n) for construct a m*n scoring matrix and\
          a traceback matrix.
    So the overall time complexity is O(N^2).
    """
    length_1 = len(seq1)  # O(1)
    length_2 = len(seq2)  # O(1)
    max_score = 0  # O(1)

    # Inside function to calculate score
    def cal_score(x: Union[int, str], y: Union[int, str]) -> int:
        if x == y:  # O(1)
            return match_score  # O(1)
        elif x == "-" or y == "-":  # O(1)
            return gap_penalty  # O(1)
        else:  # O(1)
            return mismatch_score  # O(1)

    # Initialize the scoring matrix
    scoring_matrix = np.zeros((length_1 + 1, length_2 + 1))  # O(1)

    # Initialize the traceback matrix
    Track_back_pointer = np.zeros((length_1 + 1, length_2 + 1))  # O(1)

    for i in range(1, length_1 + 1):  # O(n)
        for j in range(1, length_2 + 1):  # O(n)
            score_diagonal = scoring_matrix[i - 1][j - 1] + cal_score(
                seq1[i - 1], seq2[j - 1]
            )  # O(1)
            score_up = scoring_matrix[i][j - 1] + gap_penalty  # O(1)
            score_left = scoring_matrix[i - 1][j] + gap_penalty  # O(1)
            scoring_matrix[i][j] = max(
                0, score_left, score_up, score_diagonal
            )  # O(1)
            if scoring_matrix[i][j] == 0:  # O(1)
                Track_back_pointer[i][j] = 0  # O(1)
            if scoring_matrix[i][j] == score_left:  # O(1)
                Track_back_pointer[i][j] = 1  # O(1)
            if scoring_matrix[i][j] == score_up:  # O(1)
                Track_back_pointer[i][j] = 2  # O(1)
            if scoring_matrix[i][j] == score_diagonal:  # O(1)
                Track_back_pointer[i][j] = 3  # O(1)
            if scoring_matrix[i][j] >= max_score:  # O(1)
                max_i = i  # O(1)
                max_j = j  # O(1)
                max_score = scoring_matrix[i][j]  # O(1)

    # Initialization:
    align1, align2 = "", ""  # O(1)

    i, j = max_i, max_j  # O(1)

    while Track_back_pointer[i][j] != 0:  # O(n)
        if Track_back_pointer[i][j] == 3:  # O(1)
            align1 += seq1[i - 1]  # O(1)
            align2 += seq2[j - 1]  # O(1)
            i -= 1  # O(1)
            j -= 1  # O(1)
        elif Track_back_pointer[i][j] == 2:  # O(1)
            align1 += "-"  # O(1)
            align2 += seq2[j - 1]  # O(1)
            j -= 1  # O(1)
        elif Track_back_pointer[i][j] == 1:  # O(1)
            align1 += seq1[i - 1]  # O(1)
            align2 += "-"  # O(1)
            i -= 1  # O(1)

    # Reverse the sequences
    align1, align2 = align1[::-1], align2[::-1]  # O(n)
    i, j = 0, 0  # O(1)

    symbol = ""  # O(1)
    initial_s = 0  # O(1)
    score = 0  # O(1)
    identity = 0.0  # O(1)
    for i in range(0, len(align1)):  # O(n)
        if align1[i] == align2[i]:  # O(1)
            symbol = symbol + align1[i]  # O(1)
            identity = identity + 1  # O(1)
            score += cal_score(align1[i], align2[i])  # O(1)

        elif (
            align1[i] != align2[i] and align1[i] != "-" and align2[i] != "-"
        ):  # O(1)
            score += cal_score(align1[i], align2[i])  # O(1)
            symbol += " "  # O(1)
            initial_s = 0  # O(1)

        elif align1[i] == "-" or align2[i] == "-":  # O(1)
            symbol += " "  # O(1)
            score += gap_penalty  # O(1)

    identity = float(identity) / len(align1) * 100  # O(1)

    return int(identity), score, align1, symbol, align2  # O(1)
