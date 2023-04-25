"""Needleman-Wunsch Sequence Alignment Algorithm."""
from typing import Tuple, Union
import numpy as np


def needleman_wunsch(
    seq1: str,
    seq2: str,
    match_score: int,
    mismatch_score: int,
    gap_penalty: int,
    verbose: bool = False,
) -> Union[str, Tuple[int, str, str]]:
    """Input 2 sequences and get the needleman_wunsch match.

    N = length sequence 1
    M = length sequence 2
    The overall complextiy is O(N*M)

    The two parts with non O(1) complexity are the nested for loop and the
    while loop.

    The nested for loop is O(N*M) and the while loop is O(N+M) as explained
    below.
    """
    # get the lengths
    seq1_length = len(seq1)
    seq2_length = len(seq2)

    # create a matrix with all 0s. Needs to be the seq1 + 1 by seq2 + 1
    score_matrix = np.zeros((seq1_length + 1, seq2_length + 1), dtype=int)

    # fill in the first column and row with the appropriate gap penalties
    score_matrix[0, :] = np.arange(seq2_length + 1) * gap_penalty
    score_matrix[:, 0] = np.arange(seq1_length + 1) * gap_penalty

    # go through each row and column
    # this is a nested for loop so time complexity is O(N*M)
    # it will always be O(N*M) because we are filling in the values of a
    # NxM matrix.
    for i in range(1, seq1_length + 1):
        for j in range(1, seq2_length + 1):
            # check for match
            if seq1[i-1] == seq2[j-1]:
                # diagonal add match score
                diagonal = score_matrix[i-1][j-1] + match_score
            else:
                # diagonal add penalty
                diagonal = score_matrix[i-1][j-1] + mismatch_score
            # check insert and deletion add penalty
            delete = score_matrix[i-1][j] + gap_penalty
            insert = score_matrix[i][j-1] + gap_penalty

            # get the max value
            score_matrix[i][j] = max(diagonal, insert, delete)

    # store the alginment score from bottom right of matrix
    alignment_score: int = score_matrix[seq1_length][seq2_length]
    # placeholder strings for our alginment
    aligned_seq1 = ""
    aligned_seq2 = ""
    match_string = ""

    # we need traceback values that start at the length
    trc1 = seq1_length
    trc2 = seq2_length

    # go through the matrix tracing back until we get to 0 in
    # one of these
    #
    # the time complexity here is the length of the output sequence
    # The worst case scenario could be if we somehow we had an output with no
    # overlap leading to an output that is the size of both seqs
    # this could happen especially if the gap penalty was lower than mismatch
    # like
    # AAAAAA-----
    # ------CCCCC
    # so the time complextiy is O(N+M)
    while trc1 > 0 or trc2 > 0:
        # check if it was a diagonal move
        if (trc1 > 0 and trc2 > 0 and score_matrix[trc1][trc2] ==
            score_matrix[trc1-1][trc2-1] + (match_score if seq1[trc1-1] ==
                                            seq2[trc2-1] else mismatch_score)):
            # add values to both strings
            aligned_seq1 = seq1[trc1-1] + aligned_seq1
            aligned_seq2 = seq2[trc2-1] + aligned_seq2
            # if they are the same
            if seq1[trc1-1] == seq2[trc2-1]:
                # indicate match in match string
                match_string = '|' + match_string
            else:
                # indicate a msimatch in the match string
                match_string = '.' + match_string
            trc1 -= 1
            trc2 -= 1

        # check if it was a vertical move
        elif (trc1 > 0 and score_matrix[trc1][trc2] ==
                score_matrix[trc1-1][trc2] + gap_penalty):
            # add value to seq 1
            aligned_seq1 = seq1[trc1-1] + aligned_seq1
            # add - to seq 2
            aligned_seq2 = "-" + aligned_seq2
            # match string is open
            match_string = ' ' + match_string
            trc1 -= 1

        # else it was a horizontal move
        else:
            # add - to seq 1
            aligned_seq1 = "-" + aligned_seq1
            # add value to seq 2
            aligned_seq2 = seq2[trc2-1] + aligned_seq2
            # match string stays open
            match_string = ' ' + match_string
            trc2 -= 1

    # by default this is False
    # if verbose output the value in a nice string format
    # else output the values as a tuple
    if verbose:
        string_out = ("Alignment Score: " + str(alignment_score) + "\n" +
                      str(aligned_seq1) + "\n" + str(match_string) + "\n" +
                      str(aligned_seq2))
        return string_out
    else:
        return alignment_score, aligned_seq1, aligned_seq2
