"""Multiple Sequence Alignment Algorithm."""
from typing import Tuple, Union, List, Optional
import numpy as np
from Needleman_wunsch import needleman_wunsch


def _alignments_to_profile(sequences: list[str]) -> dict[str, list[int]]:
    """Take a set of alignments and turn it into a profile."""
    # get length
    length = len(sequences[0])

    # initialize count matrix with zeros
    counts = {
        "A": [0] * length,
        "C": [0] * length,
        "G": [0] * length,
        "T": [0] * length,
        "-": [0] * length,
    }

    # count the frequency of each nucleotide or gap
    for i in range(length):
        for seq in sequences:
            counts[seq[i]][i] += 1
    return counts


def _profile_needleman_wunsch(
    profile1: dict[str, list[int]],
    profile2: dict[str, list[int]],
    gap_penalty: int = -10,
) -> Tuple[List[List[int]], List[List[Union[int, str]]]]:
    """Modification of the needleman-wunsch algorithm for profiles."""
    # Create a matrix filled with zeros
    rowlengths: list[int] = profile1["A"]
    collengths: list[int] = profile2["A"]
    rows = len(rowlengths)
    cols = len(collengths)

    # get the number of sequences in each profile
    num_seqs_1 = 0
    for list in profile1.items():
        num_seqs_1 += list[1][0]
    num_seqs_2 = 0
    for list in profile2.items():
        num_seqs_2 += list[1][0]
    # Create the score_matrix
    score_matrix = [[0 for j in range(cols + 1)] for i in range(rows + 1)]

    # create the direction matrix
    # 1 = diagonal
    # 2 = horizontal
    # 3 = vertical
    direction_matrix: List[List[Union[int, str]]] = [
        [0 for j in range(cols + 1)] for i in range(rows + 1)
    ]

    # Fill the first row and column with gap penalties
    for i in range(1, rows + 1):
        score_matrix[i][0] = gap_penalty * i
    for j in range(1, cols + 1):
        score_matrix[0][j] = gap_penalty * j
    # fill in matrix similar to how we do it in the needleman-wunsch
    # we will use the profiles to make it so that the match score is greater
    # if there is a greater level of agreement between the sequences
    for i in range(1, rows + 1):
        for j in range(1, cols + 1):
            # start with a match score of 0
            match_score = 0
            # go through each nucleotide position and check each key
            for nucleotide in profile1.keys():
                # for each key check for match mismatch
                # current cell values (num seqs with this nucleotide)
                profile1_score = profile1[nucleotide][i - 1]
                profile2_score = profile2[nucleotide][j - 1]

                # num seqs with different nucleotide
                profile1_mismatch = num_seqs_1 - profile1_score
                profile2_mismatch = num_seqs_2 - profile2_score

                # add score for all the matches
                # remove score for all the mismatches
                match_score += profile1_score * profile2_score
                match_score += profile1_score * profile2_mismatch * -1
                match_score += profile2_score * profile1_mismatch * -1
            m_score = score_matrix[i - 1][j - 1] + match_score
            v_score = score_matrix[i - 1][j] + gap_penalty
            h_score = score_matrix[i][j - 1] + gap_penalty

            score = max(m_score, h_score, v_score)
            score_matrix[i][j] = score

            if score == m_score:
                direction_matrix[i][j] = "d"
            elif score == h_score:
                direction_matrix[i][j] = "h"
            else:
                direction_matrix[i][j] = "v"
    return score_matrix, direction_matrix


def _add_gaps(
    direction_matrix: list[List[int | str]],
    prof1_seqs: list[str],
    prof2_seqs: list[str],
) -> tuple[list[str], list[str]]:
    """Change profile sequences after needleman-wunsch for profiles."""
    # create empty strings to fill with gapped seqs
    prof1_gapped_seqs = [""] * len(prof1_seqs)
    prof2_gapped_seqs = [""] * len(prof2_seqs)

    # start in bottom right corner
    i = len(direction_matrix) - 1
    j = len(direction_matrix[0]) - 1

    # while we are not at top or left
    while direction_matrix[i][j] != 0:
        if direction_matrix[i][j] == "d":
            # Go through and add a value in profile 1 seqs
            for k in range(len(prof1_gapped_seqs)):
                prof1_gapped_seqs[k] = (
                    prof1_seqs[k][i - 1] + prof1_gapped_seqs[k]
                )
            # go through and add a value in profile 2 seqs
            for k in range(len(prof2_gapped_seqs)):
                prof2_gapped_seqs[k] = (
                    prof2_seqs[k][j - 1] + prof2_gapped_seqs[k]
                )
            i -= 1
            j -= 1
        elif direction_matrix[i][j] == "h":
            # Go through and add a value in profile 1 seqs
            for k in range(len(prof1_gapped_seqs)):
                prof1_gapped_seqs[k] = "-" + prof1_gapped_seqs[k]
            # go through and add a value in profile 2 seqs
            for k in range(len(prof2_gapped_seqs)):
                prof2_gapped_seqs[k] = (
                    prof2_seqs[k][j - 1] + prof2_gapped_seqs[k]
                )
            j -= 1
        else:  # direction_matrix[i][j] == 'v'
            # Go through and add a value in profile 1 seqs
            for k in range(len(prof1_gapped_seqs)):
                prof1_gapped_seqs[k] = (
                    prof1_seqs[k][i - 1] + prof1_gapped_seqs[k]
                )
            # go through and add a value in profile 2 seqs
            for k in range(len(prof2_gapped_seqs)):
                prof2_gapped_seqs[k] = "-" + prof2_gapped_seqs[k]
            i -= 1
    return prof1_gapped_seqs, prof2_gapped_seqs


def _n_n_alignment(
    seqs1: list[str], seqs2: list[str]
) -> tuple[list[str], list[str]]:
    """Align any n sequences with any m sequences."""
    prof1 = _alignments_to_profile(seqs1)
    prof2 = _alignments_to_profile(seqs2)

    score_matrix, direction_matrix = _profile_needleman_wunsch(prof1, prof2)
    out1, out2 = _add_gaps(direction_matrix, seqs1, seqs2)

    return out1, out2


def multiple_alignment(seqs: list[str]) -> list[str]:
    """Return an alignment of any n sequences."""
    # list to store pairwise alignments
    pairwise_alignments: list[list[Union[int, str]]] = []

    # iterate over all possible pairwise combinations of sequences
    for i in range(len(seqs)):
        for j in range(i + 1, len(seqs)):
            seq1 = seqs[i]
            seq2 = seqs[j]
            # apply Needleman-Wunsch

            # needleman_wunsch can return a string if verbose, but we set
            # verbose to false so it will always return
            # Tuple[int, str, str]
            result: Tuple[int, str, str] = needleman_wunsch(  # type: ignore
                seq1, seq2, 1, -1, -1, verbose=False
            )
            score, alignment_seq1, alignment_seq2 = result
            # add to pairwise alignments
            pairwise_alignments.append(
                [score, i, j, alignment_seq1, alignment_seq2]
            )
    # sort by scores
    sort_pairs = sorted(pairwise_alignments, key=lambda x: x[0], reverse=True)

    sort_pairs_mod = sort_pairs

    # get the base seqs and the indices
    base: list[Union[str, int]] = sort_pairs[0]
    b1 = str(base[3])
    b2 = str(base[4])
    base_seqs: list[str] = [b1, b2]
    bi1 = int(base[1])
    bi2 = int(base[2])
    base_seqs_index: list[int] = [bi1, bi2]

    # remove base profile and seqs from mod list
    del sort_pairs_mod[0]

    while len(sort_pairs_mod) > 0:
        # filter the list to only include adjacent alignments
        filtered_list = [
            sublist
            for sublist in sort_pairs_mod
            if sublist[1] in base_seqs_index or sublist[2] in base_seqs_index
        ]

        # create the profile for joining
        # we want only the seq that is not already in the base profile
        # so check which one is in base seqs index and return the other
        if filtered_list[0][1] in base_seqs_index:
            enter_seq = filtered_list[0][4]
            enter_index = filtered_list[0][2]
        else:
            enter_seq = filtered_list[0][3]
            enter_index = filtered_list[0][1]
        # do a progressive alignment with new seq
        if isinstance(enter_seq, str):
            base_aligned, enter_aligned = _n_n_alignment(
                base_seqs, [enter_seq]
            )
        # add aligned enter to base
        base_seqs = base_aligned + enter_aligned

        if isinstance(enter_index, int):
            base_seqs_index = base_seqs_index + [enter_index]
        # eliminate any alignments that would be redundant
        # we only want the alignments that contain a sequence that has not
        # already entered the base
        sort_pairs_mod = [
            sublist
            for sublist in sort_pairs_mod
            if sublist[1] not in base_seqs_index
            or sublist[2] not in base_seqs_index
        ]
    return base_seqs


def main() -> None:
    """Test the functions."""
    print(multiple_alignment(["ACTGTCA", "ACTTCA", "ACTGTA"]))


if __name__ == "__main__":
    main()
