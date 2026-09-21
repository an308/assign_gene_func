from Bio.Align import substitution_matrices
import numpy as np

GAP_PENALTY = -11 #check if this is correct

def global_alignment(seq1, seq2, scoring_function):
    """Global sequence alignment using the Needleman–Wunsch algorithm.

    Indels should be denoted with the "-" character.

    Parameters
    ----------
    seq1: str
        First sequence to be aligned.
    seq2: str
        Second sequence to be aligned.
    scoring_function: Callable

    Returns
    -------
    str
        First aligned sequence.
    str
        Second aligned sequence.
    float
        Final score of the alignment.

    Examples
    --------
    >>> global_alignment("abracadabra", "dabarakadara", lambda x, y: [-1, 1][x == y])
    ('-ab-racadabra', 'dabarakada-ra', 5.0)

    Other alignments are not possible.

    """
    num_rows = len(seq2) + 1
    num_cols = len(seq1) + 1
    nw_matrix = np.zeros((num_rows, num_cols))

    nw_matrix[0, :] = np.arange(num_cols) * GAP_PENALTY
    nw_matrix[:, 0] = np.arange(num_rows) * GAP_PENALTY

    traceback = np.zeros((num_rows, num_cols), dtype=object)

    # populate nw_matrix with scores & traceback matrix with predecessor's coordinates

    for i in range(1, num_rows):
        for j in range(1, num_cols):
            match = nw_matrix[i-1, j-1] + scoring_function(seq2[i-1], seq1[j-1])
            insertion = nw_matrix[i-1, j] + GAP_PENALTY
            deletion = nw_matrix[i, j-1] + GAP_PENALTY

            max_score = max(match, insertion, deletion)
            nw_matrix[i, j] = max_score
            
            if deletion == max_score:
                traceback[i, j] = (i, j-1)  
            elif insertion == max_score:
                traceback[i, j] = (i-1, j)
            elif match == max_score:
                traceback[i, j] = (i-1, j-1)   

    i = num_rows - 1
    j = num_cols - 1
    aligned_seq1 = []
    aligned_seq2 = []

    # traceback to get alignemnts

    while (i > 0 or j > 0):

        if i == 0: # reached top row, need to move left (deletion) to get to 0,0
            aligned_seq1.append(seq1[j-1])
            aligned_seq2.append("-")
            j -= 1
        elif j == 0: # reached leftmost (1st) col, need to move up (insertion) to get to 0,0
            aligned_seq1.append("-")
            aligned_seq2.append(seq2[i-1])
            i -= 1

        else:
            ip, jp = traceback[i, j]

            if i == ip: # deletion --> move left (i didn't change)
                aligned_seq1.append(seq1[j-1])
                aligned_seq2.append("-")
                j = jp
            elif j == jp: # insertion --> move up (j didn't change)
                aligned_seq1.append("-")
                aligned_seq2.append(seq2[i-1])
                i = ip
            else: # match/mismatch --> move diagonally (both i & j changed)
                aligned_seq1.append(seq1[j-1])
                aligned_seq2.append(seq2[i-1])
                i = ip
                j = jp

    aligned_seq1.reverse()
    aligned_seq2.reverse()
    
    aligned_seq1 = "".join(aligned_seq1)
    aligned_seq2 = "".join(aligned_seq2)

    return (aligned_seq1, aligned_seq2, nw_matrix[-1, -1])




def local_alignment(seq1, seq2, scoring_function):
    """Local sequence alignment using the Smith-Waterman algorithm.

    Indels should be denoted with the "-" character.

    Parameters
    ----------
    seq1: str
        First sequence to be aligned.
    seq2: str
        Second sequence to be aligned.
    scoring_function: Callable

    Returns
    -------
    str
        First aligned sequence.
    str
        Second aligned sequence.
    float
        Final score of the alignment.

    Examples
    --------
    >>> local_alignment("pending itch", "unending glitch", lambda x, y: [-1, 1][x == y])
    ('ending --itch', 'ending glitch', 9.0)

    Other alignments are not possible.

    """
    raise NotImplementedError()



def scoring_function(aa_i,aa_j):
    sub_matrix = substitution_matrices.load("BLOSUM62")
    return(sub_matrix[aa_i.upper(), aa_j.upper()])


# testing (GAP_PENALTY = -1 and simple_scoring_function)

# def main():
#     x, y, score = global_alignment("abracadabra", "dabarakadara", scoring_function)
#     print(x, y, score)

# if __name__ == "__main__":
#     main()
