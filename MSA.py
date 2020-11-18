# 2SP approximation for Multiple Sequence Alignment
# Rune Wind - Fall 2020
from Bio import SeqIO


#################################### FUNCTION DEFINITIONS ####################################

#########
#### Pairwise alignment
#########

### Function for reading fasta files
def read_fasta(path):
    fasta_to_lists("/Users/runewind/Downloads/brca1-testseqs.fasta")

### Function for calculating the optimal score for a given cell
def calc_score(seq1, seq2, matrix, sm, gc, i, j):
    diag = matrix[i-1][j-1] + sm[seq1[i-1]][seq2[j-1]]
    up = matrix[i-1][j] + gc
    left = matrix[i][j-1] + gc

    return min(diag, up, left)


### Function for computing the alignment matrix
def create_align_matrix(seq1, seq2, gc, sm):
    # Prime a matrix consisting of lists of lists and with gapcost in first row and column and None in all other columns
    matrix = [[i * gc if j == 0 else j * gc if i == 0 else None for j in range(len(seq2)+1)] for i in range(len(seq1)+1)]

    # Fill out the rest of the matrix with scores
    for i in range(1, len(seq1)+1):
        for j in range(1, len(seq2)+1):
            matrix[i][j] = calc_score(seq1, seq2, matrix, sm, gc, i, j)

    return matrix


### Function for finding an optimal pairwise alignment by backtracking
def back_tracking(seq1, seq2, matrix, sm, gc):
    seq1_align = []
    seq2_align = []
    
    i = len(seq1)
    j = len(seq2)

    while i > 0 or j > 0:
        # backtrack in diagonal direction
        if i > 0 and j > 0 and (matrix[i][j] == matrix[i-1][j-1] + sm[seq1[i-1]][seq2[j-1]]):
            seq1_align.append(seq1[i-1])
            seq2_align.append(seq2[j-1])
            i -= 1
            j -= 1
        # Backtrack upwards
        elif i > 0 and j >= 0 and (matrix[i][j] == matrix[i-1][j] + gc):
            seq1_align.append(seq1[i-1])
            seq2_align.append("-")
            i -= 1
        # Backtrack left
        elif i >= 0 and j > 0 and (matrix[i][j] == matrix[i][j-1] + gc):
            seq1_align.append("-")
            seq2_align.append(seq2[j-1])
            j -= 1

    return [seq1_align[::-1], seq2_align[::-1]]    


### Function calling the other functions to make alignment in one function call
def pairwise_aligment(seq1, seq2, sm, gc):
    matrix = create_align_matrix(seq1, seq2, gc, sm)
    alignments = back_tracking(seq1, seq2, matrix, sm, gc)

    return [alignments, matrix[-1][-1]]


##########
#### Multiple alignment
#########


### Function for extracting columns from list of alignments
def get_col(lst, i, j):
    return [item[i:j] for item in lst]


### Function for finding the center string
def find_center_string(S, sm , gc):
    pa_scores = [[0 for j in range(len(S))] for i in range(len(S))]

    for i in range(len(pa_scores)):
        for j in range(i):
            align_score = pairwise_aligment(S[i], S[j], sm, gc)[1]
            pa_scores[i][j] = align_score
            pa_scores[j][i] = align_score
    
    sum_scores = [sum(l) for l in pa_scores]
    
    return min(zip(sum_scores, range(len(sum_scores))))[1]
        

### Function for assembling pairwise alignments into MSA
def extend_M_with_A(M, A):
    Mi = []
    i = 0 # Index for A
    j = 0 # Index for M
    
    while i < len(A[0]) or j < len(M[0]):
        
        # In case of gap in Cs in A
        if i < len(A[0]) and A[0][i] == "-":
            if j >= len(M[0]) or M[0][j] != "-": # If there is not already a gap in the top seq, insert one
                M = [m[:j] + ["-"] + m[j:] for m in M]
            Mi.append(A[1][i])
            i += 1
            j += 1
        
        # In case of gap in M
        elif j < len(M[0]) and M[0][j] == "-":
            gap_count = 0
            while j < len(M[0]) and M[0][j] == "-":
                gap_count += 1
                j += 1
            for g in range(gap_count):
                Mi.append("-")
         
        # In case of match between A and M
        elif A[0][i] == M[0][j]:
            Mi.append(A[1][i])
            i += 1
            j += 1

    M.append(Mi)
    return M


### Function for gathering the full MSA
def MSAligner(S, sm, gc):
    # print("Finding center string...", end = "", flush=True)
    cs = find_center_string(S, sm, gc)
    # print("\tDone", flush=True)
    M = [S[cs]]

    # print("\nPairwise alignment of sequences...", flush=True)
    for i in range(len(S)):
        if i == cs:
            continue
        elif len(M) == 1:
            # print("\tCenter string and seq", i, flush=True)
            M = pairwise_aligment(S[cs], S[i], sm, gc)[0]
        else:
            # print("\n\tCenter string and seq", i, flush=True)
            A = pairwise_aligment(S[cs], S[i], sm, gc)[0]
            # print("\tExtending MSA", flush=True)
            M = extend_M_with_A(M, A)
            # print("\tDone", flush=True)
    return [M, cs]


### Function for reading fasta and outputting the sequences as lists of lists
def fasta_to_lists(file_name):
    fa = SeqIO.parse(file_name, "fasta")

    seq_list = []
    nucleic_list = ["U", "W", "S", "M", "K", "R", "Y", "B", "D", "H", "V", "N", "Z"]
    
    for record in fa:
        corrected_seq = str(record.seq)
        for symbol in nucleic_list:
            corrected_seq = corrected_seq.replace(symbol, "A")
        seq_list.append(corrected_seq)

    return [[a for a in seq]for seq in seq_list]