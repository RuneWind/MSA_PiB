from MSA import MSAligner
from MSA import pairwise_aligment
from MSA import get_col
from MSA import fasta_to_lists
from Bio import SeqIO
from operator import itemgetter
from copy import deepcopy
import os

#################################### FUNCTION DEFINITIONS ####################################

### Function for extracting gap regions
def get_gap_regions(MSA):
    gap_idx = []
    
    i = 0
    while i < len(MSA[0]):
        if MSA[0][i] != "-":
            i += 1
        else:
            start_idx = i
            i += 1
            while i < len(MSA[0]) and MSA[0][i] == "-":
                i += 1
            end_idx = i
            if end_idx - start_idx > 1:
                gap_idx.append([start_idx, end_idx])
    return gap_idx


### Function for checking if sequences in a gap region are all gaps and return their index if they are
def get_seq_idx(gap_region):
    seq_idx = []
    for i in range(len(gap_region)):
        if set(gap_region[i]) != {"-"}:
            seq_idx.append(i)
    if len(seq_idx) >= 2:
        return seq_idx
    else:
        return False


### Function for recursively aligning gap regions
def align_rec(msa, gap_region, sm, gc, start, end):
    # Get index of sequences that contain sequence characters
    seq_gap_idx = get_seq_idx(gap_region)
    if seq_gap_idx:
        gap_region = [gap_region, seq_gap_idx]

        # Extract the sequences from the gap region that contain sequence characters
        align_list = list(itemgetter(*gap_region[1])(gap_region[0]))
        # Make index map between alignment list and gap_regions
        align_map = {}
        for i in range(len(align_list)):
            for j in range(len(gap_region[0])):
                if align_list[i] == gap_region[0][j]:
                    align_map[i] = j
                    break
        # Remove gaps from sequences to be aligned
        align_list[:] = [[x for x in l if x != "-"] for l in align_list]

        ### If there are only two sequences do pairwise alignment
        if len(align_list) == 2:
            alignment = pairwise_aligment(align_list[0], align_list[1], sm, gc)[0]
            # Reinsert gap sequences into alignment
            j = 0
            result_alignment = []
            for i in range(len(gap_region[0])):
                if i in gap_region[1]:
                    result_alignment.append(alignment[j])
                    j += 1
                else:
                    result_alignment.append(["-" for a in range(len(alignment[0]))])
            # Return msa with result alignment reinserted
            return reinsert_alignment(msa, result_alignment, start, end)
        
        ### If there are more than two sequences do multiple alignment
        elif len(align_list) > 2:
            alignment, cs = MSAligner(align_list, sm, gc)
            gap_idx = get_gap_regions(alignment)
            if len(gap_idx) > 0:
                for i in range(len(gap_idx)-1, -1, -1):
                    gap = gap_idx[i]
                    sub_gap_region =  get_col(alignment[1:], gap[0], gap[1])                 
                    alignment = align_rec(alignment, sub_gap_region, sm, gc, gap[0], gap[1])
                    # Reinsert gap sequences into alignment
                    result_alignment = []
                    j = 1
                    for i in range(len(gap_region[0])):
                        if i == align_map[cs]:
                            result_alignment.append(alignment[0])
                        elif i in gap_region[1]:
                            result_alignment.append(alignment[j])
                            j += 1
                        else:
                            result_alignment.append(["-" for a in range(len(alignment[0]))])
                    # Return msa with result alignment reinserted
                    reinsert_alignment(msa, result_alignment, start, end)
            return msa
    # If there are too few sequences containing sequence characters just return the MSA
    else:
        return msa


### Function for reinserting sub alignment into larger alignment (full_align)
def reinsert_alignment(full_align, sub_align, start, end):
    # Insert extra gaps if the new sub-alignment is larger than the old gap
    for i in range(len(sub_align[0])-(end-start)):
        for l in full_align:
            l = l.insert(end, "-")

    for i in range(1, len(full_align)):
        for j in range(start, start+len(sub_align[0])):
            full_align[i][j] = sub_align[i-1][j-start]
    return full_align


### Function for gathering a MSA using the implemented functions of the extension algorithm
def MSA_ext(S, sm, gc):
    # Run first iteration on sequences without gaps
    msa = MSAligner(S, sm, gc)[0]
    
    gap_idx = get_gap_regions(msa)
    if len(gap_idx) > 0:
        for i in range(len(gap_idx)-1, -1, -1):
            gap = gap_idx[i]
            gap_region =  get_col(msa[1:], gap[0], gap[1])
            align_rec(msa, gap_region, sm, gc, gap[0], gap[1])
    return msa