import os
from MSA import MSAligner
from MSA import fasta_to_lists
from MSA_ext import MSA_ext
from msa_sp_score_3k import compute_sp_score
import pandas as pd
import time

class col:
   GREEN = '\033[92m'
   YELLOW = '\033[93m'
   RED = '\033[91m'
   BOLD = '\033[1m'
   UNDERLINE = '\033[4m'
   END = '\033[0m'

def print_alignment_to_file(seq_list, name):
	x = open(name + "_alignment.fasta", "w")
	for i in range(len(seq_list)):    
		x.write(">seq" + str(i+1) + "\n" + seq_list[i] + "\n")
	x.close()

start_time = time.time()

gc = 3
# Cost matrix for comput_sp_score function
cost = [[0, 5, 2, 5, gc],  # A
		[5, 0, 5, 2, gc],  # C
		[2, 5, 0, 5, gc],  # G
		[5, 2, 5, 0, gc],  # T
		[gc, gc, gc, gc, 0]]  #-'

# gaps cost to use
# gap_costs = [2, 3, 5]
# Diversities
d = [0.7]
# Number of sequences
k = [3, 5, 7, 10]
# Sequence lengths
n = [50, 100, 200, 300, 400, 500, 700, 1000]

# Score matrix
sm = {"A": {"A": 0, "C": 5, "G": 2, "T": 5}, 
	  "C": {"A": 5, "C": 0, "G": 5, "T": 2}, 
	  "G": {"A": 2, "C": 5, "G": 0, "T": 5}, 
	  "T": {"A": 5, "C": 2, "G": 5, "T": 0}}


score_df = pd.DataFrame({"d": [div for div in d for i in range(len(k)*len(n))],
						"m": [ns for ns in k for i in range(len(n))]*len(d),
						"n": n*len(k)*len(d),
						"n_score": None,
						"ext_score": None,
						"improve": None,
						"n_time": None,
						"ext_time": None})

print("\n######\n### Running Multiple Sequnece Alignment test\n######\n")
for i in range(len(score_df)):
	print(col.GREEN + "Running test " + col.BOLD + str(i+1) + col.END + col.GREEN + " out of " + col.BOLD + str(len(score_df)) + col.END, end="")
	print("\r", end="")
	
	div = score_df.iloc[i, 0]
	# gc = score_df.iloc[i, 0]
	n_seqs = score_df.iloc[i, 1]
	seq_len = score_df.iloc[i, 2]

	# Cost matrix for comput_sp_score function
	# cost = [[0, 5, 2, 5, gc],  # A
	# 		[5, 0, 5, 2, gc],  # C
	# 		[2, 5, 0, 5, gc],  # G
	# 		[5, 2, 5, 0, gc],  # T
	# 		[gc, gc, gc, gc, 0]]  #-'

	n_temp = []
	ext_temp = []
	n_time_temp = []
	ext_time_temp = []
	n_avg = 5
	# Do multiple runs of each instance and average
	for j in range(n_avg):
		# Simulate fasta file with specified parameters
		os.system("python simulate_fasta.py " + str(n_seqs) + " " + str(seq_len) + " " + str(div) + " > test_run.fasta")
		# Read fasta file
		S = fasta_to_lists("/Users/runewind/Documents/Skole/Uni/S9/MSA PiB/test_run.fasta")
		
		### NAIVE
		n_start = time.time()
		# Create naive MSA
		msa_org = MSAligner(S, sm, gc)[0]
		n_end = time.time()
		# Record running time
		n_time_temp.append(n_end - n_start)
		# Convert alignment representation from list of list to list of strings
		msa_org_string = ["".join(l) for l in msa_org]
		# Export alignment to file
		print_alignment_to_file(msa_org_string, "org")
		# Compute score of alignment
		score = compute_sp_score("org_alignment.fasta", cost)
		# save score in dataframe
		n_temp.append(score)
	
		### EXTENSION
		ext_start = time.time()
		# Create extension alignment
		msa_extra = MSA_ext(S, sm, gc)
		ext_end = time.time()
		# Record running time
		ext_time_temp.append(ext_end - ext_start)
		# Convert alignment representation from list of list to list of strings
		msa_e_string = ["".join(l) for l in msa_extra]
		# Export alignment to file
		print_alignment_to_file(msa_e_string, "ext")
		# Compute score of alignment
		score = compute_sp_score("ext_alignment.fasta", cost)
		ext_temp.append(score)

	# Save averaged running times in DataFrame
	score_df.at[i, "n_time"] = sum(n_time_temp) / len(n_time_temp)
	score_df.at[i, "ext_time"] = sum(ext_time_temp) / len(ext_time_temp)
	
	# Save averaged scores in DataFrame
	score_df.at[i, "n_score"] = sum(n_temp) / len(n_temp)
	score_df.at[i, "ext_score"] = sum(ext_temp) / len(ext_temp)
	# Calculate improvement of heuristic
	score_df.at[i, "improve"] = score_df.iloc[i, 3] - score_df.iloc[i, 4]


end_time = time.time()
run_time = end_time - start_time
run_mins = int(run_time // 60)
run_secs = run_time % 60

print(f"\n\nPerformed {len(score_df)} tests averaged over {n_avg} runs in {run_mins} min {run_secs:.2f} sec")
print("\nTest results:")
print(score_df)

score_df.to_csv("test_results.csv", sep="\t", index=False)




