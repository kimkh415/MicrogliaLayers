import sys
import os
import pandas as pd


def select_interactions(means_path, pvalues_path, out_name, p_thresh=0.05, m_thresh=0.5):
	
	all_pval = pd.read_csv(pvalues_path, sep='\t', header=0)
	all_means = pd.read_csv(means_path, sep='\t', header=0)

	intr = all_pval["interacting_pair"]

	pval = all_pval.iloc[:, 11:]
	means = all_means.iloc[:, 11:]

	min_pvalues = pval.min(axis=1)
	max_means = means.max(axis=1)

	selected_rows = pick_rows(min_pvalues, max_means, intr, p_thresh, m_thresh)
	#print(selected_rows)
	
	output_rowcol(selected_rows, out_name)
	

def pick_rows(min_pvals, max_means, interaction_pairs, p_thresh, m_thresh):
	
	picked = []

	for i in range(len(interaction_pairs)):
		if min_pvals[i] < p_thresh and max_means[i] > m_thresh:
			#print(interaction_pairs[i], min_pvals[i], max_means[i])
			picked.append(interaction_pairs[i])
	
	return picked


def output_rowcol(names, outname):
	out = open(outname, 'w')
	rows = '\n'.join(names)
	out.write(rows)
	out.close()


if __name__ == "__main__":
	file_path = sys.argv[1]
	#means_file = sys.argv[1]
	#pvalues_file = sys.argv[2]

	mf = os.path.join(file_path, "filtered_means.txt")
	mf_r = os.path.join(file_path, "filtered_means_r.txt")
	pf = os.path.join(file_path, "filtered_padj.txt")
	pf_r = os.path.join(file_path, "filtered_padj_r.txt")

	select_interactions(mf, pf, "rows.txt")
	select_interactions(mf_r, pf_r, "rows_r.txt")
