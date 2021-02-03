"""
Script to separate CPDB output matrices into two for visualization.

eg) 
means.txt
			clusterA_clusterB	clusterB_clusterA
	interaction1			1			4
	interaction2			2			5
	interaction3			3			6


filtered_means.txt
			clusterA_clusterB
	interaction1			1
	interaction2			2
	interaction3			3


filtered_means_r.txt
			clusterB_clusterA
	interaction1			4
	interaction2			5
	interaction3			6

"""


import sys
import pandas as pd
import os

# Path to the CPDB output directory
path = sys.argv[1]

# All neuronal clusters
neuronal = ["L4_Stellate__Rorb+", "L2/3_Ipsi__Trpc6+", "L2/3_CPN__Slc39a12+", "L2/3_CPN__Cdh13+", "L6_CThPN__Hcrtr2+", "KO_Mismatch_1", "KO_Mismatch_2", "KO_Mismatch_3", "KO_Mismatch_4", "L2/3_CPN__Qrfpr+", "L2/3_CPN__Inhba+", "KO_Mismatch_5", "L6_CThPN__Cpa6+", "L5a_CStrPN__Deptor+", "L5_CPN__Dock5+", "KO_Mismatch_6", "L5_NP__Trhr+", "L5_PT__Lgr5+", "L6b_SP__Ctgf+"]

# All microglia clusters
microglia = ["Mg_Hom_1", "Mg_Hom_2", "Mg_Hom_1_KO", "Mg_Hom_2_KO"]

# load CPDB output
mean_file = os.path.join(path, "means.txt")
pval_file = os.path.join(path, "pvalues.txt")
m = pd.read_csv(mean_file, sep='\t', header=0)
pval = pd.read_csv(pval_file, sep='\t', header=0)

# means and pvalues output matrix share the same column names
cols = m.columns.tolist()

prefix = cols[:11]  # metadata columns

cols = cols[12:]  # cluster pair columns
pairs = [x for x in cols if x.split('|')[0] in neuronal]
pairs = [x for x in pairs if x.split('|')[1] in microglia]

pairs_r = ["{}|{}".format(x.split('|')[1], x.split('|')[0]) for x in pairs]

pairs = prefix + pairs
pairs_r = prefix + pairs_r

out_m = os.path.join(path, "filtered_means.txt")
out_m_r = os.path.join(path, "filtered_means_r.txt")
out_p = os.path.join(path, "filtered_pvalues.txt")
out_p_r = os.path.join(path, "filtered_pvalues_r.txt")

m.to_csv(out_m, sep='\t', columns=pairs, index=False)
m.to_csv(out_m_r, sep='\t', columns=pairs_r, index=False)

pval.to_csv(out_p, sep='\t', columns=pairs, index=False)
pval.to_csv(out_p_r, sep='\t', columns=pairs_r, index=False)

print("DONE")
