#----------------------------------------------------------------------------------------
# This script benchmarks mutation ∆∆G results from FoldX or BindProfX calculations against
# experimental ∆∆G values in SKEMPI.
#
# Run the following script before running this script:
# - process_skempi_file.py
#----------------------------------------------------------------------------------------

import pandas as pd
import numpy as np
from scipy.stats import pearsonr
from pathlib import Path
from ddg_tools import read_protein_mutation_ddg

def main():
    
    # use ∆∆G results for crystal or homology model structures 
    # options: model, crystal
    structure = 'crystal'
    
    # method of calculating mutation ∆∆G
    # options: bindprofx, foldx
    ddg_method = 'bindprofx'
    
    # parent directory of all data files
    dataDir = Path('../data')
    
    # parent directory of all processed data files
    procDir = dataDir / 'processed'
    
    skempiDir = procDir / 'skempi'
    
    # input data files
    predictionFile = skempiDir / ('skempi_%s_mutations_%s_ddg.txt' % (structure, ddg_method))
    experimentFile = skempiDir / 'processed_skempi_mutations.txt'
    
    pred = read_protein_mutation_ddg (predictionFile, type = 'binding')
    exp = pd.read_table (experimentFile, sep='\t')
    
    ddg_pred = []
    for _, row in exp.iterrows():
        k = row.Protein, row.partners, row.Mutation_Position, row.Mut_res
        if k in pred:
            ddg_pred.append(pred[k][-1])
        else:
            ddg_pred.append(np.nan)
    exp["ddg_pred"] = ddg_pred
    exp = exp[(np.isnan(exp["ddg_pred"]) == False) & (np.isnan(exp["ddg"]) == False)]
    
    print('Number of mutations = %d' % len(exp))
    print('Pearson correlation = %f, p-value = %f' % pearsonr(exp["ddg_pred"].values, exp["ddg"].values))
    
if __name__ == '__main__':
    main()
