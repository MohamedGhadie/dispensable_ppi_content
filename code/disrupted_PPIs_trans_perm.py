
import os
import pickle
import pandas as pd
import numpy as np
from pathlib import Path
# from interactome_tools import remove_duplicate_PPIs
from protein_function import produce_illumina_expr_dict, coexpr
from stat_tools import fisher_test, sderror_on_fraction

def main():
    
    # reference interactome name
    # options: HI-II-14, IntAct, experiment
    interactome_name = 'HI-II-14'
    
    # method of calculating mutation ∆∆G
    # options: bindprofx, foldx
    ddg_method = 'bindprofx'
    
    # tissue expression database name
    # options: Illumina, GTEx, HPA, Fantom5
    expr_db = 'Illumina'
    
    # minimum number of tissue expression values required for protein pair tissue
    # co-expression to be considered
    coexprMinTissues = 5
    
    # show figures
    showFigs = False
    
    # parent directory of all data files
    dataDir = Path('../data')
    
    # directory of data files from external sources
    extDir = dataDir / 'external'
    
    # parent directory of all processed data files
    procDir = dataDir / 'processed'
    
    # directory of processed data files specific to interactome
    interactomeDir = procDir / interactome_name
    
    # directory of tissue coexpression output data files
    coexprDir = interactomeDir / 'disrupted_PPIs'
        
    # figure directory
    figDir = Path('../figures') / interactome_name / 'disrupted_PPIs'
    
    # input data files
    illuminaExprFile = extDir / 'E-MTAB-513.tsv.txt'
    gtexDir = extDir / 'GTEx_Analysis_v7_eQTL_expression_matrices'
    hpaExprFile = extDir / 'normal_tissue.tsv'
    fantomExprFile = extDir / 'hg38_fair+new_CAGE_peaks_phase1and2_tpm_ann.osc.txt'
    fantomSampleTypeFile = extDir / 'fantom5_sample_type.xlsx'
    uniprotIDmapFile = procDir / 'to_human_uniprotID_map.pkl'
    uniqueGeneSwissProtIDFile = procDir / 'uniprot_unique_gene_reviewed_human_proteome.list'
    if interactome_name is 'experiment':
        neutralPPIFile = interactomeDir / 'neutral_ppis.txt'
        deleteriousPPIFile = interactomeDir / 'deleterious_ppis.txt'
    else:
        neutralPPIFile = interactomeDir / ('neutral_ppis_%s.txt' % ddg_method)
        deleteriousPPIFile = interactomeDir / ('deleterious_ppis_%s.txt' % ddg_method)
    
    # output data files
    proteinExprFile = procDir / ('protein_expr_%s.pkl' % expr_db)
    
    # create directories if not existing
    if not procDir.exists():
        os.makedirs(procDir)
#     if not coexprDir.exists():
#         os.makedirs(coexprDir)
#     if not figDir.exists():
#         os.makedirs(figDir)
    
    #------------------------------------------------------------------------------------
    # load reference and structural interactomes
    #------------------------------------------------------------------------------------
    
#     allneutralPPIs = pd.DataFrame(columns=["Gene_1", "Gene_2", "Protein_1", "Protein_2"])
#     alldeleteriousPPIs = pd.DataFrame(columns=["Gene_1", "Gene_2", "Protein_1", "Protein_2"])
#     for interactome in interactomes:
#         for method in ddg_methods:
#             neutralPPIFile = procDir / interactome / ('neutral_ppis_%s.txt' % method)
#             deleteriousPPIFile = procDir / interactome / ('deleterious_ppis_%s.txt' % method)
#             neutralPPIs = pd.read_table (neutralPPIFile, sep='\t')
#             deleteriousPPIs = pd.read_table (deleteriousPPIFile, sep='\t')
#             neutralPPIs = neutralPPIs[["Gene_1", "Gene_2", "Protein_1", "Protein_2"]]
#             deleteriousPPIs = deleteriousPPIs[["Gene_1", "Gene_2", "Protein_1", "Protein_2"]]
#             allneutralPPIs = allneutralPPIs.append(neutralPPIs, ignore_index=True)
#             alldeleteriousPPIs = alldeleteriousPPIs.append(deleteriousPPIs, ignore_index=True)
#     
#     allneutralPPIs = remove_duplicate_PPIs (allneutralPPIs)
#     alldeleteriousPPIs = remove_duplicate_PPIs (alldeleteriousPPIs)
    
    neutralPPIs = pd.read_table (neutralPPIFile, sep='\t')
    deleteriousPPIs = pd.read_table (deleteriousPPIFile, sep='\t')
    
    print()
    print('PPIs neutral upon disruption: %d (%d proteins)' 
            % (len(neutralPPIs), len(set(neutralPPIs[["Protein_1", "Protein_2"]].values.flatten()))))
    print('PPIs deleterious upon disruption: %d (%d proteins)' 
            % (len(deleteriousPPIs), len(set(deleteriousPPIs[["Protein_1", "Protein_2"]].values.flatten()))))
    
    #------------------------------------------------------------------------------------
    # Produce tissue expression dictionary
    #------------------------------------------------------------------------------------
    
    # produce protein tissue expression profiles
    if not proteinExprFile.is_file():
        print('producing protein tissue expression dictionary')
        if expr_db is 'Illumina':
            produce_illumina_expr_dict (illuminaExprFile,
                                        uniprotIDmapFile,
                                        proteinExprFile,
                                        headers = list(range(1, 18)))
        elif expr_db is 'GTEx':
            produce_gtex_expr_dict (gtexDir,
                                    uniprotIDmapFile,
                                    proteinExprFile,
                                    uniprotIDlistFile = uniqueGeneSwissProtIDFile)
        elif expr_db is 'HPA':
            produce_hpa_expr_dict (hpaExprFile,
                                   uniprotIDmapFile,
                                   proteinExprFile)
        elif expr_db is 'Fantom5':
            produce_fantom5_expr_dict (fantomExprFile,
                                       uniprotIDmapFile,
                                       proteinExprFile,
                                       sampleTypes = 'tissues',
                                       sampleTypeFile = fantomSampleTypeFile,
                                       uniprotIDlistFile = uniqueGeneSwissProtIDFile)
    
    with open(proteinExprFile, 'rb') as f:
        expr = pickle.load(f)
    
    if expr_db is 'HPA':
        exprMap = {'Not detected':0, 'Low':1, 'Medium':2, 'High':3}
        for k, v in expr.items():
            for i, e in enumerate(v):
                v[i] = exprMap[e] if e in exprMap else np.nan
            expr[k] = np.array(v)
    
    #-----------------------------------------------------------------------------------------
    # Calculate tissue co-expression for partners in PPIs that are neutral upon disruption
    #-----------------------------------------------------------------------------------------
       
    neutralPPIs["coexpr"] = neutralPPIs.apply(lambda x: coexpr (x["Protein_1"],
                                                                x["Protein_2"],
                                                                expr,
                                                                minTissues = coexprMinTissues), axis=1)
    neutralPPIs = neutralPPIs [np.isnan(neutralPPIs["coexpr"]) == False].reset_index(drop=True)
    
    #-----------------------------------------------------------------------------------------
    # Calculate tissue co-expression for partners in PPIs that are deleterious upon disruption
    #-----------------------------------------------------------------------------------------
       
    deleteriousPPIs["coexpr"] = deleteriousPPIs.apply(lambda x: coexpr (x["Protein_1"],
                                                                        x["Protein_2"],
                                                                        expr,
                                                                        minTissues = coexprMinTissues), axis=1)
    deleteriousPPIs = deleteriousPPIs [np.isnan(deleteriousPPIs["coexpr"]) == False].reset_index(drop=True)
    
    #------------------------------------------------------------------------------------
    # Calculate the fractions of transient interactions and permanent interactions
    #------------------------------------------------------------------------------------
    
    numNeutralPPIs_transient = sum(neutralPPIs["coexpr"] < 0.1)
    numNeutralPPIs_permanent = sum(neutralPPIs["coexpr"] >= 0.1)
    numDeleteriousPPIs_transient = sum(deleteriousPPIs["coexpr"] < 0.1)
    numDeleteriousPPIs_permanent = sum(deleteriousPPIs["coexpr"] >= 0.1)
    
    # print results
    print()
    print('PPIs neutral upon disruption:')
    print('Fraction of transient PPIs: %f (%d out of %d, SE = %g)' 
            % (numNeutralPPIs_transient / len(neutralPPIs),
               numNeutralPPIs_transient,
               len(neutralPPIs),
               sderror_on_fraction (numNeutralPPIs_transient, len(neutralPPIs))))
    print('Fraction of permanent PPIs: %f (%d out of %d, SE = %g)' 
            % (numNeutralPPIs_permanent / len(neutralPPIs),
               numNeutralPPIs_permanent,
               len(neutralPPIs),
               sderror_on_fraction (numNeutralPPIs_permanent, len(neutralPPIs))))
    print('PPIs deleterious upon disruption:')
    print('Fraction of transient PPIs: %f (%d out of %d, SE = %g)' 
            % (numDeleteriousPPIs_transient / len(deleteriousPPIs),
               numDeleteriousPPIs_transient,
               len(deleteriousPPIs),
               sderror_on_fraction (numDeleteriousPPIs_transient, len(deleteriousPPIs))))
    print('Fraction of permanent PPIs: %f (%d out of %d, SE = %g)' 
            % (numDeleteriousPPIs_permanent / len(deleteriousPPIs),
               numDeleteriousPPIs_permanent,
               len(deleteriousPPIs),
               sderror_on_fraction (numDeleteriousPPIs_permanent, len(deleteriousPPIs))))
    
    print('Statistical significance for difference in transient interaction enrichment')
    fisher_test ([numNeutralPPIs_transient, numNeutralPPIs_permanent],
                 [numDeleteriousPPIs_transient, numDeleteriousPPIs_permanent])

if __name__ == '__main__':
    main()
