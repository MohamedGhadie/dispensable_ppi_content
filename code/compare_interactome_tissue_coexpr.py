#----------------------------------------------------------------------------------------
# Compare tissue coexpression of interaction partners between the structural 
# interactome, the reference interactome and random interactions.
#
# Run the following scripts before running this script:
# - produce_data_mappings.py
# - process_interactome.py
# - produce_structural_interactome.py
#----------------------------------------------------------------------------------------

import os
import pickle
import pandas as pd
import numpy as np
from pathlib import Path
from simple_tools import sample_random_pairs
from interactome_tools import read_single_interface_annotated_interactome
from protein_function import (produce_illumina_expr_dict,
                              produce_gtex_expr_dict,
                              produce_hpa_expr_dict,
                              produce_fantom5_expr_dict,
                              coexpr)
from stat_tools import bootstrap_test, sderror
from plot_tools import bar_plot

def main():
    
    # reference interactome name
    # options: HI-II-14, IntAct
    interactome_name = 'HI-II-14'
    
    # tissue expression database name
    # options: Illumina, GTEx, HPA, Fantom5
    expr_db = 'Illumina'
    
    # minimum number of tissue expression values required for protein pair tissue
    # co-expression to be considered
    coexprMinTissues = 5
    
    # number of random interactions
    numRandPairs = 10000
    
    # interactome colors
    interactome_colors = ['limegreen', 'steelblue', 'orangered']
    
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
    coexprDir = interactomeDir / 'coexpr'
        
    # figure directory
    figDir = Path('../figures') / interactome_name / 'coexpr'
    
    # input data files
    illuminaExprFile = extDir / 'E-MTAB-513.tsv.txt'
    gtexDir = extDir / 'GTEx_Analysis_v7_eQTL_expression_matrices'
    hpaExprFile = extDir / 'normal_tissue.tsv'
    fantomExprFile = extDir / 'hg38_fair+new_CAGE_peaks_phase1and2_tpm_ann.osc.txt'
    fantomSampleTypeFile = extDir / 'fantom5_sample_type.xlsx'
    uniprotIDmapFile = procDir / 'to_human_uniprotID_map.pkl'
    uniqueGeneSwissProtIDFile = procDir / 'uniprot_unique_gene_reviewed_human_proteome.list'
    interactomeFile = interactomeDir / 'human_interactome.txt'
    structuralInteractomeFile = interactomeDir / 'human_interface_annotated_interactome.txt'
    
    # output data files
    proteinExprFile = procDir / ('protein_expr_%s.pkl' % expr_db)
    coexprFile = coexprDir / ('interactome_coexpr_%s.pkl' % expr_db)
    
    # create directories if not existing
    if not coexprDir.exists():
        os.makedirs(coexprDir)
    if not figDir.exists():
        os.makedirs(figDir)
    
    #------------------------------------------------------------------------------------
    # load reference and structural interactomes
    #------------------------------------------------------------------------------------
    
    interactome = pd.read_table (interactomeFile)
    print( '\n' + 'reference interactome:' )
    print( '%d PPIs' % len( interactome ) )
    print( '%d proteins' % len( set(interactome[["Protein_1", "Protein_2"]].values.flatten()) ) )
    
    structuralInteractome = read_single_interface_annotated_interactome (structuralInteractomeFile)
    print( '\n' + 'interface-annotated interactome:' )
    print( '%d PPIs' % len( structuralInteractome ) )
    print( '%d proteins' % len( set(structuralInteractome[["Protein_1", "Protein_2"]].values.flatten()) ) )
    
    #------------------------------------------------------------------------------------
    # Produce tissue expression dictionary
    #------------------------------------------------------------------------------------
    
    # produce protein tissue expression profiles
    if not proteinExprFile.is_file():
        print('\n' + 'producing protein tissue expression dictionary')
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
    # Calculate tissue co-expression for random interactions
    #-----------------------------------------------------------------------------------------
    
    proteins = list(set(interactome[["Protein_1", "Protein_2"]].values.flatten()))
    randPairs = pd.DataFrame()
    randPairs["Protein_1"], randPairs["Protein_2"] = zip(* sample_random_pairs (proteins, numRandPairs))    
    randPairs = randPairs [randPairs["Protein_1"] != randPairs["Protein_2"]]
    randPairs["coexpr"] = randPairs.apply(lambda x: coexpr (x["Protein_1"],
                                                            x["Protein_2"],
                                                            expr,
                                                            minTissues = coexprMinTissues), axis=1)
    randPairs = randPairs [np.isnan(randPairs["coexpr"]) == False].reset_index(drop=True)
    
    #-----------------------------------------------------------------------------------------
    # Calculate tissue co-expression for all interaction partners in the reference interactome
    #-----------------------------------------------------------------------------------------
    
    refPPIs = pd.DataFrame()
    refPPIs["Protein_1"] = interactome["Protein_1"].tolist()
    refPPIs["Protein_2"] = interactome["Protein_2"].tolist()    
    refPPIs["coexpr"] = refPPIs.apply(lambda x: coexpr (x["Protein_1"],
                                                        x["Protein_2"],
                                                        expr,
                                                        minTissues = coexprMinTissues), axis=1)
    refPPIs = refPPIs [np.isnan(refPPIs["coexpr"]) == False].reset_index(drop=True)
    
    #------------------------------------------------------------------------------------------
    # Calculate tissue co-expression for all interaction partners in the structural interactome
    #------------------------------------------------------------------------------------------
    
    strucPPIs = pd.DataFrame()
    strucPPIs["Protein_1"] = structuralInteractome["Protein_1"].tolist()
    strucPPIs["Protein_2"] = structuralInteractome["Protein_2"].tolist()    
    strucPPIs["coexpr"] = strucPPIs.apply(lambda x: coexpr (x["Protein_1"],
                                                            x["Protein_2"],
                                                            expr,
                                                            minTissues = coexprMinTissues), axis=1)
    strucPPIs = strucPPIs [np.isnan(strucPPIs["coexpr"]) == False].reset_index(drop=True)
    
    #-------------------------------------------------------------------------------------
    # Save tissue co-expression results to file
    #-------------------------------------------------------------------------------------
    
    allcoexpr = {k:c for k, c in zip(['Random interactions', 'Reference interactome', 'Structural interactome'],
                                     [randPairs, refPPIs, strucPPIs])}
    with open(coexprFile, 'wb') as fout:
        pickle.dump(allcoexpr, fout)
    
    #------------------------------------------------------------------------------------
    # Compare tissue co-expression of interaction partners between reference interactome, 
    # structural interactome and random interactions
    #------------------------------------------------------------------------------------
    
    # remove NaN values
    randCoexpr = randPairs["coexpr"].tolist()
    refCoexpr = refPPIs["coexpr"].tolist()
    strucCoexpr = strucPPIs["coexpr"].tolist()
    
    # print results
    print('\n' + 'Mean tissue co-expression for interaction partners:')
    print('Random pairs: %f (SE = %g, n = %d)' % (np.mean(randCoexpr),
                                                  sderror(randCoexpr),
                                                  len(randCoexpr)))
    print('Reference interactome: %f (SE = %g, n = %d)' % (np.mean(refCoexpr),
                                                           sderror(refCoexpr),
                                                           len(refCoexpr)))
    print('Structural interactome: %f (SE = %g, n = %d)' % (np.mean(strucCoexpr),
                                                            sderror(strucCoexpr),
                                                            len(strucCoexpr)))
    print('\n' + 'Statistical significance')
    print('reference interactome vs random interactions:')
    bootstrap_test(refCoexpr, randCoexpr, iter = 100000)
    print('structural interactome vs reference interactome:')
    bootstrap_test(strucCoexpr, refCoexpr, iter = 100000)
    
    # plot results
    bar_plot([ np.mean(randCoexpr), np.mean(refCoexpr), np.mean(strucCoexpr) ],
             [ sderror(randCoexpr), sderror(refCoexpr), sderror(strucCoexpr) ],
             xlabels = ['Random\ninteractions', 'Reference\ninteractome', 'Structural\ninteractome'],
             ylabel = 'Tissue co-expression of\ninteraction partners',
             colors = interactome_colors,
             edgecolor = 'k',
             barwidth = 0.5,
             fontsize = 24,
             show = showFigs,
             figdir = figDir,
             figname = 'interactome_coexpression_%s' % expr_db)

if __name__ == '__main__':
    main()
