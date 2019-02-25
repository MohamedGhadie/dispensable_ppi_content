#----------------------------------------------------------------------------------------
# This script compares functional similarity of interaction partners between the reference 
# interactome and the structural interactome.
#----------------------------------------------------------------------------------------

import os
import pickle
import pandas as pd
import numpy as np
from pathlib import Path
from interactome_tools import read_single_interface_annotated_interactome
from protein_function import (produce_illumina_expr_dict,
                              produce_gtex_expr_dict,
                              coexpr)
from stat_tools import bootstrap_test, sderror
from plot_tools import bar_plot

def main():
    
    # reference interactome name. Options: 'HI-II-14' or 'IntAct'
    interactome_name = 'HI-II-14'
    
    # tissue expression database name. Options: 'Illumina', 'GTEx'
    expr_db = 'GTEx'
    
    # minimum number of tissue expression values required for protein pair tissue
    # co-expression to be considered
    coexprMinTissues = 5
    
    # show figures
    showFigs = False
    
    # directory of processed data files shared by all interactomes
    dataDir = Path('../data')
    
    # directory of data files from external sources
    extDir = dataDir / 'external'
    
    # directory of processed data files shared by all interactomes
    inDir = dataDir / 'processed'
    
    # directory to processed data files specific to interactome
    interactomeDir = dataDir / 'processed' / interactome_name
    
    # directory to save processed data shared by all interactomes
    outDir = dataDir / 'processed'
    
    # figure directory
    figDir = Path('../figures') / interactome_name 
    
    # create directories if not existing
    if not interactomeDir.exists():
        os.makedirs(interactomeDir)
    if not outDir.exists():
        os.makedirs(outDir)
    if not figDir.exists():
        os.makedirs(figDir)
    
    # input data files
    illuminaExprFile = extDir / 'E-MTAB-513.tsv.txt'
    gtexDir = extDir / 'GTEx_Analysis_v7_eQTL_expression_matrices'
    UniprotIDmapFile = inDir / 'to_human_uniprotID_map.pkl'
    interactomeFile = interactomeDir / 'human_interactome.txt'
    structuralInteractomeFile = interactomeDir / 'human_interface_annotated_interactome.txt'
    
    # output data files
    proteinExprFile = outDir / ('protein_expr_%s.pkl' % expr_db)
    coexprFile = interactomeDir / ('interactome_coexpr_%s.pkl' % expr_db)
    
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
    # Produce GO and tissue expression dictionaries
    #------------------------------------------------------------------------------------
    
    # produce protein tissue expression profiles
    if not proteinExprFile.is_file():
        print('\n' + 'producing protein tissue expression dictionary')
        if expr_db is 'Illumina':
            produce_illumina_expr_dict (illuminaExprFile,
                                        UniprotIDmapFile,
                                        proteinExprFile,
                                        headers = list(range(1, 18)))
        elif expr_db is 'GTEx':
            produce_gtex_expr_dict (gtexDir,
                                    UniprotIDmapFile,
                                    proteinExprFile)
    
    with open(proteinExprFile, 'rb') as f:
        expr = pickle.load(f)
    
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
    
    # save tissue co-expression results for use by other scripts
    with open(coexprFile, 'wb') as fout:
        pickle.dump([refPPIs, strucPPIs], fout)
    
    #------------------------------------------------------------------------------------
    # Compare tissue co-expression of interaction partners between reference interactome 
    # and structural interactome
    #------------------------------------------------------------------------------------
    
    # remove NaN values
    refCoexpr = refPPIs["coexpr"][np.isnan(refPPIs["coexpr"]) == False].tolist()
    strucCoexpr = strucPPIs["coexpr"][np.isnan(strucPPIs["coexpr"]) == False].tolist()
    
    # print results
    print('\n' + 'Mean tissue co-expression for interaction partners:')
    print('Reference interactome: %f (SE = %g, n = %d)' % (np.mean(refCoexpr),
                                                           sderror(refCoexpr),
                                                           len(refCoexpr)))
    print('Structural interactome: %f (SE = %g, n = %d)' % (np.mean(strucCoexpr),
                                                            sderror(strucCoexpr),
                                                            len(strucCoexpr)))
    print( 'Calculating statistical significance using bootstrap test' )
    bootstrap_test(refCoexpr, strucCoexpr, iter = 10000)
    
    # plot results
    bar_plot([ np.mean(refCoexpr), np.mean(strucCoexpr) ],
             [ sderror(refCoexpr), sderror(strucCoexpr) ],
             xlabels = ['Reference\ninteractome', 'Structural\ninteractome'],
             ylabel = 'Tissue co-expression of\ninteracting protein pairs',
             colors = ['blue', 'red'],
             edgecolor = 'k',
             barwidth = 0.5,
             fontsize = 24,
             show = showFigs,
             figdir = figDir,
             figname = 'interactome_coexpression_%s' % expr_db)

if __name__ == '__main__':
    main()
