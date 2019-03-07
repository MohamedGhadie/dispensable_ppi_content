#----------------------------------------------------------------------------------------
# Compare functional similarity of interaction partners between the reference 
# interactome and the structural interactome. Similarity is calculated as the 
# jaccard index between gene ontology (GO) terms of interaction partners.
#
# Run the following scripts before running this script:
# - process_interactome.py
# - produce_structural_interactome.py
#----------------------------------------------------------------------------------------

import os
import pickle
import pandas as pd
import numpy as np
from pathlib import Path
from interactome_tools import read_single_interface_annotated_interactome
from protein_function import produce_protein_go_dictionaries, go_sim
from stat_tools import bootstrap_test, sderror
from plot_tools import bar_plot

def main():
    
    # reference interactome name. Options: 'HI-II-14' or 'IntAct'
    interactome_name = 'IntAct'
    
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
    
    # figure directory
    figDir = Path('../figures') / interactome_name 
    
    # create directories if not existing
    if not procDir.exists():
        os.makedirs(procDir)
    if not interactomeDir.exists():
        os.makedirs(interactomeDir)
    if not figDir.exists():
        os.makedirs(figDir)
    
    # input data files
    extGOfile = extDir / 'goa_human.gaf'
    interactomeFile = interactomeDir / 'human_interactome.txt'
    structuralInteractomeFile = interactomeDir / 'human_interface_annotated_interactome.txt'
    
    # output data files
    GOfile = procDir / 'go.pkl'
    MFfile = procDir / 'gof.pkl'
    BPfile = procDir / 'gop.pkl'
    CCfile = procDir / 'goc.pkl'
    GOSimFile = interactomeDir / 'goSim_reference_structural.pkl'
    
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
    # Produce GO dictionaries
    #------------------------------------------------------------------------------------
    
    # produce protein GO association profiles
    if not GOfile.is_file():
        print( 'producing protein GO dictionaries' )
        produce_protein_go_dictionaries(extGOfile,
                                        GOfile,
                                        MFfile,
                                        BPfile,
                                        CCfile)
    
    with open(GOfile, 'rb') as f:
        GOassoc = pickle.load(f)
    with open(MFfile, 'rb') as f:
        MFassoc = pickle.load(f)
    with open(BPfile, 'rb') as f:
        BPassoc = pickle.load(f)
    with open(CCfile, 'rb') as f:
        CCassoc = pickle.load(f)
    
    # print number of GO terms for each category
    print('\n' + '%d gene ontology (GO) terms' % len({v for val in GOassoc.values() for v in val}))
    print('%d molecular function (MF) terms' % len({v for val in MFassoc.values() for v in val}))
    print('%d biological process (BP) terms' % len({v for val in BPassoc.values() for v in val}))
    print('%d cellular component (CC) terms' % len({v for val in CCassoc.values() for v in val}))
    
    #------------------------------------------------------------------------------------
    # Calculate GO similarity for all interaction partners in the reference interactome
    #------------------------------------------------------------------------------------
    
    refPPIs = pd.DataFrame()
    refPPIs["Protein_1"] = interactome["Protein_1"].tolist()
    refPPIs["Protein_2"] = interactome["Protein_2"].tolist()    
    refPPIs["goSim"] = refPPIs.apply(lambda x: go_sim(x["Protein_1"], x["Protein_2"], GOassoc), axis=1)
    refPPIs["mfSim"] = refPPIs.apply(lambda x: go_sim(x["Protein_1"], x["Protein_2"], MFassoc), axis=1)
    refPPIs["bpSim"] = refPPIs.apply(lambda x: go_sim(x["Protein_1"], x["Protein_2"], BPassoc), axis=1)
    refPPIs["ccSim"] = refPPIs.apply(lambda x: go_sim(x["Protein_1"], x["Protein_2"], CCassoc), axis=1)
    
    #------------------------------------------------------------------------------------
    # Calculate GO similarity for all interaction partners in the structural interactome
    #------------------------------------------------------------------------------------
    
    strucPPIs = pd.DataFrame()
    strucPPIs["Protein_1"] = structuralInteractome["Protein_1"].tolist()
    strucPPIs["Protein_2"] = structuralInteractome["Protein_2"].tolist()    
    strucPPIs["goSim"] = strucPPIs.apply(lambda x: go_sim(x["Protein_1"], x["Protein_2"], GOassoc), axis=1)
    strucPPIs["mfSim"] = strucPPIs.apply(lambda x: go_sim(x["Protein_1"], x["Protein_2"], MFassoc), axis=1)
    strucPPIs["bpSim"] = strucPPIs.apply(lambda x: go_sim(x["Protein_1"], x["Protein_2"], BPassoc), axis=1)
    strucPPIs["ccSim"] = strucPPIs.apply(lambda x: go_sim(x["Protein_1"], x["Protein_2"], CCassoc), axis=1)
    
    # save functional similarity results for use by other scripts
    with open(GOSimFile, 'wb') as fout:
        pickle.dump([refPPIs, strucPPIs], fout)
    
    #------------------------------------------------------------------------------------
    # Compare GO similarity of interaction partners between reference interactome and
    # structural interactome
    #------------------------------------------------------------------------------------
    
    # remove NaN values
    refGOsim = refPPIs["goSim"][np.isnan(refPPIs["goSim"]) == False].tolist()
    strucGOsim = strucPPIs["goSim"][np.isnan(strucPPIs["goSim"]) == False].tolist()
    
    # print results
    print('\nMean GO similarity for interaction partners:')
    print('Reference interactome: %f (SE = %g, n = %d)' % (np.mean(refGOsim),
                                                           sderror(refGOsim),
                                                           len(refGOsim)))
    print('Structural interactome: %f (SE = %g, n = %d)' % (np.mean(strucGOsim),
                                                            sderror(strucGOsim),
                                                            len(strucGOsim)))
    print( 'Calculating statistical significance using bootstrap test' )
    bootstrap_test(refGOsim, strucGOsim, iter = 10000)
    
    # plot results
    bar_plot([ np.mean(refGOsim), np.mean(strucGOsim) ],
             [ sderror(refGOsim), sderror(strucGOsim) ],
             xlabels = ['Reference\ninteractome', 'Structural\ninteractome'],
             ylabel = 'Functional similarity of\ninteracting protein pairs',
             colors = ['blue', 'red'],
             edgecolor = 'k',
             barwidth = 0.5,
             fontsize = 24,
             show = showFigs,
             figdir = figDir,
             figname = 'interactome_functional_similarity')

if __name__ == '__main__':
    main()
