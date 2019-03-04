#----------------------------------------------------------------------------------------
# This script compares functional similarity of interaction partners between the 
# structural interactome, the reference interactome and random interactions.
#----------------------------------------------------------------------------------------

import os
import pickle
import pandas as pd
import numpy as np
from pathlib import Path
from simple_tools import sample_random_pairs
from text_tools import produce_item_list
from interactome_tools import read_single_interface_annotated_interactome
from protein_function import produce_fastsemsim_protein_gosim_dict
from stat_tools import bootstrap_test, sderror
from plot_tools import bar_plot

def main():
    
    # reference interactome name. Options: HI-II-14, IntAct
    interactome_name = 'IntAct'
    
    # similarity measure to calculate GO similarity
    # options: Resnik, Lin, Jiang-Conrath, SimGIC, SimUI, SimIC, SimRel, Dice, SimTO
    #           SimNTO, Jaccard, Czekanowski-Dice, Cosine, GSESAME, SimICND, SimICNP
    sim_measure = 'Resnik'
    
    # mixing strategy for merging GO term semantic similarities
    # options: max, avg, BMA (best match average)
    mix_method = 'BMA'
    
    # root ontology on which GO similarity is calculated
    # options: biological_process, molecular_function, cellular_component
    ont_root = 'cellular_component'
    
    # root ontology labels used to label output files and figures
    ont_abv = {'biological_process':'bp', 'molecular_function':'mf', 'cellular_component':'cc'}
    
    # list of ontological relationships to ignore
    ont_ignore = None
    
    # list of evidence codes to ignore
    ec_ignore = None
    
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
    
    # directory of processed data files shared by all interactomes
    inDir = dataDir / 'processed'
    
    # directory to processed data files specific to interactome
    interactomeDir = dataDir / 'processed' / interactome_name
    
    # directory to save output data files
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
    
    # get root ontology label to label output files and figures
    ont_label = ont_abv [ont_root]
    
    # input data files
    ontologyFile = extDir / 'go-basic.obo'
    annotationFile = extDir / 'goa_human.gaf'
    interactomeFile = interactomeDir / 'human_interactome.txt'
    structuralInteractomeFile = interactomeDir / 'human_interface_annotated_interactome.txt'
    
    # output data files
    refPPIListFile = interactomeDir / 'refPPIs.txt'
    randPairListFile = interactomeDir / 'randPairs.txt'
    refPPIgosimParamFile = interactomeDir / ('fastsemsim_parameters_refPPIs_%s_%s' % (ont_label, sim_measure))
    randPairgosimParamFile = interactomeDir / ('fastsemsim_parameters_randPairs_%s_%s' % (ont_label, sim_measure))
    refPPIfastsemsimOutFile = interactomeDir / ('fastsemsim_output_refPPIs_%s_%s' % (ont_label, sim_measure))
    randPairfastsemsimOutFile = interactomeDir / ('fastsemsim_output_randPairs_%s_%s' % (ont_label, sim_measure))
    refPPIgosimFile = interactomeDir / ('gosim_refPPIs_%s_%s.pkl' % (ont_label, sim_measure))
    randPairgosimFile = interactomeDir / ('gosim_randPairs_%s_%s.pkl' % (ont_label, sim_measure))
    allPPIgosimFile = interactomeDir / ('allPPI_gosim_%s_%s.pkl' % (ont_label, sim_measure))
                                                   
    #------------------------------------------------------------------------------------
    # load reference and structural interactomes
    #------------------------------------------------------------------------------------
    
    interactome = pd.read_table (interactomeFile)
    interactomeProteins = list(set(interactome[["Protein_1", "Protein_2"]].values.flatten()))
    print( '\n' + 'reference interactome:' )
    print( '%d PPIs' % len(interactome) )
    print( '%d proteins' % len(interactomeProteins) )
    
    structuralInteractome = read_single_interface_annotated_interactome (structuralInteractomeFile)
    strucInteractomeProteins = list(set(structuralInteractome[["Protein_1", "Protein_2"]].values.flatten()))
    print( '\n' + 'interface-annotated interactome:' )
    print( '%d PPIs' % len(structuralInteractome) )
    print( '%d proteins' % len(strucInteractomeProteins) )
    
    #------------------------------------------------------------------------------------
    # Calculate GO similarity for random interactions
    #------------------------------------------------------------------------------------
    
    randPairs = pd.DataFrame()
    randPairs["Protein_1"], randPairs["Protein_2"] = zip(* sample_random_pairs (interactomeProteins, numRandPairs))    
    randPairs = randPairs [randPairs["Protein_1"] != randPairs["Protein_2"]]
    
    # produce protein GO similarity dictionary
    print('\n' + 'producing GO similarity dictionary for random pairs')
    randPairs[["Protein_1","Protein_2"]].to_csv(randPairListFile, index=False, header=False, sep='\t')
    produce_fastsemsim_protein_gosim_dict (randPairListFile,
                                           randPairgosimFile,
                                           sim_measure = sim_measure,
                                           mix_method = mix_method,
                                           ont_root = ont_root,
                                           ont_ignore = ont_ignore,
                                           ec_ignore = ec_ignore,
                                           ontologyFile = ontologyFile,
                                           annotationFile = annotationFile,
                                           paramOutFile = randPairgosimParamFile,
                                           fastsemsimOutFile = randPairfastsemsimOutFile)
    with open(randPairgosimFile, 'rb') as f:
        gosim = pickle.load(f)
    
    sim = []
    for p in map(tuple, map(sorted, randPairs[["Protein_1","Protein_2"]].values)):
        sim.append(gosim[p] if p in gosim else np.nan)
    randPairs["gosim"] = sim
    randPairs = randPairs [np.isnan(randPairs["gosim"]) == False].reset_index(drop=True)
    
    #------------------------------------------------------------------------------------
    # Calculate GO similarity for all interaction partners in the reference interactome
    #------------------------------------------------------------------------------------
    
    refPPIs = pd.DataFrame()
    refPPIs["Protein_1"], refPPIs["Protein_2"] = zip(* interactome[["Protein_1","Protein_2"]].values)
    
    # produce protein GO similarity dictionary
    if not refPPIgosimFile.is_file():
        print('\n' + 'producing GO similarity dictionary for reference PPIs')
        refPPIs[["Protein_1","Protein_2"]].to_csv(refPPIListFile, index=False, header=False, sep='\t')
        produce_fastsemsim_protein_gosim_dict (refPPIListFile,
                                               refPPIgosimFile,
                                               sim_measure = sim_measure,
                                               mix_method = mix_method,
                                               ont_root = ont_root,
                                               ont_ignore = ont_ignore,
                                               ec_ignore = ec_ignore,
                                               ontologyFile = ontologyFile,
                                               annotationFile = annotationFile,
                                               paramOutFile = refPPIgosimParamFile,
                                               fastsemsimOutFile = refPPIfastsemsimOutFile)
    with open(refPPIgosimFile, 'rb') as f:
        gosim = pickle.load(f)
        
    sim = []
    for p in map(tuple, map(sorted, refPPIs[["Protein_1","Protein_2"]].values)):
        sim.append(gosim[p] if p in gosim else np.nan)
    refPPIs["gosim"] = sim
    refPPIs = refPPIs [np.isnan(refPPIs["gosim"]) == False].reset_index(drop=True)
    
    #------------------------------------------------------------------------------------
    # Calculate GO similarity for all interaction partners in the structural interactome
    #------------------------------------------------------------------------------------
    
    strucPPIs = pd.DataFrame()
    strucPPIs["Protein_1"], strucPPIs["Protein_2"] = zip(* structuralInteractome[["Protein_1","Protein_2"]].values)   
    sim = []
    for p in map(tuple, map(sorted, strucPPIs[["Protein_1","Protein_2"]].values)):
        sim.append(gosim[p] if p in gosim else np.nan)
    strucPPIs["gosim"] = sim
    strucPPIs = strucPPIs [np.isnan(strucPPIs["gosim"]) == False].reset_index(drop=True)
    
    #------------------------------------------------------------------------------------
    # Save GO similarity results to file
    #------------------------------------------------------------------------------------
    
    allgosim = {k:sim for k, sim in zip(['Random interactions', 'Reference interactome', 'Structural interactome'],
                                        [randPairs, refPPIs, strucPPIs])}
    with open(allPPIgosimFile, 'wb') as fout:
        pickle.dump(allgosim, fout)
    
    #------------------------------------------------------------------------------------
    # Compare GO similarity of interaction partners between reference interactome, 
    # structural interactome and random interactions
    #------------------------------------------------------------------------------------
    
    # remove NaN values
    randGOsim = randPairs["gosim"].tolist()
    refGOsim = refPPIs["gosim"].tolist()
    strucGOsim = strucPPIs["gosim"].tolist()
    
    # print results
    print('\n' + 'Mean %s similarity for interaction partners:' % ont_root)
    print('Random interactions: %f (SE = %g, n = %d)' % (np.mean(randGOsim),
                                                         sderror(randGOsim),
                                                         len(randGOsim)))
    print('Reference interactome: %f (SE = %g, n = %d)' % (np.mean(refGOsim),
                                                           sderror(refGOsim),
                                                           len(refGOsim)))
    print('Structural interactome: %f (SE = %g, n = %d)' % (np.mean(strucGOsim),
                                                            sderror(strucGOsim),
                                                            len(strucGOsim)))
    print('\n' + 'Statistical significance')
    print('reference interactome vs random interactions:')
    bootstrap_test(refGOsim, randGOsim, iter = 10000)
    print('structural interactome vs reference interactome:')
    bootstrap_test(strucGOsim, refGOsim, iter = 10000)
    
    # plot results
    bar_plot([ np.mean(randGOsim), np.mean(refGOsim), np.mean(strucGOsim) ],
             [ sderror(randGOsim), sderror(refGOsim), sderror(strucGOsim) ],
             xlabels = ['Random\ninteractions', 'Reference\ninteractome', 'Structural\ninteractome'],
             ylabel = '%s similarity of\ninteraction partners' % ont_root,
             colors = interactome_colors,
             edgecolor = 'k',
             barwidth = 0.5,
             fontsize = 24,
             show = showFigs,
             figdir = figDir,
             figname = 'interactome_gosim_%s_%s' % (ont_label, sim_measure))

if __name__ == '__main__':
    main()
