#----------------------------------------------------------------------------------------
# Display properties of the structural interactome.
#
# Run the following scripts before running this script:
# - produce_data_mappings.py
# - produce_structural_interactome.py
#----------------------------------------------------------------------------------------

import os
import pickle
import pandas as pd
import numpy as np
from pathlib import Path
from interactome_tools import read_single_interface_annotated_interactome
from structural_annotation import read_chain_annotated_interactome
from plot_tools import bar_plot, multi_histogram_plot

def main():
    
    # reference interactome name
    # options: HI-II-14, IntAct
    interactome_name = 'IntAct'
    
    # show figures
    showFigs = False
    
    # parent directory of all data files
    dataDir = Path('../data')
    
    # parent directory of all processed data files
    procDir = dataDir / 'processed'
    
    # directory of processed data files specific to interactome
    interactomeDir = procDir / interactome_name
    
    # figure directory
    figDir = Path('../figures') / interactome_name
    
    # input data files
    proteinSeqFile = procDir / 'human_reference_sequences.pkl'
    chainListFile = procDir / 'pdb_seqres_chains.list'
    pdbChainsFile = procDir / 'pdb_seqres_chains.pkl'
    interactomeFile = interactomeDir / 'human_interactome.txt'
    chainAnnotatedInteractomeFile = interactomeDir / 'human_chain_annotated_interactome.txt'
    structuralInteractomeFile = interactomeDir / 'human_interface_annotated_interactome.txt'
    
    # output data files
    ppisPerPDBfile = interactomeDir / 'numPPIsPerPDB.pkl'
    interfaceFile = interactomeDir / 'proteinInterfaces.pkl'
    
    # create directories if not existing
    if not interactomeDir.exists():
        os.makedirs(interactomeDir)
    if not figDir.exists():
        os.makedirs(figDir)
    
    #------------------------------------------------------------------------------------
    # load reference and structural interactomes
    #------------------------------------------------------------------------------------
    
    interactome = pd.read_table (interactomeFile)
    print( '\n' + 'reference interactome:' )
    print( '%d PPIs' % len(interactome) )
    print( '%d proteins' % len(set(interactome[["Protein_1", "Protein_2"]].values.flatten())) )
    
    with open(proteinSeqFile, 'rb') as f:
        proteinSeq = pickle.load(f)
    print( '\n' + 'Protein sequences:' )
    print( '%d sequences' % len(proteinSeq.keys()) )
    
    with open(pdbChainsFile, 'rb') as f:
        pdbChains = pickle.load(f)
    with open(chainListFile, 'r') as f:
        chainIDs = set(f.read().split())
    print( '\n' + 'PDB structures available:' )
    print( '%d structures' % len(pdbChains.keys()) )
    print( '%d chains' % len(chainIDs) )
    
    chainAnnotatedInteractome = read_chain_annotated_interactome (chainAnnotatedInteractomeFile)
    print( '\n' + 'chain-annotated interactome:' )
    print( '%d PPIs' % len(chainAnnotatedInteractome) )
    print( '%d proteins' % len(set(chainAnnotatedInteractome[["Protein_1", "Protein_2"]].values.flatten())) )
    
    structuralInteractome = read_single_interface_annotated_interactome (structuralInteractomeFile)
    print( '\n' + 'interface-annotated interactome:' )
    print( '%d PPIs' % len(structuralInteractome) )
    print( '%d proteins' % len(set(structuralInteractome[["Protein_1", "Protein_2"]].values.flatten())) )

    #------------------------------------------------------------------------------------
    # Calculate distribution of number of PPIs modeled per PDB model in the
    # structural interactome
    #------------------------------------------------------------------------------------
    
    ppisPerPDB = {}
    for _, row in structuralInteractome.iterrows():
        pdbIDs = {c1.split('_')[0] for c1, _ in row.Chain_pairs}
        for pdbID in pdbIDs:
            if pdbID in ppisPerPDB:
                ppisPerPDB[pdbID].add(row.Protein_1 + '_' + row.Protein_2)
            else:
                ppisPerPDB[pdbID] = {row.Protein_1 + '_' + row.Protein_2}
    
    pdbIDs = list(ppisPerPDB.keys())
    num_ppisPerPDB = [len(ppisPerPDB[id]) for id in pdbIDs]
    with open(ppisPerPDBfile, 'wb') as fout:
        pickle.dump(num_ppisPerPDB, fout)
    
    maxVal = 10
    num_ppisPerPDB_toplot = [min(n, maxVal) for n in num_ppisPerPDB]
    num_ppisPerPDB_toplot = [num_ppisPerPDB_toplot.count(i) for i in np.arange(1, maxVal + 1)]
    
    print('\n' + 'Number of PDB models unique to 1 PPI: %d out of %d' 
          % (num_ppisPerPDB_toplot[0], sum(num_ppisPerPDB_toplot)))
    
    bar_plot(num_ppisPerPDB_toplot,
             xlabels = list(map(str, np.arange(1, maxVal))) + [ str(maxVal) + '+'],
             xlabel = 'Number of PPIs modeled by PDB structure',
             ylabel = 'Number of PDB structures',
             colors = 'red',
             fontsize = 24,
             show = showFigs,
             figdir = figDir,
             figname = 'numPPIsPerPDB')

    #------------------------------------------------------------------------------------
    # Calculate distributions for protein interface length and interface fraction in
    # the structural interactome
    #------------------------------------------------------------------------------------
         
    # compile interfaces for all interacting proteins
    interfaces = {}
    for _, row in structuralInteractome.iterrows():
        interface_1, interface_2 = row.Interfaces
        interfaces[(row.Protein_1, row.Protein_2)] = interface_1
        interfaces[(row.Protein_2, row.Protein_1)] = interface_2
    with open(interfaceFile, 'wb') as fout:
        pickle.dump(interfaces, fout)
    
    # calculate number of interface residues per interacting protein
    interfaceLength = [len(i) for i in interfaces.values()]
    print('Mean interface length = %.1f residues' % np.mean(interfaceLength))
    multi_histogram_plot (interfaceLength,
                          colors = 'b',
                          xlabel = 'Number of interfacial residues',
                          ylabel = 'Number of proteins',
                          edgecolor = 'k',
                          fontsize = 24,
                          bins = 20,
                          alpha = 1,
                          show = showFigs,
                          figdir = figDir,
                          figname = 'interface_length')
    
    # calculate ratio of interface residues per interacting protein
    with open(proteinSeqFile, 'rb') as f:
        proteinSeq = pickle.load(f)
    interfaceRatio = [len(i)/len(proteinSeq[p]) for (p, _), i in interfaces.items()]
    print('Mean interface fraction of protein sequence = %.1f%%' % (100 * np.mean(interfaceRatio)))
    multi_histogram_plot (interfaceRatio,
                          colors = 'g',
                          xlabel = 'Fraction of interfacial residues',
                          ylabel = 'Number of proteins',
                          edgecolor = 'k',
                          fontsize = 24,
                          bins = 20,
                          alpha = 1,
                          show = showFigs,
                          figdir = figDir,
                          figname = 'interface_fraction')

if __name__ == "__main__":
    main()
