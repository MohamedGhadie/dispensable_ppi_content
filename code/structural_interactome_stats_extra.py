#----------------------------------------------------------------------------------------
# This script loads the reference interactome, chain-annotated interactome, and structural 
# interactome, then compares functional similarity and tissue co-expression of interaction
# partners between the reference interactome and structural interactome. Other statistics 
# on the models used and interfaces mapped in the structural interactome are presented.
#----------------------------------------------------------------------------------------

import os
import pickle
import pandas as pd
import numpy as np
from pathlib import Path
from id_mapping import (produce_uniqueGene_swissProtIDs,
                        produce_proteinSeq_dict,
                        produce_uniprotID_dict,
                        produce_geneName_dict,
                        produce_chain_dict)
from interactome_tools import read_single_interface_annotated_interactome
from structural_annotation import read_chain_annotated_interactome
from pdb_tools import produce_chain_list
from protein_function import (produce_protein_go_dictionaries,
                              produce_protein_expr_dict,
                              go_sim,
                              coexpr)
from stat_tools import bootstrap_test, sderror
from plot_tools import bar_plot, multi_histogram_plot, heatmap_plot

def main():
    
    # reference interactome name. Options: 'HI-II-14' or 'IntAct'
    interactome_name = 'IntAct'
    
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
    GeneMapFile = inDir / 'to_human_geneName_map.pkl'
    UniqueGeneSwissProtIDFile = inDir / 'uniprot_unique_gene_reviewed_human_proteome.list'
    ProteinSeqFile = inDir / 'human_reference_sequences.pkl'
    UniprotIDmapFile = inDir / 'to_human_uniprotID_map.pkl'
    chainListFile = inDir / 'pdb_seqres_chains.list'
    pdbChainsFile = inDir / 'pdb_seqres_chains.pkl'
    interactomeFile = interactomeDir / 'human_interactome.txt'
    chainAnnotatedInteractomeFile = interactomeOutDir / 'human_chain_annotated_interactome.txt'
    structuralInteractomeFile = interactomeDir / 'human_interface_annotated_interactome.txt'
    
    # output data files
    
    
    #------------------------------------------------------------------------------------
    # load reference and structural interactomes
    #------------------------------------------------------------------------------------
    
    interactome = pd.read_table (interactomeFile)
    print( '\n' + 'reference interactome:' )
    print( '%d PPIs' % len( interactome ) )
    print( '%d proteins' % len( set(interactome[["Protein_1", "Protein_2"]].values.flatten()) ) )
    
    with open(ProteinSeqFile, 'rb') as f:
        proteinSeq = pickle.load(f)
    print( '\n' + 'Protein sequences:' )
    print( '%d sequences' % len( proteinSeq.keys() ) )
    
    with open(pdbChainsFile, 'rb') as f:
        pdbChains = pickle.load(f)
    with open(chainListFile, 'r') as f:
        chainIDs = set(f.read().split())
    print( '\n' + 'PDB structures available:' )
    print( '%d structures' % len( pdbChains.keys() ) )
    print( '%d chains' % len( chainIDs ) )
    
    chainAnnotatedInteractome = read_chain_annotated_interactome (chainAnnotatedInteractomeFile)
    print( '\n' + 'chain-annotated interactome:' )
    print( '%d PPIs' % len( chainAnnotatedInteractome ) )
    print( '%d proteins' % len( set(chainAnnotatedInteractome[["Protein_1", "Protein_2"]].values.flatten()) ) )
    
    structuralInteractome = read_single_interface_annotated_interactome (structuralInteractomeFile)
    print( '\n' + 'interface-annotated interactome:' )
    print( '%d PPIs' % len( structuralInteractome ) )
    print( '%d proteins' % len( set(structuralInteractome[["Protein_1", "Protein_2"]].values.flatten()) ) )

    #------------------------------------------------------------------------------------
    # Calculate distribution of number of PPIs modeled per PDB model in the
    # structural interactome
    #------------------------------------------------------------------------------------
    
    ppisPerPDB = {}
    mfSimPerPDB = {}
    bpSimPerPDB = {}
    ccSimPerPDB = {}
    coexprPerPDB = {}
    for i, row in structuralInteractome.iterrows():
        pdbIDs = set( [ c1.split('_')[0] for c1, _ in row.Chain_pairs] )
        for pdbID in pdbIDs:
            if pdbID in ppisPerPDB:
                ppisPerPDB[pdbID].add( row.Protein_1 + '_' + row.Protein_2 )
            else:
                ppisPerPDB[pdbID] = { row.Protein_1 + '_' + row.Protein_2 }
            
            if not np.isnan(strucPPIs.loc[i, "mfSim"]):
                if pdbID in mfSimPerPDB:
                    mfSimPerPDB[pdbID].append( strucPPIs.loc[i, "mfSim"] )
                else:
                    mfSimPerPDB[pdbID] = [ strucPPIs.loc[i, "mfSim"] ]
            
            if not np.isnan(strucPPIs.loc[i, "bpSim"]):
                if pdbID in bpSimPerPDB:
                    bpSimPerPDB[pdbID].append( strucPPIs.loc[i, "bpSim"] )
                else:
                    bpSimPerPDB[pdbID] = [ strucPPIs.loc[i, "bpSim"] ]
            
            if not np.isnan(strucPPIs.loc[i, "ccSim"]):
                if pdbID in ccSimPerPDB:
                    ccSimPerPDB[pdbID].append( strucPPIs.loc[i, "ccSim"] )
                else:
                    ccSimPerPDB[pdbID] = [ strucPPIs.loc[i, "ccSim"] ]
            
            if not np.isnan(strucPPIs.loc[i, "coexpr"]):
                if pdbID in coexprPerPDB:
                    coexprPerPDB[pdbID].append( strucPPIs.loc[i, "coexpr"] )
                else:
                    coexprPerPDB[pdbID] = [ strucPPIs.loc[i, "coexpr"] ]
    
    pdbIDs = pd.Series(list(ppisPerPDB.keys()))
    num_ppisPerPDB = pdbIDs.apply(lambda x: len(ppisPerPDB[x]))
    avg_mfSimPerPDB = pdbIDs.apply(lambda x: np.mean(mfSimPerPDB[x]) if x in mfSimPerPDB else np.nan)
    avg_bpSimPerPDB = pdbIDs.apply(lambda x: np.mean(bpSimPerPDB[x]) if x in bpSimPerPDB else np.nan)
    avg_ccSimPerPDB = pdbIDs.apply(lambda x: np.mean(ccSimPerPDB[x]) if x in ccSimPerPDB else np.nan)
    avg_coexprSimPerPDB = pdbIDs.apply(lambda x: np.mean(coexprPerPDB[x]) if x in coexprPerPDB else np.nan)
    
    outPath = interactomeOutDir / 'numPPIsPerPDB.pkl'
    with open(outPath, 'wb') as fout:
        pickle.dump(num_ppisPerPDB, fout)
    
    print( '\n' + 'Interface statistics for structural interactome:' )
    
    maxVal = 10
    num_ppisPerPDB_toplot = num_ppisPerPDB.apply(lambda x: min(x, maxVal))
    num_ppisPerPDB_toplot = [sum(num_ppisPerPDB_toplot == x) for x in np.arange(1, maxVal + 1)]
    
    print('%d out of %d PDB models used for interface annotation are unique to 1 PPI' 
          % (num_ppisPerPDB_toplot[0], sum(num_ppisPerPDB_toplot)))
    
    bar_plot(num_ppisPerPDB_toplot,
             [],
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
    interfaces = pd.DataFrame(columns = ["Protein", "Interface"])
    c = -1
    for _, row in structuralInteractome.iterrows():
        interface_1, interface_2 = row.Interfaces
        c += 1
        interfaces.loc[c] = row.Protein_1, interface_1
        c += 1
        interfaces.loc[c] = row.Protein_2, interface_2
    
    # calculate number of interface residues per interacting protein
    outPath = interactomeOutDir / 'interfaces.pkl'
    with open(outPath, 'wb') as fout:
        pickle.dump(interfaces, fout)
    
    interfaceLength = interfaces["Interface"].apply(len).tolist()
    print('Mean interface length = %.1f residues' % np.mean(interfaceLength))
    multi_histogram_plot ([interfaceLength],
                          ['b'],
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
    with open(ProteinSeqFile, 'rb') as f:
        proteinSeq = pickle.load(f)
    interfaceRatio = ( interfaces["Interface"].apply(len) / 
                       interfaces["Protein"].apply(lambda x: len(proteinSeq[x])) ).tolist()
    print('Mean interface fraction of protein sequence = %.1f%%' % (100 * np.mean(interfaceRatio)))
    multi_histogram_plot ([interfaceRatio],
                          ['g'],
                          xlabel = 'Fraction of interfacial residues',
                          ylabel = 'Number of proteins',
                          edgecolor = 'k',
                          fontsize = 24,
                          bins = 20,
                          alpha = 1,
                          show = showFigs,
                          figdir = figDir,
                          figname = 'interface_fraction')
    
    #------------------------------------------------------------------------------------
    # Plot heatmaps for number of PPIs modeled per PDB structure in relation to
    # functional similarity and tissue co-expression of interaction partners
    #------------------------------------------------------------------------------------
    
    # get all proteins in the structural interactome
    proteinList = structuralInteractome[["Protein_1", "Protein_2"]].values.flatten()
    proteins = []
    for p in proteinList:
        if p not in proteins:
            proteins.append(p)
    
    # plot heatmaps for number of PPIs per protein model along with bar plots of functional
    # similarity of all interaction partners per bin
    for aspect in ['goSim', 'mfSim', 'bpSim', 'ccSim', 'coexpr']:
        allSim = [[]] * len(proteins)
        avgSim = [-2] * len(proteins)
        for i, p in enumerate(proteins):
            sim = strucPPIs.loc[(strucPPIs[["Protein_1","Protein_2"]] == p).any(1), aspect].values
            allSim[i] = [x for x in sim if not np.isnan(x)]
            if len(allSim[i]) > 0:
                avgSim[i] = np.mean(allSim[i])
        
        sortingIndex = sorted(range(len(avgSim)),
                              reverse = True,
                              key = lambda k: avgSim[k])
        sortedProteins = [proteins[i] for i in sortingIndex if avgSim[i] > -2]
        sortedAvgSim = [avgSim[i] for i in sortingIndex if avgSim[i] > -2]
        sortedAllSim = [allSim[i] for i in sortingIndex if avgSim[i] > -2]
        
        mat = ppi_per_pdb_matrix (sortedProteins,
                                  structuralInteractome,
                                  ppisPerPDB,
                                  maximum = False)
       
        # merge number of PPIs per structural model into bins to make heatmap more clear
        cellsPerBin = 10
        mergedAvgSim, ppisPerCell, mergedMatrix = merge_matrix (mat,
                                                                sortedAllSim,
                                                                perBin = cellsPerBin,
                                                                median = True)
        
        # plot functional similarity of all interaction partners per bin
        ylabels = [0, 0.5, 1]
        if aspect == 'goSim':
            ylabel = 'Functional similarity\nwith interaction partner'
        elif aspect == 'mfSim':
            ylabel = 'Molecular function similarity\nwith interaction partner'
        elif aspect == 'bpSim':
            ylabel = 'Biological process similarity\nwith interaction partner'
        elif aspect == 'ccSim':
            ylabel = 'Cellular component similarity\nwith interaction partner'
        elif aspect == 'coexpr':
            ylabel = 'Tissue co-expression\nwith interaction partner'
            coexprMin = np.amin(mergedAvgSim)
            if coexprMin < -0.5:
                ylabels = [-1, -0.5] + ylabels
            elif coexprMin < 0:
                ylabels = [-0.5] + ylabels
        ylim = [ylabels[0], ylabels[-1]]
        bar_plot(mergedAvgSim,
                 [],
                 ylabels = ylabels,
                 ylabel = ylabel,
                 colors = 'steelblue',
                 barwidth = 0.9,
                 ylim = ylim,
                 fontsize = 24,
                 show = showFigs,
                 figdir = figDir,
                 figname = aspect + '_PPIavgSim_heatmap')
        
        # plot heatmap for number of PPIs per protein model, binned by functional similarity
        maxVal = 50
        mergedMatrix[ mergedMatrix > maxVal ] = maxVal
        heatmap_plot(mergedMatrix,
                     cmap = 'summer_r',
                     interpolation = 'nearest',
                     vmin = 1,
                     vmax = maxVal,
                     barPos = [0.85, 0.1, 0.05, 0.6],
                     barTicks = [1, maxVal / 2, maxVal],
                     barLabels = ['1',
                                  str( maxVal / 2 if maxVal % 2 > 0 else int( maxVal / 2 ) ),
                                  '≥ ' + str(maxVal)],
                     show = showFigs,
                     figdir = figDir,
                     figname = aspect + '_numPPIs_heatmap')
        
        # plot heatmap for number of PPIs per cell
        maxVal = 5
        ppisPerCell[ ppisPerCell > maxVal ] = maxVal
        heatmap_plot(ppisPerCell,
                     cmap = 'Blues',
                     interpolation = 'nearest',
                     vmin = 0,
                     vmax = maxVal,
                     barPos = [0.85, 0.1, 0.05, 0.6],
                     barTicks = [0, maxVal],
                     barLabels = ['0',
                                  #str( maxVal / 2 if maxVal % 2 > 0 else int( maxVal / 2 ) ),
                                  '≥ ' + str(maxVal)],
                     show = showFigs,
                     figdir = figDir,
                     figname = aspect + '_numPPIsPerCell_heatmap')

def ppi_per_pdb_matrix (proteinList,
                        structuralInteractome,
                        ppisPerPDB,
                        maximum = True):
    """Create a symmetric matrix of number of PPIs sharing PDB model with each PPI.

    Args:
        proteinList (list): ordered list of proteins.
        structuralInteractome (DataFrame): interactome with chain-pair annotation for each PPI.
        ppisPerPDB (dict): list of PPIs modeled by each PDB ID.
        maximum (boolean): True to select maximum number of PPIs sharing PDB model, 
                           as each PPI may have multiple models (chain-pair annotations), 
                           otherwise minimum is number of PPIs per model is selected.

    """
    mat = pd.DataFrame(np.nan,
                       index = proteinList,
                       columns = proteinList)
    for _, row in structuralInteractome.iterrows():
        if (row.Protein_1 in proteinList) and (row.Protein_2 in proteinList):
            pdbIDs = set([ c1.split('_')[0] for c1, _ in row.Chain_pairs])
            numPPIs = []
            for pdbID in pdbIDs:
                numPPIs.append(len(ppisPerPDB[pdbID]))
        
            if maximum:
                selNum = max(numPPIs)
            else:
                selNum = min(numPPIs)
        
            if proteinList.index(row.Protein_1) < proteinList.index(row.Protein_2):
                mat.loc[row.Protein_1, row.Protein_2] = selNum
            else:
                mat.loc[row.Protein_2, row.Protein_1] = selNum
    mat = mat.astype(float)
    return mat

def merge_matrix (mat,
                  sim,
                  perBin = 10,
                  median = True):
    """Merge numeric entries in a square matrix into cells.

    Args:
        mat (2d array): matrix of numbers to merge.
        sim (List): value (ex, functional similarity) associated with each row/col.
        perBin (int): number of entries per row/col bin.
        median (boolean): True to select the median of all entries per cell,  
                          otherwise the mean of all entries per cell is selected.

    """
    numRows, numCols = mat.shape
    
    i = 0
    avgSim = []
    while i < numRows:
        lastRow = min(i + perBin, numRows)
        allSim = []
        for s in sim[ i : lastRow ]:
            allSim.extend(s)
        avgSim.append( np.mean( allSim ) )
        i += perBin
    
    i = 0
    mergedMatrix = []
    numPPIs = []
    while i < numRows:
        lastRow = min(i + perBin, numRows)
        row = []
        num = []
        j = 0
        while j < numCols:
            subMat = mat.iloc[ i : lastRow, j : min(j + perBin, numCols)]
            ls = [x for x in subMat.values.flatten() if not np.isnan(x)]
            if len(ls) == 0:
                row.append( np.nan )
                num.append(0)
            elif median:
                row.append( np.median(ls) )
                num.append( len(ls) )
            else:
                row.append( np.mean(ls) )
                num.append( len(ls) )
            j += perBin
        mergedMatrix.append( row )
        numPPIs.append( num )
        i += perBin
    return avgSim, np.array( numPPIs ), np.array( mergedMatrix )

if __name__ == "__main__":
    main()
