#----------------------------------------------------------------------------------------
# This script creates plots for interface mapping description and fraction of junk PPIs
# for all structural interactomes together.
# Run script 'predict_edgotypes_and_junkPPIs.py' before running this script.
#----------------------------------------------------------------------------------------

import os
import pickle
import numpy as np
from pathlib import Path
from stat_tools import sderror
from plot_tools import (bar_plot,
                        multi_bar_plot,
                        multi_histogram_plot,
                        venn2_plot,
                        venn3_plot)

def main():
    
    # valid binary interactome names
    interactome_names = ['HI-II-14', 'IntAct']
    
    struc_interactome_names = ['Y2H-SI', 'IntAct-SI']
    
    interactome_colors = ['royalblue', 'limegreen']
    
    # show figures
    showFigs = False
    
    # directory for data files from external sources
    inDir = Path('../data/processed')
    
    # figure directory specific to interactome
    figDir = Path('../figures/combined')
    
    if not figDir.exists():
        os.makedirs(figDir)
    
    junk_dict = dict()
    bounds_dict = dict()
    for name in interactome_names:
        interactomeInDir = inDir / name
        junkPPIFile = interactomeInDir / 'fraction_junk_PPIs.pkl'
        with open(junkPPIFile, 'rb') as f:
            junk_dict[name], bounds_dict[name] = pickle.load(f)
    
    numInteractomes = len( interactome_names )
    confidenceIntervals = len( bounds_dict[ interactome_names[ 0 ] ] ) > 0
    
    # plot the fraction of junk PPIs calculated without knowledge of binding energy change
    pN_E_results = ( [ junk_dict[name][0] for name in interactome_names ]
                     + [ junk_dict[interactome_names[0]][2] ]
                     + [ junk_dict[name][3] for name in interactome_names ]
                     + [ junk_dict[interactome_names[0]][5] ] )
    
    if confidenceIntervals:
        pN_E_bounds = ( [ bounds_dict[name][0] for name in interactome_names ]
                        + [ bounds_dict[interactome_names[0]][2] ]
                        + [ bounds_dict[name][3] for name in interactome_names ]
                        + [ bounds_dict[interactome_names[0]][5] ] )
    
    if confidenceIntervals:
        upper = [ p + upper[0] for p, (lower, upper) in zip(pN_E_results, pN_E_bounds) ]
        maxY = max( upper[ : numInteractomes + 1 ] )
    else:
        maxY = max( pN_E_results[ : numInteractomes + 1 ] )
    maxY = 5 * np.ceil( maxY / 5 )
    bar_plot(pN_E_results[ : numInteractomes + 1 ],
             pN_E_bounds[  : numInteractomes + 1 ] if confidenceIntervals else [ 0 ],
             #xlabels = [ 'Prediction\n(%s)' % name for name in interactome_names ] + [ 'Experiment' ],
             xlabels = [ name for name in struc_interactome_names ] + [ 'Experiment' ],
             ylabels = np.arange(0, maxY + 5, 5),
             ylabel = ('Fraction of Junk PPIs (%)'),
             colors = 'white',
             fmt = 'k.',
             capsize = 10 if confidenceIntervals else 0,
             msize = 11,
             ewidth = 1,
             #ecolors = interactome_colors + [ 'orangered' ],
             ecolors = 'black',
             fontsize = 18,
             xlim = [0.8, 3.1],
             ylim = [0, maxY],
             yMinorTicks = True,
             adjustBottom = 0.2,
             shiftBottomAxis = -0.1,
             xbounds = (1,3),
             show = showFigs,
             figdir = figDir,
             figname = 'Fraction_junk_PPIs_geometrybased_perMut')
    if confidenceIntervals:
        upper = [ p + upper[0] for p, (lower, upper) in zip(pN_E_results, pN_E_bounds) ]
        maxY = max( upper[ numInteractomes + 1 : ] )
    else:
        maxY = max( pN_E_results[ numInteractomes + 1 : ] )
    maxY = 5 * np.ceil( maxY / 5 )
    bar_plot(pN_E_results[ numInteractomes + 1 : ],
             pN_E_bounds[  numInteractomes + 1 : ] if confidenceIntervals else [ 0 ],
             #xlabels = [ 'Prediction\n(%s)' % name for name in interactome_names ] + [ 'Experiment' ],
             xlabels = [ name for name in struc_interactome_names ] + [ 'Experiment' ],
             ylabels = np.arange(0, maxY + 5, 5),
             ylabel = ('Fraction of Junk PPIs (%)'),
             colors = 'white',
             fmt = 'k.',
             capsize = 10 if confidenceIntervals else 0,
             msize = 11,
             ewidth = 1,
             #ecolors = interactome_colors + [ 'orangered' ],
             ecolors = 'black',
             fontsize = 18,
             xlim = [0.8, 3.1],
             ylim = [0, maxY],
             yMinorTicks = True,
             adjustBottom = 0.2,
             shiftBottomAxis = -0.1,
             xbounds = (1,3),
             show = showFigs,
             figdir = figDir,
             figname = 'Fraction_junk_PPIs_geometrybased_perPPI')
    
    # plot the fraction of junk PPIs calculated using knowledge of binding energy change
    pN_E_results = ( [ junk_dict[name][1] for name in interactome_names ]
                     + [ junk_dict[interactome_names[0]][2] ]
                     + [ junk_dict[name][4] for name in interactome_names ]
                     + [ junk_dict[interactome_names[0]][5] ] )
    
    if confidenceIntervals:
        pN_E_bounds = ( [ bounds_dict[name][1] for name in interactome_names ]
                        + [ bounds_dict[interactome_names[0]][2] ]
                        + [ bounds_dict[name][4] for name in interactome_names ]
                        + [ bounds_dict[interactome_names[0]][5] ] )
    
    if confidenceIntervals:
        upper = [ p + upper[0] for p, (lower, upper) in zip(pN_E_results, pN_E_bounds) ]
        maxY = max( upper[ : numInteractomes + 1 ] )
    else:
        maxY = max( pN_E_results[ : numInteractomes + 1 ] )
    maxY = 5 * np.ceil( maxY / 5 )
    bar_plot(pN_E_results[ : numInteractomes + 1 ],
             pN_E_bounds[  : numInteractomes + 1 ] if confidenceIntervals else [ 0 ],
             #xlabels = [ 'Prediction\n(%s)' % name for name in interactome_names ] + [ 'Experiment' ],
             xlabels = [ name for name in struc_interactome_names ] + [ 'Experiment' ],
             ylabels = np.arange(0, maxY + 5, 5),
             ylabel = ('Fraction of junk PPIs (%)'),
             colors = 'white',
             fmt = 'k.',
             capsize = 10 if confidenceIntervals else 0,
             msize = 11,
             ewidth = 1,
             #ecolors = interactome_colors + [ 'orangered' ],
             ecolors = 'black',
             fontsize = 18,
             xlim = [0.8, 3.1],
             ylim = [0, maxY],
             yMinorTicks = True,
             adjustBottom = 0.2,
             shiftBottomAxis = -0.1,
             xbounds = (1,3),
             show = showFigs,
             figdir = figDir,
             figname = 'Fraction_junk_PPIs_physicsbased_perMut')
    if confidenceIntervals:
        upper = [ p + upper[0] for p, (lower, upper) in zip(pN_E_results, pN_E_bounds) ]
        maxY = max( upper[ numInteractomes + 1 : ] )
    else:
        maxY = max( pN_E_results[ numInteractomes + 1 : ] )
    maxY = 5 * np.ceil( maxY / 5 )
    bar_plot(pN_E_results[ numInteractomes + 1 : ],
             pN_E_bounds[  numInteractomes + 1 : ] if confidenceIntervals else [ 0 ],
             #xlabels = [ 'Prediction\n(%s)' % name for name in interactome_names ] + [ 'Experiment' ],
             xlabels = [ name for name in struc_interactome_names ] + [ 'Experiment' ],
             ylabels = np.arange(0, maxY + 5, 5),
             ylabel = ('Fraction of junk PPIs (%)'),
             colors = 'white',
             fmt = 'k.',
             capsize = 10 if confidenceIntervals else 0,
             msize = 11,
             ewidth = 1,
             #ecolors = interactome_colors + [ 'orangered' ],
             ecolors = 'black',
             fontsize = 18,
             xlim = [0.8, 3.1],
             ylim = [0, maxY],
             yMinorTicks = True,
             adjustBottom = 0.2,
             shiftBottomAxis = -0.1,
             xbounds = (1,3),
             show = showFigs,
             figdir = figDir,
             figname = 'Fraction_junk_PPIs_physicsbased_perPPI')
    
    #------------------------------------------------------------------------------------
    # functional similarity of interaction partners in reference and structural 
    # interactomes
    #------------------------------------------------------------------------------------
    
    # gene ontology similarity
    gosim = []
    error = []
    for name in interactome_names:
        interactomeInDir = inDir / name
        inPath = interactomeInDir / 'gosim_reference_structural.pkl'
        with open(inPath, 'rb') as fin:
            binary_gosim, structural_gosim =  pickle.load(fin)
        valid_gosim = [ [sim for sim in binary_gosim if not np.isnan(sim)],
                        [sim for sim in structural_gosim if not np.isnan(sim)] ]
        gosim.append( list( map(np.mean, valid_gosim) ) )
        error.append( list( map(sderror, valid_gosim) ) )
    
    multi_bar_plot(gosim,
                   error,
                   xlabels = ('Reference\ninteractome', 'Structural\ninteractome'),
                   ylabel = 'Functional similarity of\ninteracting protein pairs',
                   ylabels = [0, 0.1, 0.2, 0.3],
                   colors = interactome_colors,
                   edgecolor = 'k',
                   barwidth = 0.3,
                   fontsize = 24,
                   leg = ['HI-II-14', 'IntAct'],
                   show = showFigs,
                   figdir = figDir,
                   figname = 'interactome_gosim')
    
    # tissue coexpression
    coexpr = []
    error = []
    for name in interactome_names:
        interactomeInDir = inDir / name
        inPath = interactomeInDir / 'coexpr_reference_structural.pkl'
        with open(inPath, 'rb') as fin:
            binary_coexpr, structural_coexpr =  pickle.load(fin)
        valid_coexpr = [ [c for c in binary_coexpr if not np.isnan(c)],
                         [c for c in structural_coexpr if not np.isnan(c)] ]
        coexpr.append( list( map(np.mean, valid_coexpr) ) )
        error.append( list( map(sderror, valid_coexpr) ) )
    
    multi_bar_plot(coexpr,
                   error,
                   xlabels = ('Reference\ninteractome', 'Structural\ninteractome'),
                   ylabel = 'Tissue co-expression of\ninteracting protein pairs',
                   ylabels = [0, 0.1, 0.2, 0.3, 0.4, 0.5],
                   ylim = [0, 0.5],
                   colors = interactome_colors,
                   edgecolor = 'k',
                   barwidth = 0.3,
                   fontsize = 24,
                   leg = ['HI-II-14', 'IntAct'],
                   show = showFigs,
                   figdir = figDir,
                   figname = 'interactome_coexpr')
    
    #------------------------------------------------------------------------------------
    # interface description
    #------------------------------------------------------------------------------------
    
    # plot number of PPIs modeled per PDB structure
    num_ppisPerPDB = []
    for name in interactome_names:
        interactomeInDir = inDir / name
        inPath = interactomeInDir / 'numPPIsPerPDB.pkl'
        with open(inPath, 'rb') as fin:
            num_ppisPerPDB.append( pickle.load(fin) )
    
    maxVal = 10
    num_ppisPerPDB_toplot = [ numppis.apply(lambda numppis: min(numppis, maxVal)) 
                              for numppis in num_ppisPerPDB ]
    num_ppisPerPDB_toplot = [ [sum(numppis == x) for x in np.arange(1, maxVal + 1)] 
                              for numppis in num_ppisPerPDB_toplot ]
    multi_bar_plot(num_ppisPerPDB_toplot,
                   [],
                   xlabels = list(map(str, np.arange(1, maxVal))) + [ str(maxVal) + '+'],
                   xlabel = 'Number of PPIs modeled by PDB structure',
                   ylabel = 'Number of PDB structures',
                   ylabels = [0, 500, 1000, 1500, 2000],
                   ylim = [0, 2000],
                   colors = interactome_colors,
                   barwidth = 0.4,
                   fontsize = 22,
                   leg = ['HI-II-14', 'IntAct'],
                   show = showFigs,
                   figdir = figDir,
                   figname = 'numPPIsPerPDB')
    
    # compile interfaces for all interacting proteins
    interfaces = []
    for name in interactome_names:
        interactomeInDir = inDir / name
        inPath = interactomeInDir / 'interfaces.pkl'
        with open(inPath, 'rb') as fin:
            interfaces.append( pickle.load(fin) )
    
    # plot distribution of interface length
    interfaceLength = [ intf["Interface"].apply(len).tolist() for intf in interfaces ]
    multi_histogram_plot (list( reversed( interfaceLength ) ),
                          list( reversed( interactome_colors ) ),
                          xlabel = 'Number of interfacial residues',
                          ylabel = 'Number of proteins',
                          ylabels = [0, 250, 500, 750, 1000, 1250],
                          ylim = [0, 1250],
                          edgecolor = 'k',
                          fontsize = 22,
                          bins = 25,
                          alpha = 1,
                          show = showFigs,
                          figdir = figDir,
                          figname = 'interface_length')
    
    # plot distribution of fraction of interfacial residues
    ProteinSeqFile = inDir / 'human_reference_sequences.pkl'
    with open(ProteinSeqFile, 'rb') as f:
        proteinSeq = pickle.load(f)
    interfaceRatio = [ ( intf["Interface"].apply(len) / 
                         intf["Protein"].apply(lambda x: len(proteinSeq[x])) ).tolist()
                       for intf in interfaces ]
    multi_histogram_plot (list( reversed( interfaceRatio ) ),
                          list( reversed( interactome_colors ) ),
                          xlabel = 'Fraction of interfacial residues',
                          ylabel = 'Number of proteins',
                          edgecolor = 'k',
                          fontsize = 24,
                          bins = 25,
                          alpha = 1,
                          show = showFigs,
                          figdir = figDir,
                          figname = 'interface_frac')
    
    #------------------------------------------------------------------------------------
    # protein overlap with experiment
    #------------------------------------------------------------------------------------
    
    for name, col in zip(interactome_names, interactome_colors):
        proteins = []
        interactomeInDir = inDir / name
        proteinFile = interactomeInDir / 'nondisease_mut_proteins.list'
        with open(proteinFile, 'rb') as f:
            proteins.append( set( f.read().split() ) )
        proteinFile = interactomeInDir / 'nondisease_mut_proteins_exp.list'
        with open(proteinFile, 'rb') as f:
            proteins.append( set( f.read().split() ) )
        venn2_plot(proteins,
               #labels = ['Prediction', 'Experiment'],
               colors = [col, 'orangered'],
               hideNum = True,
               show = showFigs,
               figdir = figDir,
               figname = name + '_non_disease_protein_overlap')
        
        proteins = []
        interactomeInDir = inDir / name
        proteinFile = interactomeInDir / 'disease_mut_proteins.list'
        with open(proteinFile, 'rb') as f:
            proteins.append( set( f.read().split() ) )
        proteinFile = interactomeInDir / 'disease_mut_proteins_exp.list'
        with open(proteinFile, 'rb') as f:
            proteins.append( set( f.read().split() ) )
        venn2_plot(proteins,
               #labels = ['Prediction', 'Experiment'],
               colors = [col, 'orangered'],
               hideNum = True,
               show = showFigs,
               figdir = figDir,
               figname = name + '_disease_protein_overlap')
    
    proteins = []
    for name in interactome_names:
        interactomeInDir = inDir / name
        proteinFile = interactomeInDir / 'nondisease_mut_proteins.list'
        with open(proteinFile, 'rb') as f:
            proteins.append( set( f.read().split() ) )
    proteinFile = interactomeInDir / 'nondisease_mut_proteins_exp.list'
    with open(proteinFile, 'rb') as f:
        proteins.append( set( f.read().split() ) )
    
    venn3_plot(proteins,
               #labels = ['Prediction', 'Experiment'],
               colors = interactome_colors + [ 'orangered' ],
               hideNum = True,
               show = showFigs,
               figdir = figDir,
               figname = 'non_disease_protein_overlap')
    
    proteins = []
    for name in interactome_names:
        interactomeInDir = inDir / name
        proteinFile = interactomeInDir / 'disease_mut_proteins.list'
        with open(proteinFile, 'rb') as f:
            proteins.append( set( f.read().split() ) )
    proteinFile = interactomeInDir / 'disease_mut_proteins_exp.list'
    with open(proteinFile, 'rb') as f:
        proteins.append( set( f.read().split() ) )
    
    venn3_plot(proteins,
               #labels = ['Prediction', 'Experiment'],
               colors = interactome_colors + [ 'orangered' ],
               hideNum = True,
               show = showFigs,
               figdir = figDir,
               figname = 'disease_protein_overlap')
    
if __name__ == "__main__":
    main()
