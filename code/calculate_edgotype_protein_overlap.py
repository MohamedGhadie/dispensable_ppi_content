#----------------------------------------------------------------------------------------
# Calculate overlap in proteins covered by edgotyped mutations in two interactomes.
#
# Run the following scripts before running this script:
# - produce_data_mappings.py
# - calculate_junk_content_geometry.py (for prediction edgotype data)
# - calculate_junk_content_experiment.py (for experiment edgotype data)
#----------------------------------------------------------------------------------------

import os
import pickle
import pandas as pd
from pathlib import Path
from plot_tools import venn2_plot

def main():
    
    # names of reference interactomes
    # options: any pair from: HI-II-14, IntAct, experiment
    interactome_names = ['IntAct', 'experiment']
    
    # structural interactome names
    struc_name = {'HI-II-14':'Y2H-SI', 'IntAct':'IntAct-SI', 'experiment':'Experiment'}
    
    # interactome colors
    interactome_color = {'HI-II-14':'royalblue', 'IntAct':'limegreen', 'experiment':'orangered'}
    
    # name of gene name column to be used if protein column not found
    gene_col = 'Symbol'
    
    # show figures
    showFigs = False
    
    # parent directory of all data files
    dataDir = Path('../data')
    
    # parent directory of all processed data files
    procDir = dataDir / 'processed'
    
    # figure directory
    figDir = Path('../figures') / 'combined'
    
    label1 = 'experiment' if interactome_names[0] is 'experiment' else 'geometry'
    label2 = 'experiment' if interactome_names[1] is 'experiment' else 'geometry'
    
    # input data files
    uniprotIDmapFile = procDir / 'to_human_uniprotID_map.pkl'
    natMutEdgotypeFile1 = procDir / interactome_names[0] / ('nondisease_mutation_edgotype_%s.txt' % label1)
    disMutEdgotypeFile1 = procDir / interactome_names[0] / ('disease_mutation_edgotype_%s.txt' % label1)
    natMutEdgotypeFile2 = procDir / interactome_names[1] / ('nondisease_mutation_edgotype_%s.txt' % label2)
    disMutEdgotypeFile2 = procDir / interactome_names[1] / ('disease_mutation_edgotype_%s.txt' % label2)
    
    # create output directories if not existing
    if not figDir.exists():
        os.makedirs(figDir)
    
    with open(uniprotIDmapFile, 'rb') as f:
        uniprotID = pickle.load(f)
    
    #----------------------------------------------------------------------------------------
    # process proteins covered by first interactome
    #----------------------------------------------------------------------------------------

    naturalMutations1 = pd.read_table (natMutEdgotypeFile1, sep='\t')
    diseaseMutations1 = pd.read_table (disMutEdgotypeFile1, sep='\t')
    
    if 'protein' not in naturalMutations1.columns:
        naturalMutations1["protein"] = naturalMutations1[gene_col].apply(lambda x:
                                                                         uniprotID[x] if x in uniprotID 
                                                                         else x + 'unknown')
    if 'protein' not in diseaseMutations1.columns:
        diseaseMutations1["protein"] = diseaseMutations1[gene_col].apply(lambda x:
                                                                         uniprotID[x] if x in uniprotID 
                                                                         else x + 'unknown')
    natMutProteins1 = set(naturalMutations1["protein"].values)
    disMutProteins1 = set(diseaseMutations1["protein"].values)
    allMutProteins1 = natMutProteins1 | disMutProteins1
    
    print( '\n' + '%s mutations:' % struc_name[interactome_names[0]])
    print('Non-disease: %d mutations on %d proteins' % (len(naturalMutations1), len(natMutProteins1)))
    print('Disease: %d mutations on %d proteins' % (len(diseaseMutations1), len(disMutProteins1)))
    print('Total mutated proteins: %d' % len(allMutProteins1))
    
    #----------------------------------------------------------------------------------------
    # process proteins covered by second interactome
    #----------------------------------------------------------------------------------------
    
    naturalMutations2 = pd.read_table (natMutEdgotypeFile2, sep='\t')
    diseaseMutations2 = pd.read_table (disMutEdgotypeFile2, sep='\t')
    
    if 'protein' not in naturalMutations2.columns:
        naturalMutations2["protein"] = naturalMutations2[gene_col].apply(lambda x:
                                                                         uniprotID[x] if x in uniprotID 
                                                                         else x + 'unknown')
    if 'protein' not in diseaseMutations2.columns:
        diseaseMutations2["protein"] = diseaseMutations2[gene_col].apply(lambda x:
                                                                         uniprotID[x] if x in uniprotID 
                                                                         else x + 'unknown')
    natMutProteins2 = set(naturalMutations2["protein"].values)
    disMutProteins2 = set(diseaseMutations2["protein"].values)
    allMutProteins2 = natMutProteins2 | disMutProteins2
    
    print( '\n' + '%s mutations:' % struc_name[interactome_names[1]])
    print('Non-disease: %d mutations on %d proteins' % (len(naturalMutations2), len(natMutProteins2)))
    print('Disease: %d mutations on %d proteins' % (len(diseaseMutations2), len(disMutProteins2)))
    print('Total mutated proteins: %d' % len(allMutProteins2))
    
    #----------------------------------------------------------------------------------------
    # calculate overlap in protein space
    #----------------------------------------------------------------------------------------
    
    natMutProtein_overlap = natMutProteins1 & natMutProteins2
    disMutProtein_overlap = disMutProteins1 & disMutProteins2
    allMutProtein_overlap = allMutProteins1 & allMutProteins2
    
    print( '\n' + 'Overlap in proteins covered by all mutations:' )
    print( '%s proteins: %.1f %% (%d out of %d)' % (struc_name[interactome_names[0]],
                                                    100 * len(allMutProtein_overlap) / len(allMutProteins1),
                                                    len(allMutProtein_overlap),
                                                    len(allMutProteins1)) )
    print( '%s proteins: %.1f %% (%d out of %d)' % (struc_name[interactome_names[1]],
                                                    100 * len(allMutProtein_overlap) / len(allMutProteins2),
                                                    len(allMutProtein_overlap),
                                                    len(allMutProteins2)) )
    
    print( '\n' + 'Overlap in proteins covered by non-disease mutations:' )
    print( '%s proteins: %.1f %% (%d out of %d)' % (struc_name[interactome_names[0]],
                                                    100 * len(natMutProtein_overlap) / len(natMutProteins1),
                                                    len(natMutProtein_overlap),
                                                    len(natMutProteins1)) )
    print( '%s proteins: %.1f %% (%d out of %d)' % (struc_name[interactome_names[1]],
                                                    100 * len(natMutProtein_overlap)  / len(natMutProteins2),
                                                    len(natMutProtein_overlap),
                                                    len(natMutProteins2)) )
    
    print( '\n' + 'Overlap in proteins covered by disease mutations:' )
    print( '%s proteins: %.1f %% (%d out of %d)' % (struc_name[interactome_names[0]],
                                                    100 * len(disMutProtein_overlap) / len(disMutProteins1),
                                                    len(disMutProtein_overlap),
                                                    len(disMutProteins1) ) )
    print( '%s proteins: %.1f %% (%d out of %d)' % (struc_name[interactome_names[1]],
                                                    100 * len(disMutProtein_overlap) / len(disMutProteins2),
                                                    len(disMutProtein_overlap),
                                                    len(disMutProteins2) ) )
    
    # plot overlap in proteins covered by edgotyped nondisease mutations
    label = tuple([struc_name[i] for i in interactome_names])
    venn2_plot([natMutProteins1, natMutProteins2],
               labels = [struc_name[i] for i in interactome_names],
               colors = [interactome_color[i] for i in interactome_names],
               hideNum = True,
               show = showFigs,
               figdir = figDir,
               figname = 'non_disease_mut_protein_overlap_%s_%s' % label)
    # plot overlap in proteins covered by edgotyped disease mutations
    venn2_plot([disMutProteins1, disMutProteins2],
               labels = [struc_name[i] for i in interactome_names],
               colors = [interactome_color[i] for i in interactome_names],
               hideNum = True,
               show = showFigs,
               figdir = figDir,
               figname = 'disease_mut_protein_overlap_%s_%s' % label)
    # plot overlap in proteins covered by disease and non-disease edgotyped mutations
    venn2_plot([allMutProteins1, allMutProteins2],
               labels = [struc_name[i] for i in interactome_names],
               colors = [interactome_color[i] for i in interactome_names],
               hideNum = True,
               show = showFigs,
               figdir = figDir,
               figname = 'all_mut_protein_overlap_%s_%s' % label)

if __name__ == '__main__':
    main()
