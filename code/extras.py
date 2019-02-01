#----------------------------------------------------------------------------------------
# This script compares functional similarity and tissue co-expression of partners in
# perturbed interactions and non-perturbed interactions, with all interacting partners and 
# non-interacting partners in the structural interactome. Other network properties are 
# also compared.
#----------------------------------------------------------------------------------------

import os
import pickle
import pandas as pd
import numpy as np
from pathlib import Path
from random import shuffle
from id_mapping import (produce_uniqueGene_swissProtIDs,
                        produce_uniprotID_dict,
                        produce_geneName_dict,
                        produce_protein_interaction_dict)
from interactome_tools import (read_single_interface_annotated_interactome,
                               create_ppi_graph,
                               find_all_shortest_paths,
                               between)
from protein_function import (produce_protein_go_dictionaries,
                              produce_protein_expr_dict,
                              go_sim,
                              partner_sim,
                              coexpr)
from stat_tools import t_test, bootstrap_test, sderror
from plot_tools import box_plot

def main():
    
    # valid reference interactome names
    interactome_names = ['HI-II-14', 'IntAct']
    
    # choose interactome (index in interactome_names)
    interactome_choise = 0
    
    # If True, compute interaction profile similarity, centrality and betweenness for
    # protein pairs with perturbed interaction
    compute_network_properties = False
    
    # show figures
    showFigs = False
    
    # minimum number of tissue expression values required for protein pair tissue
    # co-expression to be considered
    coexprMinTissues = 5
    
    # sample size for non-interacting protein pairs
    nonPPI_sample = 10000
    
    # select reference interactome
    interactome_name = interactome_names[ interactome_choise ]
    
    # interactome colors
    interactome_colors = ['royalblue', 'limegreen']
    
    # figure color for specific interactome 
    interactome_color = interactome_colors[ interactome_choise ]
    
    # directory for data files from external sources
    inDir = Path('../data/external')
    
    # directory to save processed data shared by all interactomes
    outDir = Path('../data/processed')
    
    # directory to save processed data specific to interactome
    interactomeOutDir = outDir / interactome_name
    
    # figure directory
    figDir = Path('../figures') / interactome_name
    
    # create directories if not existing
    if not outDir.exists():
        os.makedirs(outDir)
    if not interactomeOutDir.exists():
        os.makedirs(interactomeOutDir)
    if not figDir.exists():
        os.makedirs(figDir)
    
    #------------------------------------------------------------------------------------
    # produce processed data files if any missing
    #------------------------------------------------------------------------------------
    
    GeneMapFile = outDir / 'to_human_geneName_map.pkl'
    if not GeneMapFile.is_file():
        print('producing UniProtID-to-geneName dictionary')
        produce_geneName_dict(inDir / 'HUMAN_9606_idmapping.dat',
                              inDir / 'uniprot_reviewed_human_proteome.list',
                              GeneMapFile)
    
    UniqueGeneSwissProtIDFile = outDir / 'uniprot_unique_gene_reviewed_human_proteome.list'
    if not UniqueGeneSwissProtIDFile.is_file():
        print('producing list of unique-gene UniProt IDs')
        produce_uniqueGene_swissProtIDs(inDir / 'uniprot_reviewed_human_proteome.list',
                                        GeneMapFile,
                                        UniqueGeneSwissProtIDFile)
    
    UniprotIDmapFile = outDir / 'to_human_uniprotID_map.pkl'
    if not UniprotIDmapFile.is_file():
        print('producing to-UniProt-ID dictionary')
        produce_uniprotID_dict(inDir / 'HUMAN_9606_idmapping.dat',
                               UniqueGeneSwissProtIDFile,
                               UniprotIDmapFile)
    
    #------------------------------------------------------------------------------------
    # load reference and structural interactomes
    #------------------------------------------------------------------------------------
    
    InteractomeFile = interactomeOutDir / 'human_interactome.txt'
    interactome = pd.read_table(InteractomeFile)
    print('\nreference interactome:')
    print('%d PPIs' % len(interactome))
    print('%d proteins' % len(set(interactome[["Protein_1", "Protein_2"]].values.flatten())))
    
    interfaceAnnotatedInteractomeFile = interactomeOutDir / 'human_interface_annotated_interactome.txt'
    annotatedInteractome = read_single_interface_annotated_interactome(interfaceAnnotatedInteractomeFile)
    print('\ninterface-annotated interactome:')
    print('%d PPIs' % len(annotatedInteractome))
    print('%d proteins' % len(set(annotatedInteractome[["Protein_1", "Protein_2"]].values.flatten())))
    
    ProteinPartnersFile = interactomeOutDir / 'protein_interaction_partners.pkl'
    if not ProteinPartnersFile.is_file():
        print('producing protein interaction partners dictionary')
        produce_protein_interaction_dict(InteractomeFile,
                                         ProteinPartnersFile)
    
    if compute_network_properties:
        # compute shortest path between every two proteins
        g = create_ppi_graph(ProteinPartnersFile)
        shortestPathFile = interactomeOutDir / 'shortest_paths.pkl'
        if not shortestPathFile.is_file():
            print('\nfinding shortest paths')
            find_all_shortest_paths (g, shortestPathFile)
        
        # load shortest paths
        with open(shortestPathFile, 'rb') as f:
            paths = pickle.load(f)
        
        # convert from dictionary to Series
        allpaths = []
        for p in list(paths.keys()):
            for p2 in paths.get(p, {}):
                allpaths.append(paths[p][p2])
        allpaths = pd.Series(allpaths)
        allpaths = allpaths[allpaths.apply(len) > 2].reset_index(drop=True)
    
    # produce protein GO association profiles
    GOfile = outDir / 'go.pkl'
    MFfile = outDir / 'gof.pkl'
    BPfile = outDir / 'gop.pkl'
    CCfile = outDir / 'goc.pkl'
    if not GOfile.is_file():
        print('producing protein GO dictionaries')
        produce_protein_go_dictionaries(inDir / 'goa_human.gaf',
                                        GOfile,
                                        MFfile,
                                        BPfile,
                                        CCfile)
    
    # produce protein tissue expression profiles
    proteinExprFile = outDir / 'proteinExpr.pkl'
    if not proteinExprFile.is_file():
        print('producing protein tissue expression dictionary')
        produce_protein_expr_dict(inDir / 'E-MTAB-513.tsv.txt',
                                  UniprotIDmapFile,
                                  proteinExprFile,
                                  headers = list(range(1, 18)))
    
    with open(GOfile, 'rb') as f:
        GOassoc = pickle.load(f)
    with open(MFfile, 'rb') as f:
        MFassoc = pickle.load(f)
    with open(BPfile, 'rb') as f:
        BPassoc = pickle.load(f)
    with open(CCfile, 'rb') as f:
        CCassoc = pickle.load(f)
    with open(proteinExprFile, 'rb') as f:
        expr = pickle.load(f)
    with open(ProteinPartnersFile, 'rb') as f:
        ppiPartners = pickle.load(f)
    
    # print number of GO terms for each category
    terms = set()
    for k in GOassoc.keys():
        terms.update(GOassoc[k])
    print('\n%d gene ontology (GO) terms' % len(terms))
    terms = set()
    for k in MFassoc.keys():
        terms.update(MFassoc[k])
    print('%d molecular function (MF) terms' % len(terms))
    terms = set()
    for k in BPassoc.keys():
        terms.update(BPassoc[k])
    print('%d biological process (BP) terms' % len(terms))
    terms = set()
    for k in CCassoc.keys():
        terms.update(CCassoc[k])
    print('%d cellular component (CC) terms' % len(terms))
    
    #------------------------------------------------------------------------------------
    # Calculate GO similarity, co-expression and network properties for interacting 
    # proteins perturbed by common mutation
    #------------------------------------------------------------------------------------
    
    # load mutation perturbations
    filteredMutationPerturbFile = interactomeOutDir / 'filtered_mutation_perturbs.pkl'
    if filteredMutationPerturbFile.is_file():
        with open(filteredMutationPerturbFile, 'rb') as f:
            naturalPerturbs, diseasePerturbs = pickle.load(f)
    else:
        print('\nMutation perturbations file not found')
        print('Quiting ...')
        return
    
    # select proteins and their partners whose interactions are perturbed by common mutation
    naturalPerturbed = pd.DataFrame()
    naturalPerturbed["protein"] = naturalPerturbs["Protein"].values
    naturalPerturbed["partners"] = naturalPerturbs.apply(lambda x:
                                                         [p for i, p in enumerate(x["partners"]) 
                                                          if x["perturbations"][i] > 0],
                                                         axis=1)
    naturalPerturbed = naturalPerturbed[naturalPerturbed["partners"].apply(len) > 0].reset_index(drop=True)
    
    # count PPIs perturbed by common mutation, no duplicates
    naturalPerturbedPPIs = []
    for i, row in naturalPerturbed.iterrows():
        filteredPartners = []
        for p in row.partners:
            ppi = {row.protein, p}
            if ppi not in naturalPerturbedPPIs:
                naturalPerturbedPPIs.append(ppi)
                filteredPartners.append(p)
        naturalPerturbed.loc[i, "partners"] = filteredPartners
    naturalPerturbed = naturalPerturbed[naturalPerturbed["partners"].apply(len) > 0].reset_index(drop=True)
    print('\n%d perturbed PPIs in response to non-disease mutations' % sum(naturalPerturbed["partners"].apply(len)))
    
    # Calculate functional similarity and network centrality for protein pairs whose 
    # interaction is perturbed by natural non-disease mutation
    naturalPerturbed["goSim"] = np.nan
    naturalPerturbed["mfSim"] = np.nan
    naturalPerturbed["bpSim"] = np.nan
    naturalPerturbed["ccSim"] = np.nan
    naturalPerturbed["coexpr"] = np.nan
    naturalPerturbed["partnerSim"] = np.nan
    naturalPerturbed["centrality"] = np.nan
    naturalPerturbed["partner_centrality"] = np.nan
    naturalPerturbed["betweenness"] = np.nan
    
    naturalPerturbed["goSim"] = naturalPerturbed.apply(lambda x:
                                                       [go_sim(x["protein"],
                                                               p,
                                                               GOassoc) 
                                                        for p in x["partners"]],
                                                       axis=1)
    naturalPerturbed["mfSim"] = naturalPerturbed.apply(lambda x:
                                                       [go_sim(x["protein"],
                                                               p,
                                                               MFassoc) 
                                                        for p in x["partners"]],
                                                       axis=1)
    naturalPerturbed["bpSim"] = naturalPerturbed.apply(lambda x:
                                                       [go_sim(x["protein"],
                                                               p,
                                                               BPassoc) 
                                                        for p in x["partners"]],
                                                       axis=1)
    naturalPerturbed["ccSim"] = naturalPerturbed.apply(lambda x:
                                                       [go_sim(x["protein"],
                                                               p,
                                                               CCassoc) 
                                                        for p in x["partners"]],
                                                       axis=1)
    naturalPerturbed["coexpr"] = naturalPerturbed.apply(lambda x:
                                                        [coexpr(x["protein"],
                                                                p,
                                                                expr,
                                                                minTissues = coexprMinTissues) 
                                                         for p in x["partners"]],
                                                        axis=1)
    
    if compute_network_properties:
        naturalPerturbed["partnerSim"] = naturalPerturbed.apply(lambda x:
                                                                [partner_sim(x["protein"],
                                                                             p,
                                                                             ppiPartners) 
                                                                 for p in x["partners"]],
                                                                axis=1)
        maxDegree = max(map(len, ppiPartners.values()))
        naturalPerturbed["centrality"] = naturalPerturbed.apply(lambda x:
                                                                [len(ppiPartners[x["protein"]]) / maxDegree]
                                                                if x["protein"] in ppiPartners
                                                                else [],
                                                                axis=1)
        naturalPerturbed["partner_centrality"] = naturalPerturbed.apply(lambda x:
                                                                        [len(ppiPartners[p]) / maxDegree
                                                                         for p in x["partners"] 
                                                                         if p in ppiPartners],
                                                                        axis=1)
        naturalPerturbed["betweenness"] = naturalPerturbed.apply(lambda x:
                                                                 [between(allpaths, x["protein"])],
                                                                 axis=1)
    
    #------------------------------------------------------------------------------------
    # Calculate GO similarity, co-expression and network properties for interacting 
    # proteins perturbed by disease mutation
    #------------------------------------------------------------------------------------
    
    # select proteins and their partners whose interactions are perturbed by disease mutation
    diseasePerturbed = pd.DataFrame()
    diseasePerturbed["protein"] = diseasePerturbs["Protein"].values
    diseasePerturbed["partners"] = diseasePerturbs.apply(lambda x:
                                                         [p for i, p in enumerate(x["partners"]) 
                                                          if x["perturbations"][i] > 0],
                                                         axis=1)
    diseasePerturbed = diseasePerturbed[diseasePerturbed["partners"].apply(len) > 0].reset_index(drop=True)
    
    # count PPIs perturbed by disease mutation, no duplicates
    diseasePerturbedPPIs = []
    for i, row in diseasePerturbed.iterrows():
        filteredPartners = []
        for p in row.partners:
            ppi = {row.protein, p}
            if ppi not in diseasePerturbedPPIs:
                diseasePerturbedPPIs.append(ppi)
                filteredPartners.append(p)
        diseasePerturbed.loc[i, "partners"] = filteredPartners
    diseasePerturbed = diseasePerturbed[diseasePerturbed["partners"].apply(len) > 0].reset_index(drop=True)
    print('%d perturbed PPIs in response to disease mutations' % sum(diseasePerturbed["partners"].apply(len)))
    
    diseasePerturbed["goSim"] = np.nan
    diseasePerturbed["mfSim"] = np.nan
    diseasePerturbed["bpSim"] = np.nan
    diseasePerturbed["ccSim"] = np.nan
    diseasePerturbed["coexpr"] = np.nan
    diseasePerturbed["partnerSim"] = np.nan
    diseasePerturbed["centrality"] = np.nan
    diseasePerturbed["partner_centrality"] = np.nan
    diseasePerturbed["betweenness"] = np.nan

    diseasePerturbed["goSim"] = diseasePerturbed.apply(lambda x:
                                                       [go_sim(x["protein"],
                                                               p,
                                                               GOassoc) 
                                                        for p in x["partners"]],
                                                       axis=1)
    diseasePerturbed["mfSim"] = diseasePerturbed.apply(lambda x:
                                                       [go_sim(x["protein"],
                                                               p,
                                                               MFassoc) 
                                                        for p in x["partners"]],
                                                       axis=1)
    diseasePerturbed["bpSim"] = diseasePerturbed.apply(lambda x:
                                                       [go_sim(x["protein"],
                                                               p,
                                                               BPassoc) 
                                                        for p in x["partners"]],
                                                       axis=1)
    diseasePerturbed["ccSim"] = diseasePerturbed.apply(lambda x:
                                                       [go_sim(x["protein"],
                                                               p,
                                                               CCassoc) 
                                                        for p in x["partners"]],
                                                       axis=1)
    diseasePerturbed["coexpr"] = diseasePerturbed.apply(lambda x:
                                                        [coexpr(x["protein"],
                                                                p,
                                                                expr,
                                                                minTissues = coexprMinTissues) 
                                                         for p in x["partners"]],
                                                        axis=1)
    if compute_network_properties:
        diseasePerturbed["partnerSim"] = diseasePerturbed.apply(lambda x:
                                                                [partner_sim(x["protein"],
                                                                             p,
                                                                             ppiPartners) 
                                                                 for p in x["partners"]],
                                                                axis=1)
        diseasePerturbed["centrality"] = diseasePerturbed.apply(lambda x:
                                                                [len(ppiPartners[x["protein"]]) / maxDegree]
                                                                if x["protein"] in ppiPartners
                                                                else [],
                                                                axis=1)
        diseasePerturbed["partner_centrality"] = diseasePerturbed.apply(lambda x:
                                                                        [len(ppiPartners[p]) / maxDegree
                                                                         for p in x["partners"] 
                                                                         if p in ppiPartners],
                                                                        axis=1)
        diseasePerturbed["betweenness"] = diseasePerturbed.apply(lambda x:
                                                                 [between(allpaths, x["protein"])],
                                                                 axis=1)
    
    #------------------------------------------------------------------------------------
    # Calculate GO similarity and co-expression for all interacting protein pairs in the
    # structural interactome
    #------------------------------------------------------------------------------------
    
    strucPPIs = pd.DataFrame()
    strucPPIs["Protein_1"] = annotatedInteractome["Protein_1"].copy()
    strucPPIs["Protein_2"] = annotatedInteractome["Protein_2"].copy()    
    strucPPIs["goSim"] = strucPPIs.apply(lambda x:
                                     go_sim(x["Protein_1"],
                                            x["Protein_2"],
                                            GOassoc),
                                     axis=1)
    strucPPIs["mfSim"] = strucPPIs.apply(lambda x:
                                     go_sim(x["Protein_1"],
                                            x["Protein_2"],
                                            MFassoc),
                                     axis=1)
    strucPPIs["bpSim"] = strucPPIs.apply(lambda x:
                                     go_sim(x["Protein_1"],
                                            x["Protein_2"],
                                            BPassoc),
                                     axis=1)
    strucPPIs["ccSim"] = strucPPIs.apply(lambda x:
                                     go_sim(x["Protein_1"],
                                            x["Protein_2"],
                                            CCassoc),
                                     axis=1)
    strucPPIs["coexpr"] = strucPPIs.apply(lambda x:
                                      coexpr(x["Protein_1"],
                                             x["Protein_2"],
                                             expr,
                                             minTissues = coexprMinTissues),
                                      axis=1)
    strucPPIs = strucPPIs.applymap(lambda x: [x] if isinstance(x, float) else x)
    
    #------------------------------------------------------------------------------------
    # Calculate GO similarity and co-expression for all interacting non-perturbed protein 
    # pairs in the structural interactome
    #------------------------------------------------------------------------------------
    
    nonPerturbedPPIs = strucPPIs.copy()
    keep = nonPerturbedPPIs.apply(lambda x:
                                  ({x["Protein_1"], x["Protein_2"]} not in naturalPerturbedPPIs)
                                  and 
                                  ({x["Protein_1"], x["Protein_2"]} not in diseasePerturbedPPIs),
                                  axis=1)
    nonPerturbedPPIs = nonPerturbedPPIs[keep].reset_index(drop=True)
    
    #------------------------------------------------------------------------------------
    # Calculate GO similarity and co-expression for a sample of non-interacting protein 
    # pairs in the reference interactome, proteins still need to be in the structural
    # interactome
    #------------------------------------------------------------------------------------
    
    # obtain a sample of non-interacting protein pairs
    proteins = list(set(annotatedInteractome[["Protein_1", "Protein_2"]].values.flatten()))
    pairs = []
    while len(pairs) < nonPPI_sample:
        shuffle(proteins)
        p1 = proteins[0]
        shuffle(proteins)
        p2 = proteins[0]
        pair = {p1, p2}
        if (p1 != p2) and (p1 not in ppiPartners[p2]):
            pairs.append(pair)
    
    p1list = []
    p2list = []
    for p1, p2 in pairs:
        p1list.append(p1)
        p2list.append(p2)
    nonPPIs = pd.DataFrame()
    nonPPIs["Protein_1"] = p1list
    nonPPIs["Protein_2"] = p2list   
    
    nonPPIs["goSim"] = nonPPIs.apply(lambda x:
                                     go_sim(x["Protein_1"],
                                            x["Protein_2"],
                                            GOassoc),
                                     axis=1)
    nonPPIs["mfSim"] = nonPPIs.apply(lambda x:
                                     go_sim(x["Protein_1"],
                                            x["Protein_2"],
                                            MFassoc),
                                     axis=1)
    nonPPIs["bpSim"] = nonPPIs.apply(lambda x:
                                     go_sim(x["Protein_1"],
                                            x["Protein_2"],
                                            BPassoc),
                                     axis=1)
    nonPPIs["ccSim"] = nonPPIs.apply(lambda x:
                                     go_sim(x["Protein_1"],
                                            x["Protein_2"],
                                            CCassoc),
                                     axis=1)
    nonPPIs["coexpr"] = nonPPIs.apply(lambda x:
                                      coexpr(x["Protein_1"],
                                             x["Protein_2"],
                                             expr,
                                             minTissues = coexprMinTissues),
                                      axis=1)
    nonPPIs = nonPPIs.applymap(lambda x: [x] if isinstance(x, float) else x)
    
    #------------------------------------------------------------------------------------
    # Print out GO similarity, co-expression and network properties for all pairs
    #------------------------------------------------------------------------------------
    
    bootstrapItr = 1000
    
    (naturalPerturbed_goSim_mean,
     diseasePerturbed_goSim_mean,
     nonPerturbedPPIs_goSim_mean,
     strucPPIs_goSim_mean,
     nonPPIs_goSim_mean,
     naturalPerturbed_goSim_error,
     diseasePerturbed_goSim_error,
     nonPerturbedPPIs_goSim_error,
     strucPPIs_goSim_error,
     nonPPIs_goSim_error) = compare([naturalPerturbed["goSim"].values,
                                     diseasePerturbed["goSim"].values,
                                     nonPerturbedPPIs["goSim"].values,
                                     strucPPIs["goSim"].values,
                                     nonPPIs["goSim"].values],
                                    ('PPI perturbed\nby non-disease\nmutation',
                                     'PPI perturbed\nby disease\nmutation',
                                     'PPI not perturbed\nby mutation',
                                     'All interacting\npairs',
                                     'Non-interacting\npairs'),
                                    ylabel = 'Similarity in GO association\nof pairs of proteins',
                                    ylim = [0, 1.1],
                                    xbounds = [1, 5],
                                    ybounds = [0, 1],
                                    colors = ['turquoise', 'orangered', 'green', 'steelblue', 'gold'],
                                    bootstrapItr = bootstrapItr,
                                    showfig = showFigs,
                                    figdir = figDir,
                                    figname = 'Perturbation_goSim')
    
    (naturalPerturbed_mfSim_mean,
     diseasePerturbed_mfSim_mean,
     nonPerturbedPPIs_mfSim_mean,
     strucPPIs_mfSim_mean,
     nonPPIs_mfSim_mean,
     naturalPerturbed_mfSim_error,
     diseasePerturbed_mfSim_error,
     nonPerturbedPPIs_mfSim_error,
     strucPPIs_mfSim_error,
     nonPPIs_mfSim_error) = compare([naturalPerturbed["mfSim"].values,
                                     diseasePerturbed["mfSim"].values,
                                     nonPerturbedPPIs["mfSim"].values,
                                     strucPPIs["mfSim"].values,
                                     nonPPIs["mfSim"].values],
                                    ('PPI perturbed\nby non-disease\nmutation',
                                     'PPI perturbed\nby disease\nmutation',
                                     'PPI not perturbed\nby mutation',
                                     'All interacting\npairs',
                                     'Non-interacting\npairs'),
                                    ylabel = 'Similarity in molecular function\nof pairs of proteins',
                                    ylim = [0, 1.1],
                                    xbounds = [1, 5],
                                    ybounds = [0, 1],
                                    colors = ['turquoise', 'orangered', 'green', 'steelblue', 'gold'],
                                    bootstrapItr = bootstrapItr,
                                    showfig = showFigs,
                                    figdir = figDir,
                                    figname = 'Perturbation_mfSim')
    
    (naturalPerturbed_bpSim_mean,
     diseasePerturbed_bpSim_mean,
     nonPerturbedPPIs_bpSim_mean,
     strucPPIs_bpSim_mean,
     nonPPIs_bpSim_mean,
     naturalPerturbed_bpSim_error,
     diseasePerturbed_bpSim_error,
     nonPerturbedPPIs_bpSim_error,
     strucPPIs_bpSim_error,
     nonPPIs_bpSim_error) = compare([naturalPerturbed["bpSim"].values,
                                     diseasePerturbed["bpSim"].values,
                                     nonPerturbedPPIs["bpSim"].values,
                                     strucPPIs["bpSim"].values,
                                     nonPPIs["bpSim"].values],
                                    ('PPI perturbed\nby non-disease\nmutation',
                                     'PPI perturbed\nby disease\nmutation',
                                     'PPI not perturbed\nby mutation',
                                     'All interacting\npairs',
                                     'Non-interacting\npairs'),
                                    ylabel = 'Similarity in biological process\nof pairs of proteins',
                                    ylim = [0, 1.1],
                                    xbounds = [1, 5],
                                    ybounds = [0, 1],
                                    colors = ['turquoise', 'orangered', 'green', 'steelblue', 'gold'],
                                    bootstrapItr = bootstrapItr,
                                    showfig = showFigs,
                                    figdir = figDir,
                                    figname = 'Perturbation_bpSim')
    
    (naturalPerturbed_ccSim_mean,
     diseasePerturbed_ccSim_mean,
     nonPerturbedPPIs_ccSim_mean,
     strucPPIs_ccSim_mean,
     nonPPIs_ccSim_mean,
     naturalPerturbed_ccSim_error,
     diseasePerturbed_ccSim_error,
     nonPerturbedPPIs_ccSim_error,
     strucPPIs_ccSim_error,
     nonPPIs_ccSim_error) = compare([naturalPerturbed["ccSim"].values,
                                     diseasePerturbed["ccSim"].values,
                                     nonPerturbedPPIs["ccSim"].values,
                                     strucPPIs["ccSim"].values,
                                     nonPPIs["ccSim"].values],
                                    ('PPI perturbed\nby non-disease\nmutation',
                                     'PPI perturbed\nby disease\nmutation',
                                     'PPI not perturbed\nby mutation',
                                     'All interacting\npairs',
                                     'Non-interacting\npairs'),
                                    ylabel = 'Similarity in cellular component\nof pairs of proteins',
                                    ylim = [0, 1.1],
                                    xbounds = [1, 5],
                                    ybounds = [0, 1],
                                    colors = ['turquoise', 'orangered', 'green', 'steelblue', 'gold'],
                                    bootstrapItr = bootstrapItr,
                                    showfig = showFigs,
                                    figdir = figDir,
                                    figname = 'Perturbation_ccSim')
    
    (naturalPerturbed_coexpr_mean,
     diseasePerturbed_coexpr_mean,
     nonPerturbedPPIs_coexpr_mean,
     strucPPIs_coexpr_mean,
     nonPPIs_coexpr_mean,
     naturalPerturbed_coexpr_error,
     diseasePerturbed_coexpr_error,
     nonPerturbedPPIs_coexpr_error,
     strucPPIs_coexpr_error,
     nonPPIs_coexpr_error) = compare([naturalPerturbed["coexpr"].values,
                                      diseasePerturbed["coexpr"].values,
                                      nonPerturbedPPIs["coexpr"].values,
                                      strucPPIs["coexpr"].values,
                                      nonPPIs["coexpr"].values],
                                     ('PPI perturbed\nby non-disease\nmutation',
                                      'PPI perturbed\nby disease\nmutation',
                                      'PPI not perturbed\nby mutation',
                                      'All interacting\npairs',
                                      'Non-interacting\npairs'),
                                     ylabel = 'Coexpression of pairs of proteins',
                                     ylim = [-1, 1],
                                     xbounds = [1, 5],
                                     colors = ['turquoise', 'orangered', 'green', 'steelblue', 'gold'],
                                     bootstrapItr = bootstrapItr,
                                     showfig = showFigs,
                                     figdir = figDir,
                                     figname = 'Perturbation_coexpr')
    
    if compute_network_properties:
        (naturalPerturbed_partnerSim_mean,
         diseasePerturbed_partnerSim_mean,
         naturalPerturbed_partnerSim_error,
         diseasePerturbed_partnerSim_error) = compare([naturalPerturbed["partnerSim"].values,
                                                       diseasePerturbed["partnerSim"].values],
                                                      ('PPI perturbed by\nnon-disease mutation',
                                                       'PPI perturbed by\ndisease mutation'),
                                                      ylabel = 'Fraction of partners shared\nby pairs of proteins',
                                                      colors = ['turquoise', 'orangered'],
                                                      bootstrapItr = bootstrapItr)
        
        (naturalPerturbed_centrality_mean,
         diseasePerturbed_centrality_mean,
         naturalPerturbed_centrality_error,
         diseasePerturbed_centrality_error) = compare([naturalPerturbed["centrality"].values,
                                                       diseasePerturbed["centrality"].values],
                                                      ('PPI perturbed by\nnon-disease mutation',
                                                       'PPI perturbed by\ndisease mutation'),
                                                      ylabel = 'Centrality of mutated protein (# of partners)',
                                                      colors = ['turquoise', 'orangered'],
                                                      bootstrapItr = bootstrapItr)
        
        (naturalPerturbed_partner_centrality_mean,
         diseasePerturbed_partner_centrality_mean,
         naturalPerturbed_partner_centrality_error,
         diseasePerturbed_partner_centrality_error) = compare([naturalPerturbed["partner_centrality"].values,
                                                               diseasePerturbed["partner_centrality"].values],
                                                              ('PPI perturbed by\nnon-disease mutation',
                                                               'PPI perturbed by\ndisease mutation'),
                                                              ylabel = 'Centrality of mutated protein''s partner (# of partners)',
                                                              colors = ['turquoise', 'orangered'],
                                                              bootstrapItr = bootstrapItr)
        
        (naturalPerturbed_betweenness_mean,
         diseasePerturbed_betweenness_mean,
         naturalPerturbed_betweenness_error,
         diseasePerturbed_betweenness_error) = compare([naturalPerturbed["betweenness"].apply(lambda x: [100*e for e in x]).values,
                                                        diseasePerturbed["betweenness"].apply(lambda x: [100*e for e in x]).values],
                                                       ('PPI perturbed by\nnon-disease mutation',
                                                        'PPI perturbed by\ndisease mutation'),
                                                       ylabel = 'Betweenness (% of all paths)',
                                                       colors = ['turquoise', 'orangered'],
                                                       bootstrapItr = bootstrapItr)

def compare (samples,
             groups,
             xlabel = None,
             ylabel = None,
             xlim = None,
             ylim = None,
             xbounds = None,
             ybounds = None,
             colors = None,
             bootstrapItr = None,
             showfig = True,
             figdir = None,
             figname = None):
    """Compare multiple samples, each one containing list elements.

    Args:
        samples (list): list of samples to compare.
        groups (list): list of sample names.
        xlabel (str): label for plot x-axis.
        ylabel (str): label for plot y-axis.
        xlim (list): range for plot x-axis.
        ylim (list): range for plot y-axis.
        xbounds (list): boundary for displayed x-axis.
        ybounds (list): boundary for displayed y-axis.
        colors (list): color for each box.
        bootstrapItr (int): number of resamplings for bootstrap test. If not specified, 
                            a two-sided t-test is used instead.
        showfig (boolean): True to show plot, otherwise plot is not shown.
        figdir (str): directory to save figure in.
        figname (str): name of file to save figure in. 
    """
    if ylabel is not None:
        print('\n' + ylabel.replace('\n',' '))
    numGroups = len(groups)
    processedSamples = []
    means = [np.nan] * numGroups
    errors = [np.nan] * numGroups
    labels = [''] * numGroups
    numValid = 0
    for i, sample in enumerate(samples):
        s = []
        for ls in sample:
            s.extend(ls)
        s = [e for e in s if not np.isnan(e)]
        processedSamples.append(s)    
        num = len(s)
        if num > 0:
            numValid += 1
            means[i] = np.mean(s)
            errors[i] = sderror(s)
            labels[i] = groups[i]
            print(groups[i].replace('\n',' ') + ': %.3f (SE = %g, n = %d)' % (means[i],
                                                                              errors[i],
                                                                              num))
    
    if (numValid > 0) and (figdir is not None) and (figname is not None):
        toplot = [i for i, s in enumerate(processedSamples) if len(s) > 2]
        box_plot(data = [processedSamples[i] for i in toplot],
                 xlabels = [labels[i] for i in toplot],
                 xlabel = xlabel,
                 ylabel = ylabel,
                 xlim = xlim,
                 ylim = ylim,
                 fontsize = 10,
                 adjustBottom = 0.25,
                 shiftBottomAxis = -0.1,
                 xbounds = xbounds,
                 ybounds = ybounds,
                 colors = colors,
                 show = showfig,
                 figdir = figdir,
                 figname = figname)
    if numValid > 1:
        for i, s1 in enumerate(processedSamples[:-1]):
            if len(s1) > 0:
                for j, s2 in enumerate(processedSamples[i + 1 :]):
                    if len(s2) > 0:
                        print('Statistical significance: ' + 
                              labels[i].replace('\n','') + 
                              ', ' + 
                              labels[i + j + 1].replace('\n',''))
                        if bootstrapItr is None:
                            t_test(s1, s2)
                        else:
                            bootstrap_test(s1, s2, bootstrapItr)
    means.extend(errors)
    return means

if __name__ == "__main__":
    main()
