#----------------------------------------------------------------------------------------
# This script constructs a structural interactome from a reference interactome by mapping 
# PPI binding interfaces, at atomic resolution, from experimentally determined 
# three-dimensional structural models in PDB onto PPIs in the reference interactome. 
# Next, Mendelian disease-causing mutations and common neutral mutations not associated 
# with disease are mapped onto the structural interactome, and PPI perturbations by 
# mutations are predicted. Mutation edgotypes ("edgetic" and "non-edgetic") per mutation 
# and per PPI are predicted using predicted PPI perturbations. Then, the fraction of junk 
# PPIs (PPIs neutral upon perturbation) in the structural interactome is estimated using 
# predicted mutation edgotypes, and compared to the fraction of junk PPIs estimated using 
# mutation edgotypes determined from experiments.
#
# Run script 'data_initial_processing.py' for preprocessing select data files before 
# running this script.
#----------------------------------------------------------------------------------------

import os
import io
import pickle
import pandas as pd
import numpy as np
from pathlib import Path
from text_tools import (parse_HI_II_14_interactome,
                        parse_IntAct_interactions,
                        parse_fasta_file,
                        parse_blast_file,
                        merge_refSeq_sequence_files,
                        parse_HGMD_disease_mutations,
                        filter_and_merge_dbsnp_mutations,
                        translate_HGMD_flanking_sequences,
                        remove_hgmd_synon_nonsense_mutations,
                        truncate_HGMD_flanking_sequences,
                        get_flanking_sequences,
                        match_flanking_sequences,
                        match_masked_flanking_sequences,
                        remove_dbsnp_synon_nonsense_mutations)
from id_mapping import (produce_uniqueGene_swissProtIDs,
                        produce_uniqueGene_sequences,
                        produce_proteinSeq_dict,
                        produce_uniprotID_dict,
                        produce_geneName_dict,
                        produce_protein_interaction_dict,
                        produce_protein_chain_dict,
                        produce_protein_chain_alignment_dict,
                        produce_chain_dict,
                        produce_chainSeq_dict)
from interactome_tools import (write_interactome_sequences,
                               duplicated_PPIs,
                               remove_duplicate_PPIs,
                               write_chain_annotated_interactome_to_excel,
                               read_interface_annotated_interactome,
                               write_interface_annotated_interactome_to_excel)
from structural_annotation import (produce_alignment_evalue_dict,
                                   locate_alignments,
                                   filter_chain_annotations,
                                   produce_chain_annotated_interactome,
                                   produce_interface_annotated_interactome,
                                   merge_interactome_interface_annotations,
                                   remove_duplicate_interface_annotations,
                                   read_chain_annotated_interactome,
                                   write_pdb_mapped_mutations,
                                   read_mutation_ddg)
from pdb_tools import (download_pdb_files,
                       produce_chain_list,
                       produce_pdb_ids,
                       produce_seqres_order_dict)
from protein_function import (produce_protein_go_dictionaries,
                              produce_protein_expr_dict,
                              go_sim,
                              coexpr)
from stat_tools import (t_test,
                        fisher_test,
                        bootstrap_test,
                        sderror,
                        sderror_on_fraction,
                        proportion_ratio_CI)
from plot_tools import (bar_plot,
                        multi_bar_plot,
                        box_plot,
                        multi_histogram_plot,
                        pie_plot,
                        venn2_plot,
                        network_plot)
from mutation_interface_edgotype import (mutation_PPI_interface_perturbations,
                                         energy_based_perturb_prediction)

def main():
    
    # valid reference interactome names
    interactome_names = ['HI-II-14', 'IntAct']
    
    # choose interactome (index in interactome_names)
    interactome_choise = 1
    
    # calculate confidence interval for the fraction of junk PPIs
    computeConfidenceIntervals = True
    
    # treat experiment quasi-null mutation as edgetic
    consider_exp_QuasiNull_perturbs = False
    
    # % confidence interval
    CI = 95
    
    # show figures
    showFigs = False
    
    # download missing PDB structures whose chain pairs map onto interactome
    download_missing_structures = False
    
    # select reference interactome
    interactome_name = interactome_names[ interactome_choise ]
    
    # interactome colors
    interactome_colors = ['royalblue', 'limegreen']
    
    # color specific to interactome
    interactome_color = interactome_colors[ interactome_choise ]
    
    # directory for data files from external sources
    inDir = Path('../data/external')
     
    # directory to save processed data shared by all interactomes
    outDir = Path('../data/processed')
    
    # directory to save processed data specific to interactome
    interactomeOutDir = outDir / interactome_name
    
    # figure directory
    if consider_exp_QuasiNull_perturbs:
        figDir = Path('../figures') / interactome_name / 'exp_QuasiNull_perturbs_considered'
    else:
        figDir = Path('../figures') / interactome_name
        
    # directory for PDB structure files
    pdbDir = inDir / 'pdb_files'
    
    # create output directories if not existing
    if not outDir.exists():
        os.makedirs( outDir )
    if not interactomeOutDir.exists():
        os.makedirs( interactomeOutDir )
    if not figDir.exists():
        os.makedirs( figDir )
    
    #------------------------------------------------------------------------------------
    # produce processed data files if any missing
    #------------------------------------------------------------------------------------
    
    SequenceFile = outDir / 'human_reference_sequences.txt'
    if not SequenceFile.is_file():
        print( 'Reading protein sequence fasta file' )
        parse_fasta_file(inDir / 'human_reference_sequences.fasta', 
                         SequenceFile)
    
    GeneMapFile = outDir / 'to_human_geneName_map.pkl'
    if not GeneMapFile.is_file():
        print( 'Producing UniProtID-to-geneName dictionary' )
        produce_geneName_dict(inDir / 'HUMAN_9606_idmapping.dat',
                              inDir / 'uniprot_reviewed_human_proteome.list',
                              GeneMapFile)
    
    UniqueGeneSwissProtIDFile = outDir / 'uniprot_unique_gene_reviewed_human_proteome.list'
    if not UniqueGeneSwissProtIDFile.is_file():
        print( 'Producing list of unique-gene UniProt IDs' )
        produce_uniqueGene_swissProtIDs(inDir / 'uniprot_reviewed_human_proteome.list',
                                        GeneMapFile,
                                        UniqueGeneSwissProtIDFile)
    
    UniqueGeneSequenceFile = outDir / 'human_unique_gene_reference_sequences.txt'
    if not UniqueGeneSequenceFile.is_file():
        print( 'Producing sequence file for unique-gene UniProt IDs' )
        produce_uniqueGene_sequences(SequenceFile,
                                     UniqueGeneSwissProtIDFile,
                                     GeneMapFile,
                                     UniqueGeneSequenceFile)
    
    ProteinSeqFile = outDir / 'human_reference_sequences.pkl'
    if not ProteinSeqFile.is_file():
        print( 'Producing protein sequence dictionary' )
        produce_proteinSeq_dict(inDir / 'human_reference_sequences.fasta',
                                ProteinSeqFile)
    
    UniprotIDmapFile = outDir / 'to_human_uniprotID_map.pkl'
    if not UniprotIDmapFile.is_file():
        print( 'Producing to-UniProt-ID dictionary' )
        produce_uniprotID_dict(inDir / 'HUMAN_9606_idmapping.dat',
                               UniqueGeneSwissProtIDFile,
                               UniprotIDmapFile)
    
    InteractomeFile1 = interactomeOutDir / 'human_interactome_all.txt'
    if not InteractomeFile1.is_file():
        print( 'Parsing interactome dataset' )
        if interactome_name == 'HI-II-14':
            parse_HI_II_14_interactome(inDir / 'HI-II-14.tsv',
                                       UniprotIDmapFile,
                                       InteractomeFile1,
                                       selfPPIs = False)
        elif interactome_name == 'IntAct':
            parse_IntAct_interactions(inDir / 'intact.txt', 
                                      UniqueGeneSwissProtIDFile,
                                      InteractomeFile1,
                                      geneMapFile = GeneMapFile,
                                      selfPPIs = False)
        else:
            print( '\n' + 'Error: interactome name %s not recognized. Quiting...' % interactome_name )
            return
    
    InteractomeFile = interactomeOutDir / 'human_interactome.txt'
    if not InteractomeFile.is_file():
        interactome = pd.read_table(InteractomeFile1)
        print( '\n' + 'Initial interactome size: %d PPIs' % len(interactome) )
        if interactome_name == 'IntAct':
            interactome = duplicated_PPIs ( interactome )
            print( '\n' + 'Interactome size after removing PPIs reported only once: %d PPIs' % len(interactome) )
        interactome = remove_duplicate_PPIs( interactome )
        print( '\n' + 'Interactome size after removing duplicate PPIs: %d PPIs' % len(interactome) )
        interactome.to_csv(InteractomeFile, index=False, sep='\t')
    
    InteractomeSequenceFile = interactomeOutDir / 'interactome_sequences.fasta'
    if not InteractomeSequenceFile.is_file():
        print( 'Writing interactome protein sequences' )
        write_interactome_sequences(InteractomeFile,
                                    UniqueGeneSequenceFile,
                                    InteractomeSequenceFile)
    
    ProteinPartnersFile = interactomeOutDir / 'protein_interaction_partners.pkl'
    if not ProteinPartnersFile.is_file():
        print( 'Producing protein interaction partners dictionary')
        produce_protein_interaction_dict(InteractomeFile,
                                         ProteinPartnersFile)
    
    chainSeqFile = outDir / 'chain_sequences.pkl'
    if not chainSeqFile.is_file():
        print( 'Producing PDB chain sequence dictionary from fasta records' )
        produce_chainSeq_dict(inDir / 'pdb_seqres_reduced.fasta',
                              chainSeqFile)
    
    chainListFile = outDir / 'pdb_seqres_chains.list'
    if not chainListFile.is_file():
        print( 'Producing PDB chain ID file from fasta records' )
        produce_chain_list(inDir / 'pdb_seqres_reduced.fasta',
                           chainListFile)
    
    pdbidFile = outDir / 'seqres_pdb_IDs.list'
    if not pdbidFile.is_file():
        print( 'Producing unique PDB ID file from chain list file' )
        produce_pdb_ids(chainListFile,
                        pdbidFile)
    
    pdbChainsFile = outDir / 'pdb_seqres_chains.pkl'
    if not pdbChainsFile.is_file():    
        print( 'Producing PDB chain dictionary from chain list file' )
        produce_chain_dict(chainListFile,
                           pdbChainsFile)
    
    chainResOrderFile = outDir / 'chain_seqres_order.pkl'
    if not chainResOrderFile.is_file():
        print( 'Producing PDB chain structured-residue order dictionary' )
        produce_seqres_order_dict(inDir / 'ss_dis.txt',
                                  pdbDir,
                                  chainSeqFile,
                                  chainResOrderFile)
    
    chainMapFile1 = outDir / 'human_pdb_alignment.txt'
    # pausing time in seconds after processing each 50 million lines in BLAST file
    pausetime = 0
    if not chainMapFile1.is_file():
        print( 'Parsing BLAST protein-chain alignment file' )
        parse_blast_file(inDir / 'human_pdb_e-5',
                         'us-ascii',
                         pausetime,
                         chainMapFile1)
    
    chainMapFile2 = outDir / 'human_pdb_chain_map.txt'
    # If True, chain residue is required to match aligned protein residue for their 
    # positions to be considered aligned
    resMatch = False
    # pausing time in seconds between locating alignments on queries and subjects
    pausetime = 120
    if not chainMapFile2.is_file():
        print( 'Locating aligned residues on protein and chain sequences' )
        locate_alignments(chainMapFile1,
                          chainMapFile2,
                          resMatch = resMatch,
                          pausetime = pausetime)
    
    chainMapFile = outDir / 'human_pdb_chain_map_filtered.txt'
    proteinChainsFile = outDir / 'protein_chains.pkl'
    chainAnnotatedInteractomeFile = interactomeOutDir / 'human_chain_annotated_interactome.txt'
    chainAlignFile = outDir / 'human_protein_chain_alignments.pkl'
    alignmentEvalueFile = outDir / 'human_protein_chain_min_alignment_evalues.pkl'
    # Maximum e-value cutoff to filter out protein chain annotations
    evalue = 1e-10
    # Minimum chain coverage fraction required for annotation to be kept
    chainCoverage = 0
    if not chainMapFile.is_file():
        print( 'Filtering chain annotations' )
        filter_chain_annotations(chainMapFile2,
                                 evalue,
                                 chainCoverage,
                                 chainMapFile)
    
        print( 'Producing protein chains dictionary' )
        produce_protein_chain_dict(chainMapFile,
                                   proteinChainsFile)
        
        print( 'Producing protein-chain alignment evalue dictionary' )
        produce_alignment_evalue_dict (chainMapFile,
                                       alignmentEvalueFile,
                                       method = 'min')
        
        print( 'Producing chain-annotated interactome' )
        produce_chain_annotated_interactome(InteractomeFile,
                                            proteinChainsFile,
                                            chainAnnotatedInteractomeFile,
                                            alignmentEvalueFile = alignmentEvalueFile)
        
#         print('producing protein-chain alignment dictionary')
#         produce_protein_chain_alignment_dict(chainMapFile,
#                                              chainAlignFile)
    
    chainAnnotatedInteractome = read_chain_annotated_interactome (chainAnnotatedInteractomeFile)
    uniqueChains = set()
    for ls in chainAnnotatedInteractome["Mapping_chains"].values:
        for pair in ls:
            uniqueChains.update(pair)
    uniqueChains = pd.Series(list(uniqueChains))
    uniquePDBs = uniqueChains.apply(lambda x: x[ : x.find('_') ]).drop_duplicates()
    print( '\n' + 'Interactome chain-pair annotations:' )
    print( '%d unique chains in %d unique PDB structures' % (len(uniqueChains), len(uniquePDBs)) )
    uniqueChains.to_csv(interactomeOutDir / 'interactome_chainIDs.txt', index=False)
    uniquePDBs.to_csv(interactomeOutDir / 'interactome_pdbIDs.txt', index=False)
    print( 'Unique chain IDs and unique PDB IDs written to file' )
    
    if download_missing_structures:
        print( 'Downloading missing structures for PDB IDs mapping onto interactome' )
        download_pdb_files(interactomeOutDir / 'interactome_pdbIDs.txt',
                           pdbDir)
    
    interfaceAnnotatedInteractomeFile1 = interactomeOutDir / 'human_interface_annotated_interactome_withDuplicates.txt'
    # file directory to save computed PDB chain binding interfaces
    chainInterfaceFile = outDir / 'pdb_interfaces.txt'
    # consider only interfaces with this minimum coverage fraction successfully mapped onto PPI
    mapCutoff = 0.5
    # max binding distance for interface residues in PDB structure
    bindingDist = 5
    # If True, merge all mapped interfaces of a PPI into one interface
    mergeInterfaces = True
    # max number of interfaces mapped from distinct chain-pair annotations for each PPI
    maxInterfaces = 5
    # max number of chain-pair interface calculations per PPI, including known interfaces
    maxAttempts = 100
    # If True, randomly sample from a PPI's chain pair annotations, otherwise start
    # from first annotation
    randChainPairs = False
    if not interfaceAnnotatedInteractomeFile1.is_file():
        print( 'Mapping chain interfaces onto chain-annotated interactome' )
        produce_interface_annotated_interactome(chainAnnotatedInteractomeFile,
                                                pdbDir,
                                                chainMapFile,
                                                chainInterfaceFile,
                                                chainResOrderFile,
                                                maxInterfaces,
                                                maxAttempts,
                                                randChainPairs,
                                                mapCutoff,
                                                bindingDist,
                                                interfaceAnnotatedInteractomeFile1)
    
    interfaceAnnotatedInteractomeFile = interactomeOutDir / 'human_interface_annotated_interactome.txt'
    if not interfaceAnnotatedInteractomeFile.is_file():
        if mergeInterfaces:
            print( 'Merging interface annotations for each PPI' )
            merge_interactome_interface_annotations(interfaceAnnotatedInteractomeFile1,
                                                    interfaceAnnotatedInteractomeFile)
        else:
            print( 'Removing duplicate interface annotations for each PPI without merging' )
            remove_duplicate_interface_annotations(interfaceAnnotatedInteractomeFile1,
                                                   interfaceAnnotatedInteractomeFile)
    
    #------------------------------------------------------------------------------------
    # process disease mutations
    #------------------------------------------------------------------------------------
    
    hgmdMutationsFile1 = outDir / 'HGMD_missense_nonsense_DM_mutations.txt'
    if not hgmdMutationsFile1.is_file():
        print( 'Parsing HGMD mutations file' )
        parse_HGMD_disease_mutations(inDir / 'HGMD_database_2011.txt',
                                     UniprotIDmapFile,
                                     hgmdMutationsFile1)
    
    hgmdMutationsFile2 = outDir / 'HGMD_missense_nonsense_DM_mutations_translated.txt'
    if not hgmdMutationsFile2.is_file():
        print( 'Translating HGMD mutation flanking sequences' )
        translate_HGMD_flanking_sequences(hgmdMutationsFile1,
                                          hgmdMutationsFile2)
    
    hgmdMutationsFile3 = outDir / 'HGMD_missense_DM_mutations.txt'
    if not hgmdMutationsFile3.is_file():
        print( 'Removing HGMD synonymous and nonsense mutations' )
        remove_hgmd_synon_nonsense_mutations(hgmdMutationsFile2,
                                             hgmdMutationsFile3)
    
    hgmdMutationsFile4 = outDir / 'HGMD_missense_DM_mutations_truncated.txt'
    if not hgmdMutationsFile4.is_file():
        print( 'Truncating HGMD mutation flanking sequences' )
        truncate_HGMD_flanking_sequences(hgmdMutationsFile3,
                                         hgmdMutationsFile4)
    
    hgmdMutationsFile = outDir / 'HGMD_missense_DM_mutations_matched.txt'
    if not hgmdMutationsFile.is_file():
        print( 'Mapping HGMD mutation flanking sequences onto UniProt sequences' )
        match_flanking_sequences(hgmdMutationsFile4,
                                 UniqueGeneSequenceFile,
                                 hgmdMutationsFile)
    
    #------------------------------------------------------------------------------------
    # process non-disease mutations
    #------------------------------------------------------------------------------------
    
    # number of RefSeq sequence files to merge
    numRefSeqFiles = 39
    refseqInterDir = outDir / 'refSeq_intermediate'
    refseqFile = outDir / ('refseq_human_protein_sequences.txt')
    if not refseqFile.is_file():
        print( 'Merging protein RefSeq sequences from all files' )
        merge_refSeq_sequence_files(refseqInterDir,
                                    numRefSeqFiles,
                                    refseqFile)
    
    dbsnpOutDir = outDir / 'dbsnp_intermediate'
    if not dbsnpOutDir.exists():
        os.makedirs(dbsnpOutDir)
    dbsnpMutationsFile1 = outDir / 'dbsnp_mutations1.txt'
    dbsnpMutationsFile2 = outDir / 'dbsnp_mutations2.txt'
    dbsnpMutationsFile3 = outDir / 'dbsnp_mutations3.txt'
    dbsnpMutationsFile = outDir / 'dbsnp_mutations4.txt'
    if not dbsnpMutationsFile1.is_file():
        # number of residues included in each side of mutation flanking sequence
        flankingSeqSideLen = 10
        # pausing time in seconds after processing each dbSNP chromosome file
        pausetime = 60
        print( 'Merging dbSNP mutation files' )
        filter_and_merge_dbsnp_mutations(dbsnpOutDir,
                                         UniprotIDmapFile,
                                         pausetime,
                                         dbsnpMutationsFile1)
    
        print( 'Reading dbSNP mutation flanking sequences from RefSeq transcripts' )
        get_flanking_sequences(dbsnpMutationsFile1,
                               refseqFile,
                               flankingSeqSideLen,
                               dbsnpMutationsFile2)
        
        print( 'Removing dbSNP synonymous and nonsense mutations' )
        remove_dbsnp_synon_nonsense_mutations(dbsnpMutationsFile2,
                                              dbsnpMutationsFile3)

        print( 'Mapping dbSNP flanking sequences onto UniProt sequences' )
        match_masked_flanking_sequences(dbsnpMutationsFile3,
                                        UniqueGeneSequenceFile,
                                        dbsnpMutationsFile)
    
    #------------------------------------------------------------------------------------
    # load reference and structural interactomes
    #------------------------------------------------------------------------------------
    
    interactome = pd.read_table(InteractomeFile)
    print( '\n' + 'reference interactome:' )
    print( '%d PPIs' % len(interactome) )
    print( '%d proteins' % len(set(interactome[["Protein_1", "Protein_2"]].values.flatten())) )
    
    writer = pd.ExcelWriter( str( interactomeOutDir / ( interactome_name + '_reference_interactome.xlsx' ) ) )
    interactome.to_excel(writer, sheet_name = interactome_name + '_reference_interactome')
    writer.save()
    
    with open(ProteinSeqFile, 'rb') as f:
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
    
    chainAnnotatedInteractome = read_chain_annotated_interactome(chainAnnotatedInteractomeFile)
    print( '\n' + 'chain-annotated interactome:' )
    print( '%d PPIs' % len(chainAnnotatedInteractome) )
    print( '%d proteins' % len(set(chainAnnotatedInteractome[["Protein_1", "Protein_2"]].values.flatten())) )
    
    write_chain_annotated_interactome_to_excel (
                    chainAnnotatedInteractome,
                    interactomeOutDir / ( interactome_name + '_PPI_chain_annotations.xlsx' ),
                    sheet_name = 'PPI_chain_annotations')
    
    annotatedInteractome = read_interface_annotated_interactome(interfaceAnnotatedInteractomeFile)
    print( '\n' + 'interface-annotated interactome:' )
    print( '%d PPIs' % len(annotatedInteractome) )
    print( '%d proteins' % len(set(annotatedInteractome[["Protein_1", "Protein_2"]].values.flatten())) )
    
    write_interface_annotated_interactome_to_excel (
                    annotatedInteractome,
                    interactomeOutDir / ( interactome_name + '_structural_interactome.xlsx' ),
                    sheet_name = interactome_name + '_structural_interactome')
    
    print( '\n' + 'Plotting structural interactome' )
    edges = []
    for _, row in annotatedInteractome.iterrows():
        edges.append((row.Protein_1, row.Protein_2))
    network_plot(edges,
                 show = showFigs,
                 figdir = figDir,
                 figname = 'structural_interactome')
    
    #------------------------------------------------------------------------------------
    # further filter mutations
    #------------------------------------------------------------------------------------
    
    # load processed mutations
    print( '\n' + 'Reading processed mutations' )
    naturalMutations = pd.read_table(dbsnpMutationsFile, sep='\t')
    diseaseMutations = pd.read_table(hgmdMutationsFile, sep='\t')
    
    naturalSequenceMatch = sum(naturalMutations["seq_match"])
    diseaseSequenceMatch = sum(diseaseMutations["Seq_match"])
    print( '\n' + '%d non-disease mutations matching to sequence (%.1f %%)'
            % (naturalSequenceMatch,
               100 * naturalSequenceMatch / len(naturalMutations)) )
    print( '%d disease mutations matching to sequence (%.1f %%)'
            % (diseaseSequenceMatch,
               100 * diseaseSequenceMatch / len(diseaseMutations)) )
    
    # keep only mutations whose flanking sequence matches to UniProt sequence at same 
    # position reported by mutation database
    naturalMutations = naturalMutations[naturalMutations["seq_match"]].reset_index(drop=True)
    diseaseMutations = diseaseMutations[diseaseMutations["Seq_match"]].reset_index(drop=True)
    
    # rename columns
    naturalMutations.rename(columns={"protein":"Protein",
                                     "residue":"Mut_res",
                                     "aa_position":"Mutation_Position",
                                     "context":"Mut_context",
                                     "mut_context_pos":"Mut_context_pos",
                                     "allele": "Frequency_reporting_allele",
                                     "allele.1": "Variation_allele"}, inplace=True)
    diseaseMutations.rename(columns={"Codon_number":"Mutation_Position"}, inplace=True)
    diseaseMutations["Mut_res"] = diseaseMutations.apply(lambda x:
                                                         x["Mut_context"][x["Mut_context_pos"] - 1],
                                                         axis=1)
    
    # identify common mutations overlapping in position with disease mutations
    df = diseaseMutations[ ["Protein", "Mutation_Position"] ].copy()
    df2 = naturalMutations[ ["Protein", "Mutation_Position"] ].copy()
    df2 = df2.drop_duplicates().reset_index(drop=True)
    df = df.append(df2, ignore_index=True)
    duplicates = df.duplicated(keep='first')[ len( diseaseMutations ) : ]
    duplicates = duplicates.reset_index(drop=True)
    print( '\n' + '%d unique non-disease variants overlap in location with disease mutations' % sum( duplicates ) )
    
    # remove common mutations overlapping in position with disease mutations
    for i, row in df2[ duplicates ].iterrows():
        todrop = ( ( naturalMutations["Protein"] == row.Protein ) & 
                   ( naturalMutations["Mutation_Position"] == row.Mutation_Position ) )
        naturalMutations = naturalMutations[ todrop == False ]
    print( 'Overlaping non-disease mutations removed' )
    
    numNaturalMutations = len( naturalMutations )
    numDiseaseMutations = len( diseaseMutations )
    
    # remove invalid mutations, i.e., those with mutation residue similar to wild type
    naturalMutations = naturalMutations[naturalMutations.apply(lambda x: 
                                                               x["Mut_context"][x["Mut_context_pos"]-1] != x["Mut_res"],
                                                               axis=1)]
    diseaseMutations = diseaseMutations[diseaseMutations.apply(lambda x: 
                                                               x["WT_context"][x["Mut_context_pos"]-1] != x["Mut_res"],
                                                               axis=1)]
    
    print( '\n' + 'Number of invalid mutations removed (WT = Mut)' )
    print( 'non-disease invalid mutations: %d' % ( numNaturalMutations - len( naturalMutations ) ) )
    print( 'disease invalid mutations: %d' % ( numDiseaseMutations - len( diseaseMutations ) ) )
    
    numNaturalMutations = len(naturalMutations)
    numDiseaseMutations = len(diseaseMutations)
    
    print( '\n' + 'Number of mutations after removing invalid mutations:' )
    print( 'non-disease mutations: %d' % numNaturalMutations )
    print( 'disease mutations: %d' % numDiseaseMutations )
    
    # remove duplicate mutations, by position and mutant residue
    naturalMutations = naturalMutations.drop_duplicates(subset=["Protein",
                                                                "Mutation_Position",
                                                                "Mut_res"]).reset_index(drop=True)
    diseaseMutations = diseaseMutations.drop_duplicates(subset=["Protein",
                                                                "Mutation_Position",
                                                                "Mut_res"]).reset_index(drop=True)
    
    numNaturalMutations = len(naturalMutations)
    numDiseaseMutations = len(diseaseMutations)
    
    print( '\n' + 'Number of mutations after removing duplicates by position and mutant residue:' )
    print( 'non-disease mutations: %d' % numNaturalMutations )
    print( 'disease mutations: %d' % numDiseaseMutations )
    
    # plot fractions of disease and non-disease (common) mutations
    pie_plot([numNaturalMutations, numDiseaseMutations],
             labels=['Non-disease\nmutations', 'Disease\nmutations'],
             angle=90,
             colors = ['turquoise', 'orangered'],
             show = showFigs,
             figdir = figDir,
             figname = 'mutations')
    
    #------------------------------------------------------------------------------------
    # Calculate substitution score for all mutations, including variants at same position
    #------------------------------------------------------------------------------------
    
    # substitution matrix name
    matrixName = 'PAM30'
    
    # number of resamplings for bootstrap test on average substitution scores for common
    # and disease mutations
    bootstrap_test_itr = 10000
    
    subsMatrixFile = outDir / (matrixName + '.pkl')
    if not subsMatrixFile.is_file():
        produce_substitution_matrix(matrixName,
                                    subsMatrixFile)
    with open(subsMatrixFile, 'rb') as f:
        subsTable = pickle.load(f)
    naturalChange = naturalMutations.apply(lambda x:
                                           (x["Mut_context"][x["Mut_context_pos"]-1],
                                            x["Mut_res"]),
                                           axis=1)
    naturalScore = naturalChange.apply(lambda x: subsTable[x])
    print( '\n' + 'Avg. score for natural mutations: %.1f (SE = %g, n = %d)' % (np.mean(naturalScore),
                                                                                sderror(naturalScore),
                                                                                len(naturalScore)) )
    
    diseaseChange = diseaseMutations.apply(lambda x:
                                           (x["WT_context"][x["Mut_context_pos"]-1],
                                            x["Mut_res"]),
                                           axis=1)
    diseaseScore = diseaseChange.apply(lambda x: subsTable[x])
    print( 'Avg. score for disease mutations: %.1f (SE = %g, n = %d)' % (np.mean(diseaseScore),
                                                                         sderror(diseaseScore),
                                                                         len(diseaseScore)) )
    
    t_test( naturalScore, diseaseScore )
    bootstrap_test( naturalScore, diseaseScore, bootstrap_test_itr )
    
    box_plot([naturalScore, diseaseScore],
             xlabels = ('Non-disease\nmutations','Disease\nmutations'),
             ylabel = ('PAM30 substitution score'),
             ylim = [-20, 5],
             colors = ['turquoise', 'orangered'],
             fontsize = 26,
             show = showFigs,
             figdir = figDir,
             figname = 'substitution_score')
    
    #------------------------------------------------------------------------------------
    # predict PPI perturbations
    #------------------------------------------------------------------------------------
    
    # Consider PPI perturbations only for PPIs with this maximum number of interfaces 
    # mapped from distinct PDB binding chain pairs.
    # This parameter is irrelevent if the flag "merge_interfaces" is True.
    # Set to inf for unlimited number of interfaces.
    maxInterfaces = np.inf
    
    # Predict PPI perturbation if mutation is this number of residues away from any
    # interface residue.
    # Set to 0 if mutation must be exactly at interface residue.
    numResFromInterface = 0
    
    # Consider PPI perturbations only for PPIs with this minimum number of partners
    minPartners = 1
    
    mutationPerturbFile = interactomeOutDir / 'mutation_perturbs.pkl'
    if mutationPerturbFile.is_file():
        print( '\n' + 'Loading geometry-based PPI perturbation predictions' )
        with open(mutationPerturbFile, 'rb') as f:
            naturalMutations_perturb, diseaseMutations_perturb = pickle.load(f)
    else:
        print( '\n' + 'Predicting PPI perturbations based on geometry' )
        naturalMutations_perturb = mutation_PPI_interface_perturbations(naturalMutations,
                                                                        annotatedInteractome,
                                                                        maxInterfaces,
                                                                        numResFromInterface)
        diseaseMutations_perturb = mutation_PPI_interface_perturbations(diseaseMutations,
                                                                        annotatedInteractome,
                                                                        maxInterfaces,
                                                                        numResFromInterface)
        with open(mutationPerturbFile, 'wb') as fOut:
            pickle.dump([naturalMutations_perturb, diseaseMutations_perturb], fOut)
    
    naturalPerturbs = naturalMutations.copy()
    naturalPerturbs["partners"] = naturalMutations_perturb.apply(lambda x: x[0])
    naturalPerturbs["perturbations"] = naturalMutations_perturb.apply(lambda x: x[1])
    naturalPerturbs = naturalPerturbs[naturalPerturbs["partners"].apply(len) 
                                      >= minPartners].reset_index(drop=True)
#     naturalPerturbs = naturalPerturbs[naturalPerturbs["perturbations"].apply(lambda x: sum(x > 0) < 2)].reset_index(drop=True)
    
    diseasePerturbs = diseaseMutations.copy()
    diseasePerturbs["partners"] = diseaseMutations_perturb.apply(lambda x: x[0])
    diseasePerturbs["perturbations"] = diseaseMutations_perturb.apply(lambda x: x[1])
    diseasePerturbs = diseasePerturbs[diseasePerturbs["partners"].apply(len) 
                                      >= minPartners].reset_index(drop=True)
#     diseasePerturbs = diseasePerturbs[diseasePerturbs["perturbations"].apply(lambda x: sum(x > 0) < 2)].reset_index(drop=True)   
    
    # drop duplicate mutations based on location, regardless of residue type 
    naturalPerturbs = naturalPerturbs.drop_duplicates(subset=["Protein",
                                                              "Mutation_Position"]).reset_index(drop=True)
    diseasePerturbs = diseasePerturbs.drop_duplicates(subset=["Protein",
                                                              "Mutation_Position"]).reset_index(drop=True)
    
    print( '\n' + 'Number of mutations with PPI perturbation predictions after removing duplicates by position' )
    print( 'non-disease: %d' % len( naturalPerturbs ) )
    print( 'disease: %d' % len( diseasePerturbs ) )
    
    filteredMutationPerturbFile = interactomeOutDir / 'filtered_mutation_perturbs.pkl'
    with open(filteredMutationPerturbFile, 'wb') as fOut:
        pickle.dump([naturalPerturbs, diseasePerturbs], fOut)
    
    natural_proteins = set( naturalPerturbs["Protein"] )
    disease_proteins = set( diseasePerturbs["Protein"] )
    
    # average degree for proteins carrying non-disease and disease mutations in
    # the reference interactome
    natural_degree = []
    disease_degree = []
    for p in natural_proteins:
        natural_degree.append( sum( (interactome[["Protein_1", "Protein_2"]] == p).any(1) ) )
    for p in disease_proteins:
        disease_degree.append( sum( (interactome[["Protein_1", "Protein_2"]] == p).any(1) ) )
    print( '\n' + 'Avg degree for proteins carrying mutations with predicted PPI perturbations:' )
    print( 'Reference interactome:')
    print( 'proteins carrying non-disease mutations: %.1f' % np.mean( natural_degree ) )
    print( 'proteins carrying disease mutations: %.1f' % np.mean( disease_degree ) )
    
    # average degree for proteins carrying non-disease and disease mutations in
    # the structural interactome
    natural_degree = []
    disease_degree = []
    for p in natural_proteins:
        natural_degree.append( sum( (annotatedInteractome[["Protein_1", "Protein_2"]] == p).any(1) ) )
    for p in disease_proteins:
        disease_degree.append( sum( (annotatedInteractome[["Protein_1", "Protein_2"]] == p).any(1) ) )
    print( 'Structural interactome:' )
    print( 'proteins carrying non-disease mutations: %.1f' % np.mean( natural_degree ) )
    print( 'proteins carrying disease mutations: %.1f' % np.mean( disease_degree ) )
    
    # Identify interfacial and non-interfacial mutations
    natural_interfacial = naturalPerturbs["perturbations"].apply(lambda x: any(x > 0))
    disease_interfacial = diseasePerturbs["perturbations"].apply(lambda x: any(x > 0))
    
    numNatural_interfacial = sum( natural_interfacial == True )
    numNatural_noninterfacial = sum( natural_interfacial == False )
    numDisease_interfacial = sum( disease_interfacial == True )
    numDisease_noninterfacial = sum( disease_interfacial == False )
    
    numNatural_considered = numNatural_interfacial + numNatural_noninterfacial
    numDisease_considered = numDisease_interfacial + numDisease_noninterfacial
    
    fracNatural_interfacial = numNatural_interfacial / numNatural_considered
    fracDisease_interfacial = numDisease_interfacial / numDisease_considered
    fracNatural_error = sderror_on_fraction(numNatural_interfacial, numNatural_considered)
    fracDisease_error = sderror_on_fraction(numDisease_interfacial, numDisease_considered)
    
    print( '\n' + 'Fraction of interfacial mutations:' )
    print( 'Non-disease: %.3f (SE = %g, %d out of %d)' % (fracNatural_interfacial,
                                                          fracNatural_error,
                                                          numNatural_interfacial,
                                                          numNatural_considered) )
    
    print( 'Disease: %.3f %% (SE = %g, %d out of %d)' % (fracDisease_interfacial,
                                                         fracDisease_error,
                                                         numDisease_interfacial,
                                                         numDisease_considered) )
    
    fisher_test([numNatural_interfacial, numNatural_noninterfacial],
                [numDisease_interfacial, numDisease_noninterfacial])
    
    bar_plot([fracNatural_interfacial, fracDisease_interfacial],
             error = [fracNatural_error, fracDisease_error],
             xlabels = ('Non-disease\nmutations','Disease\nmutations'),
             ylabel = ('Fraction on PPI interface'),
             colors = ['turquoise', 'orangered'],
             ylim = [0, 0.2],
             barwidth = 0.6,
             fontsize = 24,
             show = showFigs,
             figdir = figDir,
             figname = 'Fraction_interfacial_mut')
    
    # write to file interfacial non-disease mutations and perturbed protein pair after mapping 
    # mutation onto PDB chains
    natural_mutations_onchains_file = interactomeOutDir / 'nondisease_mutations_onchains.txt'
    if not natural_mutations_onchains_file.is_file():
        print( '\n' + 'Writing interfacial non-disease mutations mapped onto structures to file ' 
                + str( natural_mutations_onchains_file ) )
        write_pdb_mapped_mutations (naturalPerturbs[ natural_interfacial ],
                                    interfaceAnnotatedInteractomeFile,
                                    chainMapFile,
                                    chainSeqFile,
                                    ProteinSeqFile,
                                    chainResOrderFile,
                                    pdbDir,
                                    natural_mutations_onchains_file)
    
    # write to file interfacial disease mutations and perturbed protein pairs after mapping 
    # mutation onto PDB chains
    disease_mutations_onchains_file = interactomeOutDir / 'disease_mutations_onchains.txt'
    if not disease_mutations_onchains_file.is_file():
        print( '\n' + 'Writing interfacial disease mutations mapped onto structures to file ' 
                + str( disease_mutations_onchains_file ) )
        write_pdb_mapped_mutations (diseasePerturbs[disease_interfacial],
                                    interfaceAnnotatedInteractomeFile,
                                    chainMapFile,
                                    chainSeqFile,
                                    ProteinSeqFile,
                                    chainResOrderFile,
                                    pdbDir,
                                    disease_mutations_onchains_file)
    
    print( '\n' + 'Creating network perturbed by non-disease mutations' )
    nodes, edges, nodeColors, edgeColors = create_perturbed_network (
                                            annotatedInteractome,
                                            naturalPerturbs,
                                            interactomeOutDir / 'nondiseaseMut_perturbed_network',
                                            interactomeOutDir / 'nondiseaseMut_node_colors')
    network_plot (edges,
                  nodes = nodes,
                  nodeSizes = [20] * len(nodes),
                  edgeWidth = 1,
                  nodeColors = nodeColors,
                  edgeColors = edgeColors,
                  show = showFigs,
                  figdir = figDir,
                  figname = 'nondisease_perturbed_interactome')
    
    print( '\n' + 'Creating network perturbed by disease mutations' )
    nodes, edges, nodeColors, edgeColors = create_perturbed_network (
                                            annotatedInteractome,
                                            diseasePerturbs,
                                            interactomeOutDir / 'diseaseMut_perturbed_network',
                                            interactomeOutDir / 'diseaseMut_node_colors')
    network_plot (edges,
                  nodes = nodes,
                  nodeSizes = [20] * len(nodes),
                  edgeWidth = 1,
                  nodeColors = nodeColors,
                  edgeColors = edgeColors,
                  show = showFigs,
                  figdir = figDir,
                  figname = 'disease_perturbed_interactome')
    
    #------------------------------------------------------------------------------------
    # Fraction of predicted edgetic mutations
    #------------------------------------------------------------------------------------
    
    print( '\n' + 'Labeling mutation edgotypes' )
    naturalPerturbs["Edgotype"] = naturalPerturbs["perturbations"].apply(lambda x:
                                                                         'Edgetic' if sum(x > 0) > 0
                                                                         else 'Non-edgetic')
    diseasePerturbs["Edgotype"] = diseasePerturbs["perturbations"].apply(lambda x:
                                                                         'Edgetic' if sum(x > 0) > 0
                                                                         else 'Non-edgetic')
    
    # print predicted interactome perturbations by non-disease mutations to excel file
    naturalPerturbs_output = naturalPerturbs.copy()
    naturalPerturbs_output["partners"] = naturalPerturbs_output["partners"].apply(lambda x: ','.join(x))
    naturalPerturbs_output["perturbations"] = ( naturalPerturbs_output["perturbations"]
                                                .apply(lambda x: ','.join(map(str, [int(np.ceil(k)) for k in x]))) )
    naturalPerturbs_output.rename( columns = {"Mut_context":"WT_context"}, inplace=True )
    writer = pd.ExcelWriter( str( interactomeOutDir / (interactome_name + '_mutation_edgotypes.xlsx') ) )
    naturalPerturbs_output.to_excel (writer,
                                     sheet_name = 'Nondisease_mutations',
                                     columns = ["gene",
                                                "Protein",
                                                "ID",
                                                "fxn-class",
                                                "Mutation_Position",
                                                "WT_context",
                                                "Mut_context_pos",
                                                "Mut_res",
                                                "partners",
                                                "perturbations",
                                                "Edgotype",
                                                "MAF",
                                                "allele_origin",
                                                "alleles",
                                                "Frequency_reporting_allele",
                                                "count",
                                                "Variation_allele",
                                                "validated",
                                                "validation",
                                                "method",
                                                "assertion",
                                                "1000genome",
                                                "date",
                                                "het",
                                                "se(het)",
                                                "max_prob",
                                                "min_prob",
                                                "frame",
                                                "locus_id",
                                                "mrna_acc",
                                                "prot_acc"])
    
    # print predicted interactome perturbations by disease mutations to excel file
    diseasePerturbs_output = diseasePerturbs.copy()
    diseasePerturbs_output["partners"] = diseasePerturbs_output["partners"].apply(lambda x: ','.join(x))
    diseasePerturbs_output["perturbations"] = ( diseasePerturbs_output["perturbations"]
                                                .apply(lambda x: ','.join(map(str, [int(np.ceil(k)) for k in x]))) )
    diseasePerturbs_output.to_excel (writer,
                                     sheet_name = 'Disease_mutations',
                                     columns = ["Gene",
                                                "entrezid",
                                                "Protein",
                                                "ACC_NUM",
                                                "dbsnp",
                                                "Variant_class",
                                                "Mutation_Position",
                                                "Disease",
                                                "WT_context",
                                                "Mut_context",
                                                "Mut_context_pos",
                                                "partners",
                                                "perturbations",
                                                "Edgotype"])
    writer.save()
    
    numNaturalMut_edgetic = sum( naturalPerturbs["Edgotype"] == 'Edgetic' )
    numNaturalMut_nonedgetic = sum( naturalPerturbs["Edgotype"] == 'Non-edgetic' )
    numDiseaseMut_edgetic = sum( diseasePerturbs["Edgotype"] == 'Edgetic' )
    numDiseaseMut_nonedgetic = sum( diseasePerturbs["Edgotype"] == 'Non-edgetic' )
    
    numNaturalMut_considered = numNaturalMut_edgetic + numNaturalMut_nonedgetic
    numDiseaseMut_considered = numDiseaseMut_edgetic + numDiseaseMut_nonedgetic
    
    fracNaturalMut_edgetic = numNaturalMut_edgetic / numNaturalMut_considered
    fracDiseaseMut_edgetic = numDiseaseMut_edgetic / numDiseaseMut_considered
    fracNaturalMut_error = sderror_on_fraction(numNaturalMut_edgetic, numNaturalMut_considered)
    fracDiseaseMut_error = sderror_on_fraction(numDiseaseMut_edgetic, numDiseaseMut_considered)
    
    print( '\n' + 'Fraction of predicted edgetic mutations:' )
    print( 'Non-disease mutations: %.3f (SE = %g, %d out of %d)' % (fracNaturalMut_edgetic,
                                                                    fracNaturalMut_error,
                                                                    numNaturalMut_edgetic,
                                                                    numNaturalMut_considered) )
    
    print( 'Disease mutations: %.3f (SE = %g, %d out of %d)' % (fracDiseaseMut_edgetic,
                                                                fracDiseaseMut_error,
                                                                numDiseaseMut_edgetic,
                                                                numDiseaseMut_considered) )
    
    fisher_test([numNaturalMut_edgetic, numNaturalMut_nonedgetic],
                [numDiseaseMut_edgetic, numDiseaseMut_nonedgetic])
    
    pie_plot([numNaturalMut_nonedgetic, numNaturalMut_edgetic],
             angle = 90,
             colors = ['thistle', 'mediumslateblue'],
             show = showFigs,
             figdir = figDir,
             figname = 'non_disease_edgetic_mutations_geometrybased')
    pie_plot([numDiseaseMut_nonedgetic, numDiseaseMut_edgetic],
             angle=90,
             colors = ['thistle', 'mediumslateblue'],
             show = showFigs,
             figdir = figDir,
             figname = 'disease_edgetic_mutations_geometrybased')
    
    #------------------------------------------------------------------------------------
    # Fraction of predicted mono-edgetic mutations
    #------------------------------------------------------------------------------------
    
    originalNumNatural = len( naturalPerturbs )
    originalNumDisease = len( diseasePerturbs )
    
    naturalPerturbNum = naturalPerturbs["perturbations"].apply(lambda x: sum( x > 0 ))
    diseasePerturbNum = diseasePerturbs["perturbations"].apply(lambda x: sum( x > 0 ))
    
    multi_histogram_plot ([naturalPerturbNum[ naturalPerturbNum > 0 ],
                           diseasePerturbNum[ diseasePerturbNum > 0 ]],
                          ['g', 'r'],
                          xlabel = 'Number of PPIs perturbed',
                          ylabel = 'Frequency among edgetic mutations',
                          edgecolor = 'k',
                          fontsize = 14,
                          bins = 20,
                          alpha = 0.1,
                          leg = ['Non-disease', 'Disease'],
                          show = showFigs,
                          figdir = figDir,
                          figname = 'prediction_numPerturbs_geometrybased')
    
    naturalPerturbs["single_perturbation"] = naturalPerturbs["perturbations"].apply(lambda x: sum(x > 0) == 1)
    diseasePerturbs["single_perturbation"] = diseasePerturbs["perturbations"].apply(lambda x: sum(x > 0) == 1)
    
    print( '\n' + 'Fraction of predicted mono-edgetic mutations among edgetic mutations:' )
    print( 'non-disease mutation: %.3f (%d out of %d)' 
            % ( sum( naturalPerturbs["single_perturbation"] ) / sum( naturalPerturbs["Edgotype"] == 'Edgetic'),
                sum( naturalPerturbs["single_perturbation"] ),
                sum( naturalPerturbs["Edgotype"] == 'Edgetic') ) )
    print( 'disease mutation: %.3f (%d out of %d)' 
            % ( sum( diseasePerturbs["single_perturbation"] ) / sum( diseasePerturbs["Edgotype"] == 'Edgetic'),
                sum( diseasePerturbs["single_perturbation"] ),
                sum( diseasePerturbs["Edgotype"] == 'Edgetic') ) )
    
    numNatural_edgetic = sum( naturalPerturbs["single_perturbation"] )
    numNatural_nonedgetic = sum( naturalPerturbs["single_perturbation"] == False )
    numDisease_edgetic = sum( diseasePerturbs["single_perturbation"] )
    numDisease_nonedgetic = sum( diseasePerturbs["single_perturbation"] == False )
    
    numNatural_considered = numNatural_edgetic + numNatural_nonedgetic
    numDisease_considered = numDisease_edgetic + numDisease_nonedgetic
    
    fracNatural_edgetic = numNatural_edgetic / numNatural_considered
    fracDisease_edgetic = numDisease_edgetic / numDisease_considered
    fracNatural_error = sderror_on_fraction(numNatural_edgetic, numNatural_considered)
    fracDisease_error = sderror_on_fraction(numDisease_edgetic, numDisease_considered)
    
    print( '\n' + 'Fraction of mono-edgetic mutations:' )
    print( 'Non-disease mutations: %.3f (SE = %g, %d out of %d)' % (fracNatural_edgetic,
                                                                    fracNatural_error,
                                                                    numNatural_edgetic,
                                                                    numNatural_considered) )
    
    print( 'Disease mutations: %.3f (SE = %g, %d out of %d)' % (fracDisease_edgetic,
                                                                fracDisease_error,
                                                                numDisease_edgetic,
                                                                numDisease_considered) )
    
    fisher_test([numNatural_edgetic, numNatural_nonedgetic],
                [numDisease_edgetic, numDisease_nonedgetic])
    
    #labels = ['n = ' + str(numNatural_nonedgetic), 'n = ' + str(numNatural_edgetic)]
    pie_plot([numNatural_nonedgetic, numNatural_edgetic],
             angle = 90,
             colors = ['thistle', 'mediumslateblue'],
             show = showFigs,
             figdir = figDir,
             figname = 'non_disease_monoedgetic_mutations_geometrybased')
    pie_plot([numDisease_nonedgetic, numDisease_edgetic],
             angle=90,
             colors = ['thistle', 'mediumslateblue'],
             show = showFigs,
             figdir = figDir,
             figname = 'disease_monoedgetic_mutations_geometrybased')
    
    #------------------------------------------------------------------------------------
    # Fraction of mutation-targeted PPIs with ∆∆G exceeding a specified cutoff
    #------------------------------------------------------------------------------------
    
    # Minimum reduction in binding free energy DDG required for interaction perturbation.
    ddgCutoff = 0.5
    
    # read loss in binding free energy for interfacial mutations mapped on PDB chains
    naturalMutationsDDG = read_mutation_ddg(interactomeOutDir / 'nondisease_mutations_onchains_ddg.txt')
    diseaseMutationsDDG = read_mutation_ddg(interactomeOutDir / 'disease_mutations_onchains_ddg.txt')
    
    writer = pd.ExcelWriter( str( interactomeOutDir / ( interactome_name + '_mutation_∆∆G.xlsx' ) ) )
    naturalMutationsDDG.to_excel(writer, sheet_name = 'nondisease_mutations_∆∆G')
    diseaseMutationsDDG.to_excel(writer, sheet_name = 'disease_mutations_∆∆G')
    writer.save()
    
    numNatural_ddg_considered = len( naturalMutationsDDG["DDG"] )
    numDisease_ddg_considered = len( diseaseMutationsDDG["DDG"] )
    
    print( '\n' + 'Avg. change in binding energy (∆∆G) for mutation-targeted PPIs:' )
    print( 'Non-disease: %.1f (SE = %g, n = %d)' % (np.mean( naturalMutationsDDG["DDG"] ),
                                                    sderror( naturalMutationsDDG["DDG"] ),
                                                    numNatural_ddg_considered) )
    
    print( 'Disease: %.1f (SE = %g, n = %d)' % (np.mean( diseaseMutationsDDG["DDG"] ),
                                                sderror( diseaseMutationsDDG["DDG"] ),
                                                numDisease_ddg_considered) )
    t_test( naturalMutationsDDG["DDG"], diseaseMutationsDDG["DDG"] )
    
    multi_histogram_plot ([diseaseMutationsDDG["DDG"], naturalMutationsDDG["DDG"]],
                          ['red', 'green'],
                          xlabel = 'Change in binding free energy (∆∆G)',
                          ylabel = 'Number of mutations',
                          leg = ['Disease interfacial mutations', 'Non-disease interfacial mutations'],
                          bins = 25,
                          alpha = 0.7,
                          fontsize = 24,
                          show = showFigs,
                          figdir = figDir,
                          figname = 'mut_ddg_histogram')
    
    numNatural_ddg = sum( naturalMutationsDDG["DDG"] > ddgCutoff )
    numDisease_ddg = sum( diseaseMutationsDDG["DDG"] > ddgCutoff )
    fracNatural_ddg = numNatural_ddg / numNatural_ddg_considered
    fracDisease_ddg = numDisease_ddg / numDisease_ddg_considered
    fracNatural_ddg_error = sderror_on_fraction( numNatural_ddg, numNatural_ddg_considered )
    fracDisease_ddg_error = sderror_on_fraction( numDisease_ddg, numDisease_ddg_considered )
    
    print( '\n' + 'Fraction of mutation-targeted PPIs with ∆∆G > %.1f kcal/mol:' % ddgCutoff )
    print( 'Non-disease: %.3f (SE = %g, ntot = %d)' % (fracNatural_ddg,
                                                       fracNatural_ddg_error,
                                                       numNatural_ddg_considered) )
    
    print( 'Disease: %.3f (SE = %g, ntot = %d)' % (fracDisease_ddg,
                                                   fracDisease_ddg_error,
                                                   numDisease_ddg_considered) )
    
    fisher_test([numNatural_ddg, numNatural_ddg_considered - numNatural_ddg],
                [numDisease_ddg, numDisease_ddg_considered - numDisease_ddg])
    
    bar_plot([fracNatural_ddg, fracDisease_ddg],
             error = [fracNatural_ddg_error, fracDisease_ddg_error],
             xlabels = ('PPIs with\nnon-disease\nmutations\nat interface',
                        'PPIs with\ndisease\nmutations\nat interface'),
             ylabel = ('Fraction with ∆∆G > %.1f kcal/mol' % ddgCutoff),
             ylabels = [0, 0.2, 0.4, 0.6, 0.8],
             ylim = [0, 0.8],
             colors = ['turquoise', 'orangered'],
             barwidth = 0.6,
             fontsize = 24,
             show = showFigs,
             figdir = figDir,
             figname = 'Sample_fraction_ddg_>%.1f' % ddgCutoff)
    
    #------------------------------------------------------------------------------------
    # Fractions of predicted edgetic and mono-edgetic mutations using physics-based 
    # prediction of PPI perturbation (i.e., using binding energy change)
    #------------------------------------------------------------------------------------
    
    ddgMutationPerturbFile = interactomeOutDir / 'ddg_mutation_perturbs.pkl'
    if ddgMutationPerturbFile.is_file():
        print( '\n' + 'Loading physics-based mutation edgotype predictions' )
        with open(ddgMutationPerturbFile, 'rb') as f:
            ( fracNaturalMut_edgetic_ddg,
              fracNatural_edgetic_ddg,
              fracDiseaseMut_edgetic_ddg,
              fracDisease_edgetic_ddg ) = pickle.load(f)
        print( 'Fraction of predicted non-disease edgetic mutations: %f' % fracNaturalMut_edgetic_ddg )
        print( 'Fraction of predicted disease edgetic mutations: %f' % fracDiseaseMut_edgetic_ddg )
        print( 'Fraction of predicted non-disease mono-edgetic mutations: %f' % fracNatural_edgetic_ddg )
        print( 'Fraction of predicted disease mono-edgetic mutations: %f' % fracDisease_edgetic_ddg )
    else:
        numItr = 1000
        print( '\n' + 'Performing physics-based edgotype prediction for non-disease mutations' )
        perturbations, knownDDG, unknownDDG = energy_based_perturb_prediction (naturalPerturbs,
                                                                               naturalMutationsDDG,
                                                                               ddgCutoff)
        if unknownDDG > 0:
            print( '\n' + 'Perturbation probability will be used to re-predict perturbation ' +
                   'for %d PPIs with unknown ∆∆G out of %d' % (unknownDDG, knownDDG + unknownDDG) )
            fracEdgetic = []
            fracSingleEdgetic = []
            for i in range( numItr ):
                perturbations, knownDDG, unknownDDG = energy_based_perturb_prediction (naturalPerturbs,
                                                                                       naturalMutationsDDG,
                                                                                       ddgCutoff)
                edgetic = perturbations.apply( lambda x: sum(x > 0) > 0 )
                singleEdgetic = perturbations.apply( lambda x: sum(x > 0) == 1 )
                fracEdgetic.append( sum( edgetic ) / len( edgetic ) )
                fracSingleEdgetic.append( sum( singleEdgetic ) / len( singleEdgetic ) )
            fracNaturalMut_edgetic_ddg = np.mean( fracEdgetic )
            fracNatural_edgetic_ddg = np.mean( fracSingleEdgetic )
            print( 'Avg. fraction of predicted non-disease edgetic mutations: %f (SEM = %g)' 
                    % (fracNaturalMut_edgetic_ddg, sderror( fracEdgetic )) )
            print( 'Avg. fraction of predicted non-disease mono-edgetic mutations: %f (SEM = %g)' 
                    % (fracNatural_edgetic_ddg, sderror( fracSingleEdgetic )) )
        else:
            print( '∆∆G known for all PPIs' )
            edgetic = perturbations.apply( lambda x: sum(x > 0) > 0 )
            singleEdgetic = perturbations.apply( lambda x: sum(x > 0) == 1 )
            fracNaturalMut_edgetic_ddg = sum( edgetic ) / len( edgetic )
            fracNatural_edgetic_ddg = sum( singleEdgetic ) / len( singleEdgetic )
            print( 'Fraction of predicted non-disease edgetic mutations: %f (SEF = %g, ntot = %d)' 
                    % (fracNaturalMut_edgetic_ddg,
                       sderror_on_fraction( sum( edgetic ), len( edgetic ) ),
                       len( edgetic ) ) )
            print( 'Fraction of predicted non-disease mono-edgetic mutations: %f (SEF = %g, ntot = %d)' 
                    % (fracNatural_edgetic_ddg,
                       sderror_on_fraction( sum( singleEdgetic ), len( singleEdgetic ) ),
                       len( singleEdgetic ) ) )
    
        print( '\n' + 'Performing physics-based edgotype prediction for disease mutations' )
        perturbations, knownDDG, unknownDDG = energy_based_perturb_prediction (diseasePerturbs,
                                                                               diseaseMutationsDDG,
                                                                               ddgCutoff)
        if unknownDDG > 0:
            print( '\n' + 'Perturbation probability will be used to re-predict perturbation ' +
                   'for %d PPIs with unknown ∆∆G out of %d' % (unknownDDG, knownDDG + unknownDDG) )
            fracEdgetic = []
            fracSingleEdgetic = []
            for i in range( numItr ):
                perturbations, knownDDG, unknownDDG = energy_based_perturb_prediction (diseasePerturbs,
                                                                                       diseaseMutationsDDG,
                                                                                       ddgCutoff)
                edgetic = perturbations.apply( lambda x: sum(x > 0) > 0 )
                singleEdgetic = perturbations.apply( lambda x: sum(x > 0) == 1 )
                fracEdgetic.append( sum( edgetic ) / len( edgetic ) )
                fracSingleEdgetic.append( sum( singleEdgetic ) / len( singleEdgetic ) )
            fracDiseaseMut_edgetic_ddg = np.mean( fracEdgetic )
            fracDisease_edgetic_ddg = np.mean( fracSingleEdgetic )
            print( 'Avg. fraction of predicted disease edgetic mutations: %f (SEM = %g)' 
                    % (fracDiseaseMut_edgetic_ddg, sderror( fracEdgetic )) )
            print( 'Avg. fraction of predicted disease mono-edgetic mutations: %f (SEM = %g)' 
                    % (fracDisease_edgetic_ddg, sderror( fracSingleEdgetic )) )
        else:
            print( '∆∆G known for all PPIs' )
            edgetic = perturbations.apply( lambda x: sum(x > 0) > 0 )
            singleEdgetic = perturbations.apply( lambda x: sum(x > 0) == 1 )
            fracDiseaseMut_edgetic_ddg = sum( edgetic ) / len( edgetic )
            fracDisease_edgetic_ddg = sum( singleEdgetic ) / len( singleEdgetic )
            print( 'Fraction of predicted disease edgetic mutations: %f (SEF = %g, ntot = %d)' 
                    % (fracDiseaseMut_edgetic_ddg,
                       sderror_on_fraction( sum( edgetic ), len( edgetic ) ),
                       len( edgetic ) ) )
            print( 'Fraction of predicted disease mono-edgetic mutations: %f (SEF = %g, ntot = %d)' 
                    % (fracDisease_edgetic_ddg,
                       sderror_on_fraction( sum( singleEdgetic ), len( singleEdgetic ) ),
                       len( singleEdgetic ) ) )
        with open(ddgMutationPerturbFile, 'wb') as fOut:
            pickle.dump( [ fracNaturalMut_edgetic_ddg,
                           fracNatural_edgetic_ddg,
                           fracDiseaseMut_edgetic_ddg,
                           fracDisease_edgetic_ddg ], fOut )
            
    #------------------------------------------------------------------------------------
    # Load experimental edgotype data
    #------------------------------------------------------------------------------------
    
    expMutations = pd.read_excel( inDir / "Sahni_2015_Table_S3.xlsx",
                                  sheet_name = 'Table S3C' )
    expMutations_ppi = pd.read_excel( inDir / "Sahni_2015_Table_S3.xlsx",
                                      sheet_name = 'Table S3A' )
    
    expMutations = expMutations[ expMutations["Edgotype_class"]
                                 .apply(lambda x: x in {'Edgetic',
                                                        'Quasi-null',
                                                        'Quasi-wild-type'})
                               ].reset_index( drop = True )
    
    expMutations["partners"] = expMutations["Allele_ID"].apply(
        lambda x: expMutations_ppi.loc[ expMutations_ppi["Allele_ID"] == x,
                                        "Interactor_Gene_ID" ].values )
    expMutations["Mut_interaction"] = expMutations["Allele_ID"].apply(
        lambda x: expMutations_ppi.loc[ expMutations_ppi["Allele_ID"] == x,
                                        "Y2H_score" ].values )
    
    expMutations["WT_interaction"] = expMutations.apply(
        lambda x: np.array( [ expMutations_ppi.loc[ (expMutations_ppi["Category"] == 'Wild-type') & 
                                                    (expMutations_ppi["Entrez_Gene_ID"] == x["Entrez_Gene_ID"]) & 
                                                    (expMutations_ppi["Interactor_Gene_ID"] == p),
                                                    "Y2H_score"].item()
                              for p in x["partners"] ] ),
        axis=1 )
    
    expMutations["perturbations"] = expMutations.apply(
        lambda x:  0 + ( ( x["WT_interaction"] == 1 ) & ( x["Mut_interaction"] == 0 ) ), axis=1 )
    
    expNaturalMutations = expMutations[ expMutations["Category"]
                                        == 'Non-disease variant' ].reset_index( drop = True )
    expDiseaseMutations = expMutations[ expMutations["Category"]
                                        == 'Disease mutation' ].reset_index( drop = True )
    
    #------------------------------------------------------------------------------------
    # Fraction of edgetic mutations in experiments
    #------------------------------------------------------------------------------------
    
    if consider_exp_QuasiNull_perturbs:
        # Number of edgetic (E) and quasi-null (Q) mutations in experiment
        numNaturalMut_edgetic_exp = sum( (expNaturalMutations["Edgotype_class"] == 'Edgetic') |
                                         (expNaturalMutations["Edgotype_class"] == 'Quasi-null') )
        numDiseaseMut_edgetic_exp = sum( (expDiseaseMutations["Edgotype_class"] == 'Edgetic') |
                                         (expDiseaseMutations["Edgotype_class"] == 'Quasi-null') )
        
        # Number of non-edgetic mutations (quasi-wild-type) in experiment
        numNaturalMut_nonedgetic_exp = sum( expNaturalMutations["Edgotype_class"] == 'Quasi-wild-type' )
        numDiseaseMut_nonedgetic_exp = sum( expDiseaseMutations["Edgotype_class"] == 'Quasi-wild-type' )
    else:
        # Number of edgetic (E) mutations in experiment
        numNaturalMut_edgetic_exp = sum( expNaturalMutations["Edgotype_class"] == 'Edgetic' )
        numDiseaseMut_edgetic_exp = sum( expDiseaseMutations["Edgotype_class"] == 'Edgetic' )
        
        # Number of non-edgetic mutations (quasi-wild-type or quasi-null) in experiment
        numNaturalMut_nonedgetic_exp = sum( (expNaturalMutations["Edgotype_class"] == 'Quasi-null') |
                                            (expNaturalMutations["Edgotype_class"] == 'Quasi-wild-type') )
        numDiseaseMut_nonedgetic_exp = sum( (expDiseaseMutations["Edgotype_class"] == 'Quasi-null') |
                                            (expDiseaseMutations["Edgotype_class"] == 'Quasi-wild-type') )    
    # Total number of mutations tested in experiment
    numNaturalMut_considered_exp = numNaturalMut_edgetic_exp + numNaturalMut_nonedgetic_exp
    numDiseaseMut_considered_exp = numDiseaseMut_edgetic_exp + numDiseaseMut_nonedgetic_exp
    
    # Fraction of mutations that are edgetic (E) in experiment
    fracNaturalMut_edgetic_exp = numNaturalMut_edgetic_exp / numNaturalMut_considered_exp
    fracDiseaseMut_edgetic_exp = numDiseaseMut_edgetic_exp / numDiseaseMut_considered_exp
    
    # Standard error of the fraction of edgetic mutations in experiment
    fracNaturalMut_error_exp = sderror_on_fraction(numNaturalMut_edgetic_exp, numNaturalMut_considered_exp)
    fracDiseaseMut_error_exp = sderror_on_fraction(numDiseaseMut_edgetic_exp, numDiseaseMut_considered_exp)
    
    print( '\n' + 'Fraction of experimental edgetic mutations:' )
    print('Non-disease mutations: %.3f (SE = %g, %d out of %d)' % (fracNaturalMut_edgetic_exp,
                                                                   fracNaturalMut_error_exp,
                                                                   numNaturalMut_edgetic_exp,
                                                                   numNaturalMut_considered_exp))
    
    print('Disease mutations: %.3f (SE = %g, %d out of %d)' % (fracDiseaseMut_edgetic_exp,
                                                               fracDiseaseMut_error_exp,
                                                               numDiseaseMut_edgetic_exp,
                                                               numDiseaseMut_considered_exp))
    
    fisher_test([numNaturalMut_edgetic_exp, numNaturalMut_nonedgetic_exp],
                [numDiseaseMut_edgetic_exp, numDiseaseMut_nonedgetic_exp])
    
    pie_plot([numNaturalMut_nonedgetic_exp, numNaturalMut_edgetic_exp],
             angle = 90,
             #labels = ['Nonedgetic', 'Edgetic'],
             colors = ['thistle', 'mediumslateblue'],
             show = showFigs,
             figdir = figDir,
             figname = 'non_disease_edgetic_mutations_exp')
    pie_plot([numDiseaseMut_nonedgetic_exp, numDiseaseMut_edgetic_exp],
             angle = 90,
             #labels = ['Nonedgetic', 'Edgetic'],
             colors = ['thistle', 'mediumslateblue'],
             show = showFigs,
             figdir = figDir,
             figname = 'disease_edgetic_mutations_exp')
    
    #------------------------------------------------------------------------------------
    # Fraction of mono-edgetic mutations in experiments
    #------------------------------------------------------------------------------------
    
    if consider_exp_QuasiNull_perturbs:
        naturalPerturbNum_exp = expNaturalMutations.loc[ (expNaturalMutations["Edgotype_class"] == 'Edgetic') |
                                                         (expNaturalMutations["Edgotype_class"] == 'Quasi-null'),
                                                         "perturbations" ].apply(lambda x: sum( x > 0 ))
        diseasePerturbNum_exp = expDiseaseMutations.loc[ (expDiseaseMutations["Edgotype_class"] == 'Edgetic') |
                                                         (expNaturalMutations["Edgotype_class"] == 'Quasi-null'),
                                                         "perturbations" ].apply(lambda x: sum( x > 0 ))
    else:
        naturalPerturbNum_exp = expNaturalMutations.loc[ expNaturalMutations["Edgotype_class"] == 'Edgetic',
                                                         "perturbations" ].apply(lambda x: sum( x > 0 ))
        diseasePerturbNum_exp = expDiseaseMutations.loc[ expDiseaseMutations["Edgotype_class"] == 'Edgetic',
                                                         "perturbations" ].apply(lambda x: sum( x > 0 ))
    
    multi_histogram_plot ([naturalPerturbNum_exp,
                           diseasePerturbNum_exp],
                          ['g', 'r'],
                          xlabel = 'Number of PPIs perturbed',
                          ylabel = 'Frequency among edgetic mutations',
                          edgecolor = 'k',
                          fontsize = 14,
                          bins = 20,
                          alpha = 0.1,
                          leg = ['Non-disease', 'Disease'],
                          show = showFigs,
                          figdir = figDir,
                          figname = 'experiment_numPerturbs')
    
    # Number of edgetic (E) perturbations in response to mutations in experiment
    if consider_exp_QuasiNull_perturbs:
        numNatural_edgetic_exp = sum( expNaturalMutations["perturbations"].apply(lambda x: sum( x > 0 ) == 1 ) )
        numNatural_nonedgetic_exp = sum( expNaturalMutations["perturbations"].apply(lambda x: sum( x > 0 ) != 1) )
        numDisease_edgetic_exp = sum( expDiseaseMutations["perturbations"].apply(lambda x: sum( x > 0 ) == 1) )
        numDisease_nonedgetic_exp = sum( expDiseaseMutations["perturbations"].apply(lambda x: sum( x > 0 ) != 1) )
    else:
        numNatural_edgetic_exp = sum( (expNaturalMutations["Edgotype_class"] == 'Edgetic') &
                                      expNaturalMutations["perturbations"].apply(lambda x: sum(x > 0) == 1) )
        numNatural_nonedgetic_exp = sum( (expNaturalMutations["Edgotype_class"] != 'Edgetic') |
                                         expNaturalMutations["perturbations"].apply(lambda x: sum(x > 0) != 1) )
        numDisease_edgetic_exp = sum( (expDiseaseMutations["Edgotype_class"] == 'Edgetic') &
                                      expDiseaseMutations["perturbations"].apply(lambda x: sum(x > 0) == 1) )
        numDisease_nonedgetic_exp = sum( (expDiseaseMutations["Edgotype_class"] != 'Edgetic') |
                                         expDiseaseMutations["perturbations"].apply(lambda x: sum(x > 0) != 1) )
    
    print( '\n' + 'Fraction of experimental mono-edgetic mutations among edgetic mutations:' )
    print( 'non-disease mutation: %.3f (%d out of %d)' 
            % ( numNatural_edgetic_exp / numNaturalMut_edgetic_exp,
                numNatural_edgetic_exp,
                numNaturalMut_edgetic_exp ) )
    print( 'disease mutation: %.3f (%d out of %d)' 
            % ( numDisease_edgetic_exp / numDiseaseMut_edgetic_exp,
                numDisease_edgetic_exp,
                numDiseaseMut_edgetic_exp ) )
      
    # Total number of perturbations and non-perturbations in experiment
    numNatural_considered_exp = numNatural_edgetic_exp + numNatural_nonedgetic_exp
    numDisease_considered_exp = numDisease_edgetic_exp + numDisease_nonedgetic_exp
    
    # Fraction of edgetic perturbations in experiment
    fracNatural_edgetic_exp = numNatural_edgetic_exp / numNatural_considered_exp
    fracDisease_edgetic_exp = numDisease_edgetic_exp / numDisease_considered_exp
    
    # Standard error of the fraction of edgetic perturbations in experiment
    fracNatural_error_exp = sderror_on_fraction(numNatural_edgetic_exp, numNatural_considered_exp)
    fracDisease_error_exp = sderror_on_fraction(numDisease_edgetic_exp, numDisease_considered_exp)
    
    print( '\n' + 'Fraction of experimental mono-edgetic mutations:' )
    print('non-disease mutations: %.3f (SE = %g, %d out of %d)' % (fracNatural_edgetic_exp,
                                                                   fracNatural_error_exp,
                                                                   numNatural_edgetic_exp,
                                                                   numNatural_considered_exp))
    
    print('disease mutations: %.3f (SE = %g, %d out of %d)' % (fracDisease_edgetic_exp,
                                                               fracDisease_error_exp,
                                                               numDisease_edgetic_exp,
                                                               numDisease_considered_exp))
    
    fisher_test([numNatural_edgetic_exp, numNatural_nonedgetic_exp],
                [numDisease_edgetic_exp, numDisease_nonedgetic_exp])
    
    pie_plot([numNatural_nonedgetic_exp, numNatural_edgetic_exp],
             angle = 90,
             #labels = ['Nonedgetic', 'Edgetic'],
             colors = ['thistle', 'mediumslateblue'],
             show = showFigs,
             figdir = figDir,
             figname = 'non_disease_monoedgetic_mutations_exp')
    pie_plot([numDisease_nonedgetic_exp, numDisease_edgetic_exp],
             angle = 90,
             #labels = ['Nonedgetic', 'Edgetic'],
             colors = ['thistle', 'mediumslateblue'],
             show = showFigs,
             figdir = figDir,
             figname = 'disease_monoedgetic_mutations_exp')
        
    #----------------------------------------------------------------------------------------
    # Calculate overlap in proteins covered by predictions and experiments
    #----------------------------------------------------------------------------------------
    
    natural_proteins = set(naturalPerturbs["Protein"])
    disease_proteins = set(diseasePerturbs["Protein"])
    all_proteins = natural_proteins | disease_proteins
    outPath = interactomeOutDir / 'nondisease_mut_proteins.list'
    with open(outPath, 'w') as fout:
        for p in sorted(natural_proteins):
            fout.write(p + '\n')
    outPath = interactomeOutDir / 'disease_mut_proteins.list'
    with open(outPath, 'w') as fout:
        for p in sorted(disease_proteins):
            fout.write(p + '\n')
    
    print( '\n' + 'Our mutations:' )
    print('Non-disease: %d mutations in %d proteins' % (len(naturalPerturbs),
                                                        len(natural_proteins)))
    print('Disease: %d mutations in %d proteins' % (len(diseasePerturbs),
                                                    len(disease_proteins)))
    print('Total mutated proteins: %d' % len(natural_proteins | disease_proteins))
    
    natural_genes_exp = set(expNaturalMutations["Symbol"])
    disease_genes_exp = set(expDiseaseMutations["Symbol"])
    
    with open(UniprotIDmapFile, 'rb') as f:
        uniprotID = pickle.load(f)
    natural_proteins_exp = set([uniprotID[g] for g in natural_genes_exp if g in uniprotID])
    disease_proteins_exp = set([uniprotID[g] for g in disease_genes_exp if g in uniprotID])
    all_proteins_exp = natural_proteins_exp | disease_proteins_exp
    outPath = interactomeOutDir / 'nondisease_mut_proteins_exp.list'
    with open(outPath, 'w') as fout:
        for p in sorted(natural_proteins_exp):
            fout.write(p + '\n')
    outPath = interactomeOutDir / 'disease_mut_proteins_exp.list'
    with open(outPath, 'w') as fout:
        for p in sorted(disease_proteins_exp):
            fout.write(p + '\n')
    
    print( '\n' + 'Experimental mutations:' )
    print('Non-disease: %d mutations in %d genes, %d proteins' % (numNaturalMut_considered_exp,
                                                                  len(natural_genes_exp),
                                                                  len(natural_proteins_exp)))
    print('Disease: %d mutations in %d genes, %d proteins' % (numDiseaseMut_considered_exp,
                                                              len(disease_genes_exp),
                                                              len(disease_proteins_exp)))
    print('Total mutated proteins: %d' % len(natural_proteins_exp | disease_proteins_exp))
    
    natural_overlap_proteins = natural_proteins & natural_proteins_exp
    disease_overlap_proteins = disease_proteins & disease_proteins_exp
    all_proteins_overlap = all_proteins & all_proteins_exp
    
    print( '\n' + 'Overlap in proteins covered by all mutations:' )
    print( 'prediction proteins: %.1f %% (%d out of %d)' % ( 100 * len(all_proteins_overlap) 
                                                                 / len(all_proteins),
                                                             len(all_proteins_overlap),
                                                             len(all_proteins) ) )
    print( 'experiment proteins: %.1f %% (%d out of %d)' % ( 100 * len(all_proteins_overlap) 
                                                                 / len(all_proteins_exp),
                                                             len(all_proteins_overlap),
                                                             len(all_proteins_exp) ) )
    
    print( '\n' + 'Overlap in proteins covered by non-disease mutations:' )
    print( 'prediction proteins: %.1f %% (%d out of %d)' % ( 100 * len(natural_overlap_proteins) 
                                                                 / len(natural_proteins),
                                                             len(natural_overlap_proteins),
                                                             len(natural_proteins) ) )
    print( 'experiment proteins: %.1f %% (%d out of %d)' % ( 100 * len(natural_overlap_proteins) 
                                                                 / len(natural_proteins_exp),
                                                             len(natural_overlap_proteins),
                                                             len(natural_proteins_exp) ) )
    
    print( '\n' + 'Overlap in proteins covered by disease mutations:' )
    print( 'prediction proteins: %.1f %% (%d out of %d)' % ( 100 * len(disease_overlap_proteins) 
                                                                 / len(disease_proteins),
                                                             len(disease_overlap_proteins),
                                                             len(disease_proteins) ) )
    print( 'experiment proteins: %.1f %% (%d out of %d)' % ( 100 * len(disease_overlap_proteins) 
                                                                 / len(disease_proteins_exp),
                                                             len(disease_overlap_proteins),
                                                             len(disease_proteins_exp) ) )
    
    venn2_plot([natural_proteins, natural_proteins_exp],
               #labels = ['Prediction', 'Experiment'],
               colors = [interactome_color, 'orangered'],
               hideNum = True,
               show = showFigs,
               figdir = figDir,
               figname = 'non_disease_protein_overlap')
    venn2_plot([disease_proteins, disease_proteins_exp],
               #labels = ['Prediction', 'Experiment'],
               colors = [interactome_color, 'orangered'],
               hideNum = True,
               show = showFigs,
               figdir = figDir,
               figname = 'disease_protein_overlap')
    # plot legend circles
    venn2_plot([{1}, {2}],
               #labels = ['Prediction', 'Experiment'],
               colors = [interactome_color, 'orangered'],
               hideNum = True,
               show = showFigs,
               figdir = figDir,
               figname = 'legend_overlap')
    
    #------------------------------------------------------------------------------------
    # Compare gene ontology (GO) association of proteins in each dataset
    #------------------------------------------------------------------------------------
    
    # minimum number of tissue expression values required for protein pair tissue
    # co-expression to be considered
    coexprMinTissues = 5
    
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
    
    # Compute similarity in GO association between proteins used in predictions and 
    # proteins used in experiment, carrying non-disease mutations
    naturalProteins = list(natural_proteins)
    naturalProteinsExp = list(natural_proteins_exp)
    MFsim_natural = []
    BPsim_natural = []
    CCsim_natural = []
    for p1 in naturalProteins:
        for p2 in naturalProteinsExp:
            MFsim_natural.append(go_sim (p1, p2, MFassoc))
            BPsim_natural.append(go_sim (p1, p2, BPassoc))
            CCsim_natural.append(go_sim (p1, p2, CCassoc))
    MFsim_natural = [x for x in MFsim_natural if not np.isnan(x)]
    BPsim_natural = [x for x in BPsim_natural if not np.isnan(x)]
    CCsim_natural = [x for x in CCsim_natural if not np.isnan(x)]
    
    # Compute similarity in GO association between proteins used in predictions and 
    # proteins used in experiment, carrying disease mutations
    diseaseProteins = list(disease_proteins)
    diseaseProteinsExp = list(disease_proteins_exp)
    MFsim_disease = []
    BPsim_disease = []
    CCsim_disease = []
    for p1 in diseaseProteins:
        for p2 in diseaseProteinsExp:
            MFsim_disease.append(go_sim (p1, p2, MFassoc))
            BPsim_disease.append(go_sim (p1, p2, BPassoc))
            CCsim_disease.append(go_sim (p1, p2, CCassoc))
    MFsim_disease = [x for x in MFsim_disease if not np.isnan(x)]
    BPsim_disease = [x for x in BPsim_disease if not np.isnan(x)]
    CCsim_disease = [x for x in CCsim_disease if not np.isnan(x)]
    
    multi_bar_plot([[np.mean(MFsim_natural), np.mean(BPsim_natural), np.mean(CCsim_natural)],
                    [np.mean(MFsim_disease), np.mean(BPsim_disease), np.mean(CCsim_disease)]],
                   [[sderror(MFsim_natural), sderror(BPsim_natural), sderror(CCsim_natural)],
                    [sderror(MFsim_disease), sderror(BPsim_disease), sderror(CCsim_disease)]],
                   xlabels = ('Molecular function', 'Biological process', 'Cellular component'),
                   ylabel = ('Similarity (Jaccard index)'),
                   colors = ['turquoise', 'orangered'],
                   leg = ['Non-disease proteins', 'Disease proteins'],
                   fontsize = 24,
                   ylim = [0, 1],
                   show = showFigs,
                   figdir = figDir,
                   figname = 'protein_GO_overlap')
    
    # Compute tissue co-expression between proteins used in predictions and 
    # proteins used in experiment, carrying non-disease mutations
    coexpr_natural = []
    for p1 in naturalProteins:
        for p2 in naturalProteinsExp:
            coexpr_natural.append(coexpr (p1, p2, expr, minTissues=5))
    coexpr_natural = [x for x in coexpr_natural if not np.isnan(x)]
    
    # Compute tissue co-expression between proteins used in predictions and 
    # proteins used in experiment, carrying disease mutations
    coexpr_disease = []
    for p1 in diseaseProteins:
        for p2 in diseaseProteinsExp:
            coexpr_disease.append(coexpr (p1, p2, expr, minTissues=5))
    coexpr_disease = [x for x in coexpr_disease if not np.isnan(x)]
    
    bar_plot([np.mean(coexpr_natural),
              np.mean(coexpr_disease)],
             error = [sderror(coexpr_natural),
                      sderror(coexpr_disease)],
             xlabels = ('Non-disease proteins',
                        'Disease proteins'),
             ylabel = ('Co-expression'),
             colors = ['turquoise', 'orangered'],
             fontsize = 24,
             ylim = [0, 1],
             show = showFigs,
             figdir = figDir,
             figname = 'protein_coexpr_overlap')
    
    #------------------------------------------------------------------------------------
    # apply Bayes' theorem to calculate the fraction of PPIs that are junk, i.e., neutral 
    # under perturbation
    #------------------------------------------------------------------------------------
    
    junkPPIFile = interactomeOutDir / 'fraction_junk_PPIs.pkl'
    
    # Probability for new missense mutations to be neutral (N)
    pN = 0.27

    # Probability for new missense mutations to be mildly deleterious (M)
    pM = 0.53

    # Probability for new missense mutations to be strongly detrimental (S)
    pS = 0.20

    # Probability for strongly detrimental mutations (S) to be edgetic (E)
    pE_S = 0

    pN_E_results = []
    pN_E_bounds = []
    
    for approach in ['geometrybased_prediction_perMut',
                     'physicsbased_prediction_perMut',
                     'experiment_perMut',
                     'geometrybased_prediction_perPPI',
                     'physicsbased_prediction_perPPI',
                     'experiment_perPPI']:
        if approach == 'geometrybased_prediction_perMut':
            pE_N = fracNaturalMut_edgetic
            pE_M = fracDiseaseMut_edgetic
            numEdgetic = [ numNaturalMut_edgetic, numDiseaseMut_edgetic ]
            numMut = [ numNaturalMut_considered, numDiseaseMut_considered ]
            print( '\n' + 'Using geometry-based prediction of edgetic mutations' )
        elif approach == 'physicsbased_prediction_perMut':
            pE_N = fracNaturalMut_edgetic_ddg
            pE_M = fracDiseaseMut_edgetic_ddg
            numNaturalMut_edgetic_ddg = round( fracNaturalMut_edgetic_ddg * numNaturalMut_considered )
            numDiseaseMut_edgetic_ddg = round( fracDiseaseMut_edgetic_ddg * numDiseaseMut_considered )
            numEdgetic = [ numNaturalMut_edgetic_ddg, numDiseaseMut_edgetic_ddg ]
            numMut = [numNaturalMut_considered, numDiseaseMut_considered]
            print( '\n' + 'Using physics-based predicted edgetic mutations' )
        elif approach == 'experiment_perMut':
            pE_N = fracNaturalMut_edgetic_exp
            pE_M = fracDiseaseMut_edgetic_exp
            numEdgetic = [ numNaturalMut_edgetic_exp, numDiseaseMut_edgetic_exp ]
            numMut = [ numNaturalMut_considered_exp, numDiseaseMut_considered_exp ]
            print( '\n' + 'Using edgetic mutations from experiments' )
        elif approach == 'geometrybased_prediction_perPPI':
            pE_N = fracNatural_edgetic
            pE_M = fracDisease_edgetic
            numEdgetic = [numNatural_edgetic, numDisease_edgetic]
            numMut = [numNatural_considered, numDisease_considered]
            print( '\n' + 'Using geometry-based prediction of mono-edgetic mutations' )
        elif approach == 'physicsbased_prediction_perPPI':
            pE_N = fracNatural_edgetic_ddg
            pE_M = fracDisease_edgetic_ddg
            numNatural_edgetic_ddg = round( fracNatural_edgetic_ddg * numNatural_considered )
            numDisease_edgetic_ddg = round( fracDisease_edgetic_ddg * numDisease_considered )
            numEdgetic = [ numNatural_edgetic_ddg, numDisease_edgetic_ddg ]
            numMut = [ numNatural_considered, numDisease_considered]
            print( '\n' + 'Using physics-based prediction of mono-edgetic mutations' )
        elif approach == 'experiment_perPPI':
            pE_N = fracNatural_edgetic_exp
            pE_M = fracDisease_edgetic_exp
            numEdgetic = [numNatural_edgetic_exp, numDisease_edgetic_exp]
            numMut = [numNatural_considered_exp, numDisease_considered_exp]
            print( '\n' + 'Using mono-edgetic mutations from experiments' )
        
        # Probability for a new missense mutation to be edgetic
        pE = (pE_N * pN) + (pE_M * pM) + (pE_S * pS)
        
        # Probability for edgetic mutations to be effectively neutral
        pN_E = pE_N * pN / pE
        
        pN_E_results.append( 100 * pN_E )
        
        print( 'P(N) = %.1f %%' % (100 * pN) )
        print( 'P(M) = %.1f %%' % (100 * pM) )
        print( 'P(S) = %.1f %%' % (100 * pS) )
        print( 'P(E|N) = %.1f %%' % (100 * pE_N) )
        print( 'P(E|M) = %.1f %%' % (100 * pE_M) )
        print( 'P(E|S) = %.1f %%' % (100 * pE_S) )
        print( 'P(E) = P(E|N)*P(N) + P(E|M)*P(M) + P(E|S)*P(S) = %.1f %%' % (100 * pE) )
        print( 'Fraction of junk PPIs P(N|E) = P(E|N)*P(N)/P(E) = %.1f %%' % (100 * pN_E) )
        
        # calculate 95% confidence interval
        if computeConfidenceIntervals:
            n_N, n_M = numMut
            k_obs_N, k_obs_M = numEdgetic
            pE_M_pE_N_lower, pE_M_pE_N_upper = proportion_ratio_CI (k_obs_M, n_M, k_obs_N, n_N)
            if ( np.isnan( pE_M_pE_N_lower ) or np.isnan( pE_M_pE_N_upper ) ):
                print( '%.1f%% confidence interval for P(N|E) = NA' % CI )
                pN_E_bounds.append( [ 0 , 0 ] )
            else:
                pN_E_lower = 1 / ( pE_M_pE_N_upper * (pM / pN) + 1 )
                pN_E_upper = 1 / ( pE_M_pE_N_lower * (pM / pN) + 1 )
                print( '%.1f%% confidence interval for P(N|E) = (%f, %f)' % (CI,
                                                                             100 * pN_E_lower,
                                                                             100 * pN_E_upper) )
                pN_E_bounds.append( [100 * (pN_E - pN_E_lower), 100 * (pN_E_upper - pN_E)] )
    with open(junkPPIFile, 'wb') as fOut:
        pickle.dump([pN_E_results, pN_E_bounds], fOut)
    
    # plot fraction of junk PPIs from predictions and experiments
    if computeConfidenceIntervals:
        upper = [ p + upper for p, (lower, upper) in zip(pN_E_results, pN_E_bounds) ]
        maxY = max( upper[ : 3 ] )
    else:
        maxY = max( pN_E_results[ : 3 ] )
    maxY = 5 * np.ceil( maxY / 5 )
    bar_plot(pN_E_results[ : 3 ],
             pN_E_bounds[ : 3 ] if computeConfidenceIntervals else [ 0 ],
             xlabels = ('Geometry-based\nprediction',
                        'Physics-based\nprediction',
                        'Experiment'),
             ylabels = np.arange(0, maxY + 5, 5),
             ylabel = ('Fraction of junk PPIs (%)'),
             colors = 'white',
             fmt = 'k.',
             capsize = 10 if computeConfidenceIntervals else 0,
             msize = 11,
             ewidth = 1,
             #ecolors = [interactome_color, interactome_color, 'orangered'],
             ecolors = 'black',
             fontsize = 14,
             xlim = [0.8, 3.1],
             ylim = [0, maxY],
             yMinorTicks = True,
             adjustBottom = 0.2,
             shiftBottomAxis = -0.1,
             xbounds = (1,3),
             show = showFigs,
             figdir = figDir,
             figname = 'Fraction_junk_PPIs_perMut')
    if computeConfidenceIntervals:
        upper = [ p + upper for p, (lower, upper) in zip(pN_E_results, pN_E_bounds) ]
        maxY = max( upper[ 3 : ] )
    else:
        maxY = max( pN_E_results[ 3 : ] )
    maxY = 5 * np.ceil( maxY / 5 )
    bar_plot(pN_E_results[ 3 : ],
             pN_E_bounds[ 3 : ] if computeConfidenceIntervals else [ 0 ],
             xlabels = ('Geometry-based\nprediction',
                        'Physics-based\nprediction',
                        'Experiment'),
             ylabels = np.arange(0, maxY + 5, 5),
             ylabel = ('Fraction of junk PPIs (%)'),
             colors = 'white',
             fmt = 'k.',
             capsize = 10 if computeConfidenceIntervals else 0,
             msize = 11,
             ewidth = 1,
             #ecolors = [interactome_color, interactome_color, 'orangered'],
             ecolors = 'black',
             fontsize = 14,
             xlim = [0.8, 3.1],
             ylim = [0, maxY],
             yMinorTicks = True,
             adjustBottom = 0.2,
             shiftBottomAxis = -0.1,
             xbounds = (1,3),
             show = showFigs,
             figdir = figDir,
             figname = 'Fraction_junk_PPIs_perPPI')

def create_perturbed_network (interactome,
                              interactome_perturbations,
                              network_outPath,
                              nodeColor_outPath,
                              other_perturbations = None):
    """Create network from perturbed interactome for plotting by Cytoscape software.

    Args:
        interactome (DataFrame): interactome with all edges, inlcuding perturbed.
        interactome_perturbations (DataFrame): interactome perturbations.
        network_outPath (str): file directory to save network edges with their colors.
        nodeColor_outPath (str): file directory to save node colors.
    
    """
    network = interactome[ ["Protein_1", "Protein_2"] ].copy()
    network["Edge_color"] = 'black'
#     edges = []
#     for _, row in network.iterrows():
#         edges.append( (row.Protein_1, row.Protein_2) )
#     edgeColors = ['black'] * len(edges)
    nodes = list( set( network[ ["Protein_1", "Protein_2"] ].values.flatten() ) )
    numMut = [ 0 ] * len(nodes)
    for _, row in interactome_perturbations.iterrows():
        ind = nodes.index(row.Protein)
        numMut[ ind ] += 1
        for partner, perturbation in zip(row.partners, row.perturbations):
            if perturbation > 0:
#                 try:
#                     ind = edges.index( (row.Protein, partner) )
#                 except ValueError:
#                     ind = edges.index( (partner, row.Protein) )
#                 edgeColors[ind] = 'red'
#                 network.loc[network.apply(lambda x: 
#                                           {x["Protein_1"], x["Protein_2"]}, axis=1) 
#                             == set( edges[ind] ),
#                             "Edge_color"] = 'red'
                ind = ( ( (network["Protein_1"] == row.Protein) & (network["Protein_2"] == partner) ) |
                        ( (network["Protein_1"] == partner) & (network["Protein_2"] == row.Protein) ) )
                network.loc[ind, "Edge_color"] = 'red'
    
    if other_perturbations is not None:
        for _, row in other_perturbations.iterrows():
#         ind = nodes.index(row.Protein)
#         numMut[ind] += 1
            for partner, perturbation in zip(row.partners, row.perturbations):
                if perturbation > 0:
                    ind = ( ( (network["Protein_1"] == row.Protein) & (network["Protein_2"] == partner) ) |
                            ( (network["Protein_1"] == partner) & (network["Protein_2"] == row.Protein) ) )
                    perturbed = ( network["Edge_color"] == 'red' )
                    network.loc[ind & perturbed, "Edge_color"] = 'blue'
    
    network["Edge"] = network.apply(lambda x: x["Protein_1"] + ' (pp) ' + x["Protein_2"], axis=1)
    network.to_csv(network_outPath, index=False, sep='\t')

    # assign colors to network nodes
    nodeColors = []
    for num in numMut:
        if num > 9:
            nodeColors.append('blue')
        elif num > 2:
            nodeColors.append('mediumpurple')
        elif num > 0:
            nodeColors.append('magenta')
        else:
            nodeColors.append('white')
    
    # write network node colors to file
    with io.open(nodeColor_outPath, "w") as fout:
        fout.write('Protein' + '\t' + 'Color' + '\n')
        for node, color in zip(nodes, nodeColors):
            fout.write(node + '\t' + color + '\n')
    
#     return nodes, edges, nodeColors, edgeColors
    return ( nodes,
             network.apply(lambda x: ( x["Protein_1"], x["Protein_2"] ), axis=1).tolist(),
             nodeColors,
             network["Edge_color"].tolist() )

if __name__ == "__main__":
    main()
