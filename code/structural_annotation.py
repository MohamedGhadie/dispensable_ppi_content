import os
import io
import time
import pickle
import pandas as pd
import numpy as np
from pathlib import Path
from random import sample
from Bio.PDB import PDBParser
from simple_tools import create_dir, isolate_pairs, merge_list_pairs
from text_tools import read_list_table
from interactome_tools import (write_chain_annotated_interactome,
                               read_chain_annotated_interactome,
                               write_unmerged_interface_annotated_interactome,
                               read_unmerged_interface_annotated_interactome,
                               write_single_interface_annotated_interactome,
                               read_interface_annotated_interactome,
                               write_interface_annotated_interactome)
from pdb_tools import (allow_pdb_downloads,
                       suppress_pdb_warnings,
                       clear_structures,
                       return_structure,
                       load_pdbtools_chain_strucRes_labels,
                       load_pdbtools_chain_sequences,
                       return_chain_sequence,
                       ordered_chain_residues,
                       return_chain_res_posToID,
                       get_interface_by_chainIDs)

known_interfaces = {}

def load_dictionaries (chainSequenceFile = None, chainStrucResLabelFile = None):
    
    if chainSequenceFile:
        load_pdbtools_chain_sequences (chainSequenceFile)
    if chainStrucResLabelFile:
        load_pdbtools_chain_strucRes_labels (chainStrucResLabelFile)

def clear_interfaces ():
    """Clear dict of known chain interfaces
    
    """
    global known_interfaces
    print('\tclearing known interfaces dictionary')
    known_interfaces.clear()

def produce_alignment_evalue_dict (inPath, outPath, method = 'min'):
    
    chainMap = pd.read_table(inPath, sep='\t')
    evalueDict = {}
    if method is 'min':
        for i, row in chainMap.iterrows():
            k = row.Query + '-' + row.Subject
            if k in evalueDict:
                if (row.Expect < evalueDict[k]):
                    evalueDict[k] = row.Expect
            else:
                evalueDict[k] = row.Expect
    elif method is 'mean':
        for i, row in chainMap.iterrows():
            k = row.Query + '-' + row.Subject
            if k in evalueDict:
                evalueDict[k].append( row.Expect )
            else:
                evalueDict[k] = [ row.Expect ]
        for k in evalueDict:
            evalueDict[k] = np.mean( evalueDict[k] )
    with open(outPath, 'wb') as fOut:
        pickle.dump(evalueDict, fOut)
    
def locate_alignments (inPath,
                       outPath,
                       resMatch = False,
                       pausetime = 0):
    """Transform alignement matching sequences to query and subject
        sequence residue positions.

    Args:
        inPath (str): directory containing alignment table.
        pausetime (int): pausing time between locating alignments on queries and subjects.
        outPath (str): directory to save alignment residue positions to.

    """
    alignments = pd.read_table(inPath, sep="\t")
    print('\t%d alignments to locate' % len(alignments))
    
    print('\tlocating alignments on query sequences')
    alignments["Qpos"] = alignments.apply(lambda x: locate_alignment(x["Qseq"],
                                                                     x["Sseq"],
                                                                     x["Qstart"],
                                                                     resMatch = resMatch), axis=1)
    
    print('\tpausing for %d seconds before locating alignments on PDB sequences' % pausetime)
    time.sleep(pausetime)
    
    print('\tlocating alignments on PDB sequences')
    alignments["Spos"] = alignments.apply(lambda x: locate_alignment(x["Sseq"],
                                                                     x["Qseq"],
                                                                     x["Sstart"],
                                                                     resMatch = resMatch), axis=1)
    
    alignments.drop(["Qseq","Sseq","Match"], axis=1, inplace=True)
    alignments["Qpos"] = alignments["Qpos"].apply(lambda x: ','.join(map(str, x)))
    alignments["Spos"] = alignments["Spos"].apply(lambda x: ','.join(map(str, x)))
    alignments.to_csv(outPath, index=False, sep='\t')

def locate_alignment (Qseq,
                      Sseq,
                      Qstart,
                      ResMatch = False):
    """Transform a single alignement matching sequence to query sequence residue positions.

    Args:
        Qseq (str): query sequence.
        Sseq (str): subject sequence.
        Qstart (int): alignment start position on query sequence.
        ResMatch (boolean): True if aligned positions must have same residue, otherwise False.
    
    Returns:
        List: matching residue positions in query sequence.

    """
    if ResMatch:
        matchPos = [i for i, ch in enumerate(Qseq) if (ch != '-') and (Sseq[i] == ch)]
    else:
        matchPos = [i for i, ch in enumerate(Qseq) if (ch != '-') and (Sseq[i] != '-')]
    gapPos = [i for i, ch in enumerate(Qseq) if ch == '-']
    if len(gapPos) == 0:
        return [pos + Qstart for pos in matchPos]
    else:
        numGaps = [len([g for g in gapPos if g < pos]) for pos in matchPos]
        return [pos + Qstart - gaps for pos, gaps in zip(matchPos, numGaps)]
                                 
def map_positions (refPos, alignPos, pos):
    """Map a list of reference positions to new positions given an allignment.

    Args:
        refPos (list): aligned reference positions.
        alignPos (list): aligned new positions.
        pos (list): reference positions to map.
    
    Returns:
        tuple: mapped positions and fraction of those mapped through allignment.

    """
    mappedPos = [alignPos[refPos.index(i)] for i in pos if i in refPos]
    frac = round( len(mappedPos) / len(pos), 2)
    return mappedPos, frac

def produce_chain_annotated_interactome (inPath,
                                         proteinChainsFile,
                                         outPath,
                                         alignmentEvalueFile = None):
    """Annotate interacting proteins with pairs of chains from the same structure

    Args:
        inPath(): file directory containing table of protein-protein interactions
        chainMapFile (str): file directory containing dict of chains mapping to each protein
        outPath (str): file directory to save chain-annotated interactome.

    """
    annotatedInteractome = pd.read_table(inPath, sep="\t")
    with open(proteinChainsFile, 'rb') as f:
        proteinChains = pickle.load(f)
    if alignmentEvalueFile:
        with open(alignmentEvalueFile, 'rb') as f:
            alignmentEvalues = pickle.load(f)
    else:
        alignmentEvalues = None
    
    annotatedInteractome["Mapping_chains"] = annotatedInteractome.apply(lambda x:
                                                    sameStruct_mapping_chains(x["Protein_1"],
                                                                              x["Protein_2"],
                                                                              proteinChains,
                                                                              alignmentEvalues = alignmentEvalues),
                                                    axis=1)
    
    annotatedInteractome = annotatedInteractome[annotatedInteractome["Mapping_chains"].apply(len)
                                                > 0].reset_index(drop=True)
    print('\t%d interactions annotated with same-structure chains' % len(annotatedInteractome))
    write_chain_annotated_interactome (annotatedInteractome, outPath)
    
def sameStruct_mapping_chains (protein1,
                               protein2,
                               proteinChains,
                               alignmentEvalues = None):
    """Identify same-structure chain pairs mapping to two interacting proteins

    Args:
        protein1 (str): ID of first interacting protein
        protein2 (str): ID of second interacting protein
        proteinChains (dict): dict of chains mapping to each protein

    Returns:
        List: list of same-structure chain pairs mapping to the protein-protein interaction

    """
    chainPairs = []
    if (protein1 in proteinChains) and (protein2 in proteinChains):
        p1chains = sorted(set(proteinChains[protein1]))
        p2chains = sorted(set(proteinChains[protein2]))
    
        pdbID1 = [c.split('_')[0] for c in p1chains]
        pdbID2 = [c.split('_')[0] for c in p2chains]
        
        for i, c1 in enumerate(p1chains):
            for j, c2 in enumerate(p2chains):
                if (c1 != c2) and (pdbID1[i] == pdbID2[j]):
                    chainPairs.append( (c1, c2) )
        
        if alignmentEvalues:
            avgEvalues = []
            for c1, c2 in chainPairs:
                evalue1 = alignmentEvalues[protein1 + '-' + c1]
                evalue2 = alignmentEvalues[protein2 + '-' + c2]
                avgEvalues.append( np.mean( [evalue1, evalue2] ) )
            sortedPairs = [ chainPairs[i] for i, _ in sorted(enumerate(avgEvalues), key=lambda x: x[1]) ]
            return sortedPairs   
    return chainPairs

def filter_chain_annotations(inPath,
                             evalue,
                             chainCoverage,
                             outPath):
    
    chainMap = pd.read_table(inPath, sep='\t')
    chainMap = chainMap[chainMap["Expect"] < evalue]
    chainMap = chainMap[(chainMap["Send"] - chainMap["Sstart"])
                        >= chainMap["Slen"] * chainCoverage].reset_index(drop=True)
    chainMap = chainMap.sort_values("Expect", axis=0, ascending=True)
    chainMap = chainMap.drop_duplicates(subset=['Query','Subject'], keep='first')
    chainMap = chainMap.sort_values(['Query','Subject'], axis=0, ascending=True)
    chainMap.to_csv(outPath, index=False, sep='\t')

def filter_chain_annotations_by_protein (inPath, proteins, outPath):
    
    with io.open(inPath, "r", errors='ignore') as f, io.open(outPath, "w") as fout:
        headers = f.readline()
        fout.write( headers )
        headerSplit = headers.strip().split('\t')
        numCol = len( headerSplit )
        queryPos = headerSplit.index('Query')
    with io.open(inPath, "r", errors='ignore') as f, io.open(outPath, "a") as fout:
        next(f)
        for line in f:
            linesplit = line.strip().split('\t')
            if len( linesplit ) == numCol:
                if linesplit[ queryPos ] in proteins:
                    fout.write( line )

def produce_interface_annotated_interactome (inPath,
                                             pdbDir,
                                             chainMapFile,
                                             interfaceFile,
                                             chainStrucResFile,
                                             maxInterfaces,
                                             maxAttempts,
                                             rnd,
                                             mapCutoff,
                                             bindingDist,
                                             outPath,
                                             downloadPDB = True,
                                             suppressWarnings = False):
    """Map interaction interfaces onto protein-protein interaction network.

    Args:
        inPath (str): file directory containing protein interaction network.
        pdbDir (str): directory containing pdb structure files.
        chainMapFile (str): file directory containing table of chains mapping onto each protein.
        interfaceFile (str): file directory containing known interfaces of each chain with its
                                binding chain. Newly calculated interfaces are added to this file.
        chainStrucResFile (str): file directory containing dict of labels for structured 
                                    residues that are part of SEQRES.
        numinf (int): number of chaninPairs to use for interface mapping.
        rnd (boolean): True to randomly select chain pairs from PPI annotation list 
                        for interface mapping.
        outPath (str): file directory to save interface-annotated interactome to.

    """
    global known_interfaces
    clear_structures()
    clear_interfaces()
    allow_pdb_downloads (downloadPDB)
    suppress_pdb_warnings (suppressWarnings)
    load_dictionaries (chainSequenceFile=chainSeqFile, chainStrucResLabelFile=chainStrucResFile)
    print('\treading chain-annotated interactome')
    interactome = read_chain_annotated_interactome(inPath)
    print('\treading protein-chain mapping table')
    chainMap = read_list_table(chainMapFile, ["Qpos", "Spos"], [int, int], '\t')
    annotatedPPIs, numPPIs = 0, len(interactome)
    
    if interfaceFile.is_file():
        print('\tloading known chain interfaces')
        read_chain_interfaces(interfaceFile)
    else:
        print('creating interface file directory')
        filedir, _ = os.path.split(interfaceFile)
        create_dir(filedir)
    
    print('\tresolving interfaces for all %d PPIs' % numPPIs)
    interfaces = []
    for i, ppi in interactome.iterrows():
        print('\t' + 'PPI: ' + protein1 + ' - ' + protein2)
        sel_interfaces, sel_frac, sel_chainpairs = map_multi_interfaces (pdbDir,
                                                                         ppi.Protein_1,
                                                                         ppi.Protein_2,
                                                                         ppi.Mapping_chains,
                                                                         chainMap,
                                                                         maxInterfaces,
                                                                         maxAttempts,
                                                                         rnd,
                                                                         mapCutoff,
                                                                         bindingDist,
                                                                         interfaceFile,
                                                                         PPIcount = i + 1,
                                                                         numPPIs = numPPIs)
        interface.append( (sel_interfaces, sel_frac, sel_chainpairs) )
        if sel_chainpairs:
            annotatedPPIs += 1
        print('\t' + '%d PPIs resolved' % (i + 1))
        print('\t' + '%d PPIs annotated' % annotatedPPIs)
    interactome["Interfaces"] = interfaces
    
    print('\tfinalizing interface-annotated interactome')
    interactome = interactome[interactome["Interfaces"].apply(lambda x: len(x[0]) > 0)].reset_index(drop=True)
    interactome.drop("Mapping_chains", axis=1, inplace=True)
    interactome["Chain_pairs"] = interactome["Interfaces"].apply(lambda x: x[2])
    interactome["Mapping_frac"] = interactome["Interfaces"].apply(lambda x: x[1])
    interactome["Interfaces"] = interactome["Interfaces"].apply(lambda x: x[0])
    write_unmerged_interface_annotated_interactome (interactome, outPath)

def map_multi_interfaces (pdbDir,
                          protein1,
                          protein2,
                          chainPairs,
                          chainMap,
                          maxInterfaces,
                          maxAttempts,
                          rnd,
                          mapCutoff,
                          bindingDist,
                          interfaceFile,
                          PPIcount = None,
                          numPPIs = None):
    """Map the first n interaction interfaces passing a mapping fraction cutoff for two 
        interacting proteins from mapping chain pairs from distinct PDB models.

    Args:
        pdbDir (str): directory containing pdb structure files.
        protein1 (str): ID of first interacting protein.
        protein2 (str): ID of second interacting protein.
        chainPairs (list): list of chain pairs mapping onto the two interacting proteins.
        chainMap (dataframe): table of chains mapping onto each protein.
        numinf (int): number of chaninPairs to use for interface mapping.
        rnd (boolean): True to randomly select pairs from chain pairs list for interface mapping.
        mapCutoff (float): minimum mapping fraction needed for interface from a chain pair
                            to be selected as the only interface.
        interfaceFile (str): file directory containing known interfaces of each chain with its
                                binding chain. Newly calculated interfaces are added to this file.
    
    Returns:
        list: all interfaces between the two proteins in tuples.
        list: mapping fraction for each two interfaces in tuples.

    """
    numAttempts, interfacesFound, numChainPairs = 0, 0, len(chainPairs)
    sel_interfaces, sel_frac, sel_chainpairs, sel_pdbs = [], [], [], []
    if rnd:
        chPairs = sample(chainPairs, numChainPairs)
    else:
        chPairs = chainPairs
    for i, chainPair in enumerate(chPairs):
        if (numAttempts == maxAttempts) or (interfacesFound == maxInterfaces):
            break
        chainID1, chainID2 = chainPair
        pdbID, _ = chainID1.split('_')
        if ppiCount and numPPIs:
            print('\t\t' + 'PPI: %s - %s (# %d out of %d) (%.1f %%)' 
                    % (protein1, protein2, PPIcount, numPPIs, 100 * PPIcount / numPPIs))
        else:
            print('\t\t' + 'PPI: %s - %s' % (protein1, protein2))
        print('\t\t' + '%d attempts (max %d)' % (numAttempts, maxAttempts))
        print('\t\t' + '%d interfaces mapped (seeking %d)' % (interfacesFound, maxInterfaces))
        print('\t\t' + 'Chain pair: %s - %s (# %d out of %d)' % (chainID1, chainID2, i + 1, numChainPairs))
        if pdbID not in sel_pdbs:
            interfaces1, interfaces2 = map_twoside_interfaces (pdbDir,
                                                               protein1,
                                                               protein2,
                                                               chainPair,
                                                               bindingDist,
                                                               chainMap,
                                                               interfaceFile)
            passFrac1 = [f for _, f in interfaces1 if f >= mapCutoff]
            passFrac2 = [f for _, f in interfaces2 if f >= mapCutoff]
            if (len(passFrac1) > 0) and (len(passFrac2) > 0):
                inf1 = [i for i, f in interfaces1 if f >= mapCutoff]
                inf2 = [i for i, f in interfaces2 if f >= mapCutoff]
                sel_interfaces.append((inf1, inf2))
                sel_frac.append((passFrac1, passFrac2))
                sel_chainpairs.append(([chainID1]*len(inf1), [chainID2]*len(inf2)))
                sel_pdbs.append(pdbID)
                interfacesFound += 1
                print('\t\t\t' + 'interface %d found' % interfacesFound)
            numAttempts += 1
        else:
            print('\t\t\t' + 'interface already mapped from PDB structure')
    return sel_interfaces, sel_frac, sel_chainpairs

def map_twoside_interfaces (pdbDir,
                            protein1,
                            protein2,
                            chainPair,
                            bindingDist,
                            chainMap,
                            interfaceFile):
    
    global known_interfaces
    chain1, chain2 = chainPair
    chainKey1 = chain1 + '-' + chain2
    chainKey2 = chain2 + '-' + chain1
    if (chainKey1 not in known_interfaces) or (chainKey2 not in known_interfaces):
        print('\t\t\tcomputing interface for chain pair ' + chainKey1)
        known_interfaces[chainKey1], known_interfaces[chainKey2] = get_interface_by_chainIDs (pdbDir,
                                                                                              chain1,
                                                                                              chain2,
                                                                                              maxDist = bindingDist,
                                                                                              suppressWarnings = True)
        print('\t\t\t' + 'writing chain interfaces to file')
        write_chain_interfaces(interfaceFile)
    else:
        print('\t\t\tinterface known for chain pair ' + chainKey1)
    if (len(known_interfaces[chainKey1]) > 0) and (len(known_interfaces[chainKey2]) > 0):
        alignment1 = chainMap[(chainMap["Query"]==protein1) &
                              (chainMap["Subject"]==chain1)]
        alignment2 = chainMap[(chainMap["Query"]==protein2) &
                              (chainMap["Subject"]==chain2)]
        interfaces1 = alignment1.apply(lambda x:
                                             map_positions(x["Spos"],
                                                           x["Qpos"],
                                                           known_interfaces[chainKey1]),
                                              axis=1)
        interfaces2 = alignment2.apply(lambda x:
                                             map_positions(x["Spos"],
                                                           x["Qpos"],
                                                           known_interfaces[chainKey2]),
                                              axis=1)
        interfaces1 = [(i, f) for i, f in interfaces1.values if f > 0]
        interfaces2 = [(i, f) for i, f in interfaces2.values if f > 0]
        return interfaces1, interfaces2
    else:
        return [], []

def merge_interactome_interface_annotations (inPath, outPath):
    """Merge interface annotations for each PPI in an interactome into one tuple.

    Args:
        inPath (str): file directory containing interface-annotated interactome.
        outPath (str): file directory to save merged-interface-annotated interactome to.
    
    """
    interactome = read_unmerged_interface_annotated_interactome ( inPath )
    interactome["Interfaces"] = interactome["Interfaces"].apply( isolate_pairs )
    interactome["Interfaces"] = interactome["Interfaces"].apply( merge_list_pairs )
    interactome["Chain_pairs"] = interactome["Chain_pairs"].apply( isolate_pairs )
    interactome = interactome.drop("Mapping_frac", axis=1)
    write_single_interface_annotated_interactome (interactome, outPath)

def remove_duplicate_interface_annotations (inPath, outPath):
    """Remove duplicate interfaces annotations from interactome.

    Args:
        inPath (str): file directory containing interface-annotated interactome.
        outPath (str): file directory to save interface-annotated interactome to.
    
    """
    interactome = read_unmerged_interface_annotated_interactome ( inPath )
    interactome["Interfaces"] = interactome["Interfaces"].apply( isolate_pairs )
    interactome["Mapping_frac"] = interactome["Mapping_frac"].apply( isolate_pairs )
    
    uniqueInterfaces = interactome.apply(lambda x:
                                           drop_duplicate_interfaces(x["Interfaces"],
                                                                     x["Mapping_frac"]),
                                           axis=1)
    interactome["Interfaces"] = uniqueInterfaces.apply(lambda x: x[0])
    interactome = interactome.drop("Mapping_frac", axis=1)
    write_interface_annotated_interactome( outPath )

def drop_duplicate_interfaces(interfaces, fractions):
    """Remove duplicate interfaces (tuples) from a list of interfaces.

    Args:
        interfaces (list): list of interfaces in tuple form.
        fractions (list): mapping fractions associated with the interfaces.
    
    Returns:
        list: unique interfaces.
        list: mapping fractions associated with interfaces.
    """
    uniqueInt = []
    uniqueFrac = []
    for i, intfc in enumerate(interfaces[:-1]):
        if all([intfc != intfc2 for intfc2 in interfaces[i+1:]]):
            uniqueInt.append(intfc)
            uniqueFrac.append(fractions[i])
    uniqueInt.append(interfaces[-1])
    uniqueFrac.append(fractions[-1])
    return [uniqueInt, uniqueFrac]

def filter_interfaces_by_frac (interfaces,
                               fractions,
                               cutoff):
    """Filter interfaces (tuples) from a list of interfaces by mapping fraction. Keep
        interfaces whose both side map with a fraction higher than a cutoff.

    Args:
        interfaces (list): list of interfaces in tuple form.
        fractions (list): mapping fractions associated with the interfaces.
        cutoff (float): cutoff on mapping fraction.
    
    Returns:
        list: unique interfaces.
        list: mapping fractions associated with interfaces.
    
    """
    selInt = [interfaces[i] for i, (a, b) in enumerate(fractions) if (a >= cutoff) and (b >= cutoff)]
    selFrac = [(a, b) for (a, b) in fractions if (a >= cutoff) and (b >= cutoff)]
    return [selInt, selFrac]

def write_chain_interfaces (outPath):
    
    with io.open(outPath, "w") as fout:
        fout.write('\t'.join(["Chain_1", "Chain_2", "Chain1_interface"]) + '\n')
        for chainKey in sorted(known_interfaces.keys()):
            chain1, chain2 = chainKey.split('-')
            if len(known_interfaces[chainKey]) > 0:
                interface = ','.join([str(elm) for elm in known_interfaces[chainKey]])
            else:
                interface = '-1'
            fout.write('\t'.join([chain1, chain2, interface]) + '\n')

def read_chain_interfaces (inPath):
    
    if inPath.is_file():
        known_interfaces_df = read_list_table(inPath, "Chain1_interface", int, '\t')
        for _, row in known_interfaces_df.iterrows():
            chainKey = row.Chain_1 + '-' + row.Chain_2
            if -1 in row.Chain1_interface:
                known_interfaces[chainKey] = []
            else:
                known_interfaces[chainKey] = row.Chain1_interface

def write_pdb_mapped_mutations_old (mutations,
                                interactomeFile,
                                chainMapFile,
                                chainSeqFile,
                                proteinSeqFile,
                                chainResOrderFile,
                                pdbDir,
                                outPath):
    """Map mutations onto PDB chains and write to file the mutation, host protein, partner
    protein, position on host protein, position on chain, PDB ID, chain_mutation 
    (WT res, chain ID, pos on chain, mut res), partner chain.

    Args:
        mutations (DataFrame): mutations to be mapped onto PDB chains and writen to file.
        interactomeFile (str): file directory containing interactome network.
        chainMapFile (str): file directory containing table of chains mapping onto each protein.
        chainSeqFile (str): file directory containing dictionary of chain sequences.
        proteinSeqFile (str): file directory containing dictionary of protein sequences.
        chainResOrderFile (str): file directory containing dict of SEQRES residue positions
                                having structural coordinates, in same order.
        pdbDir (str): file directory containing PDB structure files.
        outPath (str): file directory to save mapped mutations to. 

    """
    print('\t' + 'loading interface-annotated interactome')
    interactome = read_interface_annotated_interactome(interactomeFile)
    
    print('\t' + 'loading protein-chain alignment table')
    chainMap = read_list_table(chainMapFile, ['Qpos', 'Spos'], [int, int], '\t')
    
    print('\t' + 'loading chain sequences')
    with open(chainSeqFile, 'rb') as f:
        chainSeq = pickle.load(f)
    
    print('\t' + 'loading protein sequences file')
    with open(proteinSeqFile, 'rb') as f:
        proteinSeq = pickle.load(f)
    
    print('\t' + 'loading chain structured residue order file')
    with open(chainResOrderFile, 'rb') as f:
        chainResOrder = pickle.load(f)
    
    print('\t' + 'writing mutations')
    with io.open(outPath, "a") as fout:
        fout.write('\t'.join(['protein',
                              'partner',
                              'protein_pos',
                              'chain_pos',
                              'pdb_id',
                              'Chain_mutation',
                              'Partner_chain']) + '\n')
        
        for i, mut in mutations.iterrows():
            print('\t' + 'mutation # %d' % i)
            perturbing = False
            pos = mut.Mutation_Position
            perturbedPartners = [p for p, perturb in zip(mut.partners, mut.perturbations) if perturb > 0]
            ppis = interactome[ (interactome[ ["Protein_1", "Protein_2"] ] == mut.Protein).any(1) ]
            for j, ppi in ppis.iterrows():
                if (ppi.Protein_1 in perturbedPartners) or (ppi.Protein_2 in perturbedPartners):
                    perturbing = True
                    print('\t\t' + 'PPI # %d' % j)
                    mapped = False
                    chainpairs = ppi.Chain_pairs
                    if ppi.Protein_2 == mut.Protein:
                        chainpairs = [tuple(reversed(x)) for x in chainpairs]
                        partner = ppi.Protein_1
                    else:
                        partner = ppi.Protein_2
                    for ch1, ch2 in chainpairs:
                        print('\t\t\t' + 'chain pair: %s-%s ' % (ch1, ch2))
                        pdbid, ch1_id = ch1.split('_')
                        _, ch2_id = ch2.split('_')
                        structureFile = pdbDir / ('pdb' + pdbid + '.ent')
                        if structureFile.is_file():
                            struc = PDBParser(QUIET=True).get_structure(pdbid, str(structureFile))
                            model = struc[0]
                            resOrder = chainResOrder[ ch1 ]
                            chain1 = model[ ch1_id ]
                            residues = [res for res in chain1.get_residues()][ : len(resOrder) ]
                            mappings = chainMap[ (chainMap["Query"] == mut.Protein) 
                                                 & (chainMap["Subject"] == ch1) ]
                            for _, mapping in mappings.iterrows():
                                try:
                                    mapPos = mapping.Spos[ mapping.Qpos.index(pos) ]
                                    wtRes = proteinSeq[ mut.Protein ][ pos - 1 ]
                                    chainRes = chainSeq[ ch1 ][ mapPos - 1 ]
                                    if chainRes != mut.Mut_res:
                                        try:
                                            ind = resOrder.index( mapPos )
                                            res = residues[ ind ]
                                            _, resNum, _ = res.get_id()
                                            fout.write('\t'.join([mut.Protein,
                                                                  partner,
                                                                  str(pos),
                                                                  str(mapPos),
                                                                  pdbid,
                                                                  chainRes + ch1_id
                                                                           + str(resNum)
                                                                           + mut.Mut_res,
                                                                  ch2_id]) +  '\n')
                                            print('\t\t\t' + 'mutation mapping writen to file')
                                            mapped = True
                                        except ValueError:
                                            print('\t\t\t' + 'chain residue coordinates not known')
                                    else:
                                        print('\t\t\t' + 'chain residue same as mutation residue')
                                except ValueError:
                                    print('\t\t\t' + 'mutation position not part of protein-chain alignment')
                        else:
                            print('\t\t\t' + 'structure file for chain ' + ch1 + ' not found')
                    if not mapped:
                        print('\t\t\t' + 'mutation mapping not successfull')
                    fout.write('\n')
            if not perturbing:
                print('\t\t' + 'perturbed PPI not found for mutation')

def write_mutation_structure_maps (mutations,
                                   interactomeFile,
                                   chainMapFile,
                                   chainSeqFile,
                                   proteinSeqFile,
                                   chainStrucResFile,
                                   pdbDir,
                                   outPath,
                                   chainInterfaceFile = None,
                                   downloadPDB = True,
                                   suppressWarnings = False):
    """Map mutations onto PDB chains and write to file the mutation, host protein, partner
    protein, position on host protein, position on chain, PDB ID, chain_mutation 
    (WT res, chain ID, pos on chain, mut res), partner chain.

    Args:
        mutations (DataFrame): mutations to be mapped onto PDB chains and writen to file.
        interactomeFile (str): file directory containing interactome network.
        chainMapFile (str): file directory containing table of chains mapping onto each protein.
        chainSeqFile (str): file directory containing dictionary of chain sequences.
        proteinSeqFile (str): file directory containing dictionary of protein sequences.
        pdbDir (str): file directory containing PDB structure files.
        outPath (str): file directory to save mapped mutations to. 

    """
    global known_interfaces
    clear_structures()
    clear_interfaces()
    allow_pdb_downloads (downloadPDB)
    suppress_pdb_warnings (suppressWarnings)
    load_dictionaries (chainSequenceFile=chainSeqFile, chainStrucResLabelFile=chainStrucResFile)
    interactome = read_interface_annotated_interactome(interactomeFile)
    chainMap = read_list_table(chainMapFile, ['Qpos', 'Spos'], [int, int], '\t')
    
    check_interface = False
    if chainInterfaceFile:
        if chainInterfaceFile.is_file():
            print('\t' + 'loading chain interfaces')
            read_chain_interfaces(chainInterfaceFile)
            check_interface = True
        else:
            warnings.warn('Chain interface file not found. Mutations will be mapped onto structure without checking interface')
    
    print('\t' + 'writing mutations')
    with io.open(outPath, "a") as fout:
        # write header line to file
        fout.write('\t'.join(['protein',
                              'partner',
                              'protein_pos',
                              'pdb_id',
                              'chain_id',
                              'chain_pos',
                              'chain_mutation',
                              'partner_chain']) + '\n')
        # go through each mutation
        for i, mut in mutations.iterrows():
            print('\t' + 'mutation index: %d' % i)
            perturbing = False
            pos = mut.Mutation_Position
            perturbedPartners = [p for p, perturb in zip(mut.partners, mut.perturbations) if perturb > 0]
            ppis = interactome[ (interactome[ ["Protein_1", "Protein_2"] ] == mut.Protein).any(1) ]
            
            # go through each interaction (PPI) perturbed by the mutation
            for j, ppi in ppis.iterrows():
                if (ppi.Protein_1 in perturbedPartners) or (ppi.Protein_2 in perturbedPartners):
                    perturbing, mapped = True, False
                    print('\t\t' + 'PPI # %d' % j)
                    
                    # get chain pairs (models) used to map interface for this PPI
                    if ppi.Protein_2 == mut.Protein:
                        chainPairs = [tuple(reversed(x)) for x in ppi.Chain_pairs]
                        partner = ppi.Protein_1
                    else:
                        chainPairs = ppi.Chain_pairs
                        partner = ppi.Protein_2
                    
                    # go through each pair of chains (model)
                    for ch1, ch2 in chainPairs:
                        print('\t\t\t' + 'chain pair: %s-%s ' % (ch1, ch2))
                        (pdbid, ch1_id), (_, ch2_id) = ch1.split('_'), ch2.split('_')
                        struc = return_structure (pdbid, pdbDir)
                        if struc:
                            model = struc[0]
                            # get structured residues that are part of the chain SEQRES
                            residues = ordered_chain_residues (pdbid, model, ch1_id, pdbDir)
                            if residues:
                                # map mutation position back onto chain pair (model) through sequence alignment
                                mappings = mutation_structure_map (chainMap, mut.Protein, ch1, pos)
                                # if mutation position maps through any protein-chain alignment 
                                if mappings is not None:
                                    # make sure chain wildtype residue is different than mutation residue
                                    mappings["chainRes"] = mappings["posMaps"].apply(lambda x: return_chain_sequence(ch1)[x-1])
                                    mappings = mappings[mappings["chainRes"] != mut.Mut_res]
                                    if check_interface:
                                        k = ch1 + '-' + ch2
                                        interfacial = mappings["posMaps"].apply(lambda x: x in known_interfaces[k]
                                                                                          if k in known_interfaces
                                                                                          else False)
                                        mappings = mappings[interfacial]
                                    if not mappings.empty:
                                        # select the first position map from alignment table
                                        chainRes, mapPos = mappings[["chainRes","posMaps"]].iloc[0]
                                        # map chain residue position to residue ID
                                        resID = return_chain_res_posToID (pdbid, ch1_id, mapPos, pdbDir)
                                        if resID:
                                            # write mutation structure mapping to file
                                            _, resNum, _ = resID
                                            fout.write('\t'.join([mut.Protein,
                                                                  partner,
                                                                  str(pos),
                                                                  pdbid,
                                                                  ch1_id,
                                                                  str(mapPos),
                                                                  ''.join([chainRes,
                                                                           ch1_id,
                                                                           str(resNum),
                                                                           mut.Mut_res]),
                                                                  ch2_id]) +  '\n')
                                            print('\t\t\t' + 'mutation mapping writen to file')
                                            mapped = True
                                        else:
                                            print('\t\t\t' + 'chain residue coordinates not known')
                                    else:
                                        print('\t\t\t' + 'chain residue same as mutation residue')
                                else:
                                    print('\t\t\t' + 'mutation position not part of protein-chain alignment')
                            else:
                                print('\t\t\t' + 'chain residues not found')
                        else:
                            print('\t\t\t' + 'no structure found for chain ' + ch1)
                    if not mapped:
                        print('\t\t\t' + 'mutation mapping not successfull')
                    fout.write('\n')
            if not perturbing:
                print('\t\t' + 'perturbed PPI not found for mutation')

def mutation_structure_map (strucMap,
                            protein,
                            chainID,
                            pos,
                            firstMap = False):
    
    mappings = strucMap [(strucMap["Query"] == protein) & 
                         (strucMap["Subject"] == chainID)].reset_index(drop=True)
    mappings["posMaps"] = position_map (pos, mappings)
    mappings = mappings [np.isnan(mappings["posMaps"]) == False]
    if firstMap and not mappings.empty:
        return mappings.iloc[0]
    else:
        return mappings if not mappings.empty else None

def position_map (resPos, strucMap):
    
    posMap = []
    for Qpos, Spos in strucMap[["Qpos", "Spos"]].values:
        if resPos in Qpos:
            posMap.append( Spos[Qpos.index(resPos)] )
        else:
            posMap.append(np.nan)
    return np.array(posMap)

# def read_mutation_ddg (inPath):
#     
#     """Read mutation change (loss) in binding free energy from file writen by 
#     function "write_pdb_mapped_mutations".
# 
#     Args:
#         inPath (str): file directory containing change in binding free energy. 
# 
#     """
#     ddg = pd.DataFrame(columns=["protein",
#                                 "partner",
#                                 "protein_pos",
#                                 "chain_pos",
#                                 "pdb_id",
#                                 "chain_mutation",
#                                 "partner_chain",
#                                 "submitted",
#                                 "DDG"])
#     c = -1
#     with io.open(inPath, "r", encoding="utf-8") as f:
#         next(f)
#         for line in f:
#             strsplit = list( map ( str.strip, line.split('\t') ) )
#             if len(strsplit) == 9:
#                 c += 1
#                 ddg.loc[c] = strsplit
#     ddg = ddg.drop_duplicates(subset = ["protein", "partner", "protein_pos"],
#                               keep = 'first').reset_index(drop=True)
#     ddg["protein_pos"] = ddg["protein_pos"].apply(int)
#     ddg["chain_pos"] = ddg["chain_pos"].apply(int)
#     ddg["DDG"] = ddg["DDG"].apply(float)
#     return ddg
