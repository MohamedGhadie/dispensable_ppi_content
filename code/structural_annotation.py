import os
import io
import time
import pickle
import pandas as pd
import numpy as np
from pathlib import Path
from random import sample
from Bio.PDB import PDBParser
from interactome_tools import (write_chain_annotated_interactome,
                               read_chain_annotated_interactome,
                               write_unmerged_interface_annotated_interactome,
                               read_unmerged_interface_annotated_interactome,
                               write_single_interface_annotated_interactome,
                               read_interface_annotated_interactome,
                               write_interface_annotated_interactome)
from simple_tools import (create_dir,
                          isolate_pairs,
                          merge_list_pairs)
from text_tools import (read_list_table,
                        write_hpc_job)
from pdb_tools import (clear_structures,
                       write_partial_structure,
                       load_chainResOrder,
                       get_interface_by_chainIDs)

known_interfaces = {}
knownCCIs = newCCIs = numPPIs = PPIcount = annotatedPPIs = 0

def clear_globals ():
    """Clear global variables
    
    """
    global knownCCIs, newCCIs, numPPIs, PPIcount, annotatedPPIs
    print('\tclearing global variables used in structural annotation')
    knownCCIs = newCCIs = numPPIs = PPIcount = annotatedPPIs = 0

def clear_interfaces ():
    """Clear dict of known chain interfaces
    
    """
    global known_interfaces
    print('\tclearing known interfaces dictionary')
    known_interfaces.clear()

def produce_alignment_evalue_dict (inPath,
#                                    proteinChainsFile,
                                   outPath,
                                   method = 'min'):
    
    chainMap = pd.read_table(inPath, sep='\t')
    print(len(chainMap))
    time.sleep(10)
#     with open(proteinChainsFile, 'rb') as f:
#         proteinChains = pickle.load(f)
    evalueDict = {}
#     for i, p in enumerate(proteinChains.keys()):
#         print(i)
#         proteinRows = (chainMap['Query'] == p)
#         for j, c in enumerate(proteinChains[p]):
#             print('\t%d' % j)
#             evalues = chainMap.loc[ proteinRows & (chainMap['Subject'] == c), "Expect"].values
#             if len(evalues) == 1:
#                 evalueDict [p + '-' + c] = evalues[0]
#             elif method == 'min':
#                 evalueDict [p + '-' + c] = min( evalues )
#             elif method == 'mean':
#                 evalueDict [p + '-' + c] = np.mean( evalues )
    if method == 'min':
        for i, row in chainMap.iterrows():
            print(i)
            k = row.Query + '-' + row.Subject
            if k in evalueDict:
                if (row.Expect < evalueDict[k]):
                    evalueDict[k] = row.Expect
            else:
                evalueDict[k] = row.Expect
    elif method == 'mean':
        for i, row in chainMap.iterrows():
            print(i)
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
                                 
def map_positions (refPos,
                   alignPos,
                   pos):
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
    alignmentEvalues = None
    if alignmentEvalueFile is not None:
        with open(alignmentEvalueFile, 'rb') as f:
            alignmentEvalues = pickle.load(f)
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
        
        if alignmentEvalues is not None:
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

def produce_interface_annotated_interactome (inPath,
                                             pdbDir,
                                             chainMapFile,
                                             interfaceFile,
                                             chainResOrderFile,
                                             maxInterfaces,
                                             maxAttempts,
                                             rnd,
                                             mapCutoff,
                                             bindingDist,
                                             outPath):
    """Map interaction interfaces onto protein-protein interaction network.

    Args:
        inPath (str): file directory containing protein interaction network.
        pdbDir (str): directory containing pdb structure files.
        chainMapFile (str): file directory containing table of chains mapping onto each protein.
        interfaceFile (str): file directory containing known interfaces of each chain with its
                                binding chain. Newly calculated interfaces are added to this file.
        chainResOrderFile (str): file directory containing dict of residue order in list
                                    of residues with known coordinates.
        numinf (int): number of chaninPairs to use for interface mapping.
        rnd (boolean): True to randomly select chain pairs from PPI annotation list 
                        for interface mapping.
        outPath (str): file directory to save interface-annotated interactome to.

    """
    global known_interfaces, knownCCIs, newCCIs, numPPIs, PPIcount
    clear_globals ()
    clear_structures()
    clear_interfaces()
    load_chainResOrder(chainResOrderFile)
    print('\treading chain-annotated interactome')
    interactome = read_chain_annotated_interactome(inPath)
    print('\treading protein-chain mapping table')
    chainMap = read_list_table(chainMapFile, ["Qpos", "Spos"], [int, int], '\t')
    numPPIs = len(interactome)
    
    if interfaceFile.is_file():
        print('\tloading known chain interfaces')
        read_chain_interfaces(interfaceFile)
    else:
        print('creating interface file directory')
        filedir, _ = os.path.split(interfaceFile)
        create_dir(filedir)
    
    print('\tloading chain residue order dictionary')
    with open(chainResOrderFile, 'rb') as f:
        resOrder = pickle.load(f)
    print('\tresolving interfaces for all %d PPIs' % numPPIs)
    interactome["Interfaces"] = interactome.apply(lambda x:
                                                  map_multi_interfaces(pdbDir,
                                                                       x["Protein_1"],
                                                                       x["Protein_2"],
                                                                       x["Mapping_chains"],
                                                                       chainMap,
                                                                       resOrder,
                                                                       maxInterfaces,
                                                                       maxAttempts,
                                                                       rnd,
                                                                       mapCutoff,
                                                                       bindingDist,
                                                                       interfaceFile),
                                                  axis=1)
    
    
    print('\tfinalizing interface-annotated interactome')
    interactome = interactome[interactome["Interfaces"].apply(lambda x:
                                                              len(x[0]) > 0)
                             ].reset_index(drop=True)
    interactome.drop("Mapping_chains", axis=1, inplace=True)
    interactome["Chain_pairs"] = interactome["Interfaces"].apply(lambda x: x[2])
    interactome["Mapping_frac"] = interactome["Interfaces"].apply(lambda x: x[1])
    interactome["Interfaces"] = interactome["Interfaces"].apply(lambda x: x[0])
    write_unmerged_interface_annotated_interactome (interactome, outPath)
    print('\t%d known CCIs' % knownCCIs)
    print('\t%d new CCIs' % newCCIs)

def map_multi_interfaces (pdbDir,
                          protein1,
                          protein2,
                          chainPairs,
                          chainMap,
                          resOrder,
                          maxInterfaces,
                          maxAttempts,
                          rnd,
                          mapCutoff,
                          bindingDist,
                          interfaceFile):
    """Map the first n interaction interfaces passing a mapping fraction cutoff for two 
        interacting proteins from mapping chain pairs from distinct PDB models.

    Args:
        pdbDir (str): directory containing pdb structure files.
        protein1 (str): ID of first interacting protein.
        protein2 (str): ID of second interacting protein.
        chainPairs (list): list of chain pairs mapping onto the two interacting proteins.
        chainMap (dataframe): table of chains mapping onto each protein.
        resOrder (dict): dict of residue order in list of residues with known coordinates.
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
    global numPPIs, PPIcount, annotatedPPIs
    PPIcount += 1
    print('\t' + 'PPI: ' + protein1 + ' - ' + protein2)
    
    numAttempts = interfacesFound = 0
    numChainPairs = len(chainPairs)
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
        print('\t\t' + 'PPI: ' + protein1 + ' - ' + protein2 + ' (# %d out of %d) (%.1f %%)' 
                % (PPIcount, numPPIs, 100 * PPIcount / numPPIs))
        print('\t\t' + '%d attempts (max %d)' % (numAttempts, maxAttempts))
        print('\t\t' + '%d interfaces mapped (seeking %d)' % (interfacesFound, maxInterfaces))
        print('\t\t' + 'Chain pair: ' + chainID1 + ' - ' + chainID2 + ' (# %d out of %d)' % (i + 1, numChainPairs))
        if pdbID not in sel_pdbs:
            if (chainID1 in resOrder) and (chainID2 in resOrder):
                interfaces1, interfaces2 = map_twoside_interfaces (pdbDir,
                                                                   protein1,
                                                                   protein2,
                                                                   chainPair,
                                                                   bindingDist,
                                                                   chainMap,
                                                                   interfaceFile)
#                 if (len(interfaces1) > 0):
#                     ch1interfaces, ch2interfaces, ch1frac, ch2frac = interfaces
#                     passes1 = [(i, f) for i, f in zip(ch1interfaces, ch1frac) if f >= mapCutoff]
#                     passes2 = [(i, f) for i, f in zip(ch2interfaces, ch2frac) if f >= mapCutoff]
                passFrac1 = [f for _, f in interfaces1 if f >= mapCutoff]
                passFrac2 = [f for _, f in interfaces2 if f >= mapCutoff]
                if (len(passFrac1) > 0) and (len(passFrac2) > 0):
                    inf1 = [i for i, f in interfaces1 if f >= mapCutoff]
                    inf2 = [i for i, f in interfaces2 if f >= mapCutoff]
#                         passInf1, passFrac1 = zip(*passes1)
#                         passInf2, passFrac2 = zip(*passes2)
                    sel_interfaces.append((inf1, inf2))
                    sel_frac.append((passFrac1, passFrac2))
                    sel_chainpairs.append(([chainID1]*len(inf1), [chainID2]*len(inf2)))
                    sel_pdbs.append(pdbID)
                    interfacesFound += 1
                    print('\t\t\t' + 'interface %d found' % interfacesFound)
                numAttempts += 1
            else:
                print('\t\t\t' + 'chain residue order not known')
        else:
            print('\t\t\t' + 'interface already mapped from PDB structure')
    
    annotatedPPIs += (len(sel_chainpairs) > 0)
    print('\t' + '%d PPIs resolved' % PPIcount)
    print('\t' + '%d PPIs annotated' % annotatedPPIs)
    return (sel_interfaces, sel_frac, sel_chainpairs)

def map_twoside_interfaces (pdbDir,
                            protein1,
                            protein2,
                            chainPair,
                            bindingDist,
                            chainMap,
                            interfaceFile):
    global known_interfaces, knownCCIs, newCCIs
    
    chain1, chain2 = chainPair
    chainKey1 = chain1 + '-' + chain2
    chainKey2 = chain2 + '-' + chain1
    if (chainKey1 not in known_interfaces) or (chainKey2 not in known_interfaces):
        newCCIs += 1
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
        knownCCIs += 1
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
#         interfaces1 = list( interfaces1[ interfaces1.apply(lambda x: x[1] > 0) ].values )
#         interfaces2 = list( interfaces2[ interfaces2.apply(lambda x: x[1] > 0) ].values )
#         interfaces1 = interfaces1.values
#         interfaces2 = interfaces2.values 
#         frac1 = [f for _, f in interfaces1 if f > 0]
#         interfaces1 = [i for i, f in interfaces1 if f >0]
#         frac2 = [f for _, f in interfaces2 if f > 0]
#         interfaces2 = [i for i, f in interfaces2 if f > 0]
#         return (interfaces1, interfaces2, frac1, frac2)
        interfaces1 = [(i, f) for i, f in interfaces1.values if f > 0]
        interfaces2 = [(i, f) for i, f in interfaces2.values if f > 0]
        return (interfaces1, interfaces2)
    else:
        return ([],[])

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
    #interactome["Mapping_frac"] = uniqueInterfaces.apply(lambda x: x[1])
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

def write_pdb_mapped_mutations (mutations,
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

def read_mutation_ddg (inPath):
    
    """Read mutation change (loss) in binding free energy from file writen by 
    function "write_pdb_mapped_mutations".

    Args:
        inPath (str): file directory containing change in binding free energy. 

    """
    ddg = pd.DataFrame(columns=["protein",
                                "partner",
                                "protein_pos",
                                "chain_pos",
                                "pdb_id",
                                "chain_mutation",
                                "partner_chain",
                                "submitted",
                                "DDG"])
    c = -1
    with io.open(inPath, "r", encoding="utf-8") as f:
        next(f)
        for line in f:
            strsplit = list( map ( str.strip, line.split('\t') ) )
            if len(strsplit) == 9:
                c += 1
                ddg.loc[c] = strsplit
    ddg = ddg.drop_duplicates(subset = ["protein", "partner", "protein_pos"],
                              keep = 'first').reset_index(drop=True)
    ddg["protein_pos"] = ddg["protein_pos"].apply(int)
    ddg["chain_pos"] = ddg["chain_pos"].apply(int)
    ddg["DDG"] = ddg["DDG"].apply(float)
    return ddg

def copy_mutation_ddg (inPath1, inPath2, outPath):
    
    """Copy mutation change in binding free energy from one file to another.

    Args:
        inPath1 (str): file directory containing change in binding free energy.
        inPath2 (str): file directory containing mutations with unknown change in binding free energy.
        outPath (str): file directory to output mutations with change in binding free energy.

    """
    ddg = {}
    with io.open(inPath1, "r", encoding="utf-8") as f:
        next(f)
        for line in f:
            strsplit = list( map ( str.strip, line.split('\t') ) )
            if len(strsplit) > 7:
                pdbid, mut, partnerChain = strsplit[ 4 : 7 ]
                k = (pdbid,) + tuple(sorted([mut[1], partnerChain])) + (mut,)
                ddg[k] = strsplit[-1]
    write_mutation_ddg_tofile (ddg, inPath2, outPath)

def write_mutation_ddg_tofile (ddg, inPath, outPath):
    
    with io.open(inPath, "r", encoding="utf-8") as f, io.open(outPath, "w") as fout:
        for line in f:
            strsplit = list( map ( str.strip, line.split('\t') ) )
            if len(strsplit) >= 7:
                pdbid, mut, partnerChain = strsplit[ 4 : 7 ]
                k = (pdbid,) + tuple(sorted([mut[1], partnerChain])) + (mut,)
                if len(strsplit) == 7:
                    if k in ddg:
                        if ddg[k] in ['X', 'S']:
                            strsplit.append( ddg[k] )
                        else:
                            strsplit.append('S')
                            strsplit.append( ddg[k] )
                elif strsplit[-1] == 'S':
                    if k in ddg:
                        if ddg[k] not in ['X', 'S']:
                            strsplit.append( ddg[k] )
            fout.write( '\t'.join(map(str, strsplit)) +  '\n' )

def produce_bindprofx_jobs (mutations,
                            pdbDir,
                            outDir,
                            write_hpc_jobfiles = True,
                            nodes = 1,
                            ppn = 1,
                            pmem = 7700,
                            walltime = '1:00:00:00',
                            rapid = None,
                            username = '',
                            hpcCommands = None,
                            serverDataDir = '../data'):
    
    clear_structures()
    
    dataDir = outDir / 'data'
    jobDir = outDir / 'jobs'
    if not dataDir.exists():
        os.makedirs(dataDir)
    if not jobDir.exists():
        os.makedirs(jobDir)
    
    for struc, mutList in mutations.items():
        strucid = '_'.join(struc)
        strucDir = dataDir / strucid
        if not strucDir.exists():
            os.makedirs(strucDir)
        mutListFile = strucDir / 'mutList.txt'
        mutList = [ '%s;' % mutList.pop(0) ] + ['\n%s;' % mut for mut in mutList]
        with io.open(mutListFile, "w") as fout:
            for mut in mutList:
                fout.write(mut)
        pdbid, chainID1, chainID2 = struc[:3]
        write_partial_structure (pdbid,
                                 [chainID1, chainID2],
                                 pdbDir,
                                 strucDir / 'complex.pdb')
        
        if write_hpc_jobfiles:
            commands = [ '../bin/get_final_score.py %s/%s' % (serverDataDir, strucid) ]
            if hpcCommands:
                commands = hpcCommands + commands
            write_hpc_job (jobDir / (strucid + '_job.txt'),
                           nodes = nodes,
                           ppn = ppn,
                           pmem = pmem,
                           walltime = walltime,
                           outputfile = '%s/%s/outputfile' % (serverDataDir, strucid),
                           errorfile = '%s/%s/errorfile' % (serverDataDir, strucid),
                           rapid = rapid,
                           jobid = username + '_bindprofx_' + strucid,
                           commands = commands)
