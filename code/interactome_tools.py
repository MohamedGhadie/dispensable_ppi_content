import pickle
import pandas as pd
import numpy as np
import networkx as nx
from text_tools import write_fasta_file
from simple_tools import (str_to_tuples,
                          str_to_list,
                          str_to_list_new,
                          list_to_str_new)

def create_ppi_graph (inPath):
    
    with open(inPath, 'rb') as f:
        partners = pickle.load(f)
    g = nx.Graph()
    proteins = list(partners.keys())
    g.add_nodes_from(proteins)
    for p in proteins:
        neighbors = partners[p]
        for nb in neighbors:
            g.add_edge(p, nb)
    return g

def find_all_shortest_paths (g, outPath):
    
    paths = dict(nx.all_pairs_shortest_path(g))
    with open(outPath, 'wb') as f:
        pickle.dump(paths, f)

def between (paths, p):
    subPaths = paths[paths.apply(lambda x: (x[0] != p) and (x[-1] != p))]
    return sum(subPaths.apply(lambda x: p in x[1:-1])) / len(subPaths)
    
def shortest_path (p1, p2, partners, path=[]):
    
    if p1 == p2:
        return [p1]
    else:
        minPath = []
        minLen = np.inf
        neighbors = partners[p1]
        for nb in neighbors:
            if nb not in path:
                pathExtension = [p1] + shortest_path(nb, p2, partners, path + [p1])
                if 1 < len(pathExtension) < minLen:
                    minPath = pathExtension
                    minLen = len(minPath)
        return minPath

def write_chain_annotated_interactome_to_excel (interactome,
                                                outPath,
                                                sheet_name = 'Sheet1'):
    
    write_chain_annotated_interactome ( interactome,
                                        outPath,
                                        toExcel = True,
                                        sheet_name = sheet_name )

def write_chain_annotated_interactome ( interactome,
                                        outPath,
                                        toExcel = False,
                                        sheet_name = 'Sheet1' ):
    
    interactomeCp = interactome.copy()
    interactomeCp["Mapping_chains"] = interactomeCp["Mapping_chains"].apply(lambda x:
                                                                            ', '.join( ['-'.join(pair)
                                                                                        for pair in x] ))
    if toExcel:
        writer = pd.ExcelWriter( str( outPath ) )
        interactomeCp.to_excel ( writer, sheet_name )
        writer.save()
    else:
        interactomeCp.to_csv ( outPath, index=False, sep='\t' )

def read_chain_annotated_interactome (inPath):
    """Read interactome annotated with chain pairs from a file.

    Args:
        inPath (str): file directory containing chain-annotated interactome.
    
    Returns:
        DataFrame: interactome with mapping chain pairs in tuple form.

    """
    interactome = pd.read_table(inPath, sep='\t')
    interactome["Mapping_chains"] = interactome["Mapping_chains"].apply( str_to_tuples )
    return interactome

def write_unmerged_interface_annotated_interactome ( interactome, outPath ):
    
    interactomeCp = interactome.copy()
    if "Mapping_frac" in interactomeCp.columns:
        interactomeCp["Mapping_frac"] = interactomeCp["Mapping_frac"].apply( lambda x:
                                                                             list_to_str_new( x, ['|', '+', '/'] ) )
    
    write_interface_annotated_interactome ( interactomeCp,
                                            outPath,
                                            delm = ['|', '+', '/', ','],
                                            chain_delm = ['|', '+', '/'] )

def read_unmerged_interface_annotated_interactome ( inPath ):
    """Read interactome annotated with interface residues from a file using a list of
        specified delimiters.

    Args:
        inPath (str): file directory containing interface-annotated interactome.
    
    Returns:
        DataFrame: interactome with interfaces in list form.

    """
    delm = ['|', '+', '/', ',']
    interactome = read_interface_annotated_interactome ( inPath, delm = delm)
    if "Mapping_frac" in interactome.columns:
        interactome["Mapping_frac"] = interactome["Mapping_frac"].apply(lambda x: 
                                                                        str_to_list_new(x, delm[ : -1 ], float))
    return interactome

def write_single_interface_annotated_interactome_to_excel (interactome,
                                                           outPath,
                                                           sheet_name = 'Sheet1'):
    
    write_interface_annotated_interactome_to_excel ( interactome,
                                                     outPath,
                                                     delm = ['+', ','],
                                                     sheet_name = sheet_name )

def write_single_interface_annotated_interactome (interactome, outPath):
    
    write_interface_annotated_interactome ( interactome,
                                            outPath,
                                            delm = ['+', ','] )

def read_single_interface_annotated_interactome ( inPath ):
    """Read interactome annotated with single interfaces.

    Args:
        inPath (str): file directory containing interface-annotated interactome.
    
    Returns:
        DataFrame: interactome with interfaces in tuple form.

    """
    interactome = read_interface_annotated_interactome ( inPath )
    interactome["Interfaces"] = interactome["Interfaces"].apply( lambda x: x[ 0 ] )
    return interactome

def write_interface_annotated_interactome_to_excel ( interactome,
                                                     outPath,
                                                     delm = None,
                                                     sheet_name = 'Sheet1' ):
    
    write_interface_annotated_interactome ( interactome,
                                            outPath,
                                            delm = delm,
                                            toExcel = True,
                                            sheet_name = sheet_name )

def write_interface_annotated_interactome ( interactome,
                                            outPath,
                                            delm = None,
                                            chain_delm = None,
                                            toExcel = False,
                                            sheet_name = 'Sheet1' ):
    
    if delm is None:
        delm = ['|', '+', ',']
    if chain_delm is None:
        chain_delm = ['|', '+']
    interactomeCp = interactome.copy()
    interactomeCp["Interfaces"] = interactomeCp["Interfaces"].apply( lambda x:
                                                                     list_to_str_new( x, delm ) )
    if "Chain_pairs" in interactomeCp.columns:
        interactomeCp["Chain_pairs"] = interactomeCp["Chain_pairs"].apply( lambda x:
                                                                           list_to_str_new( x, chain_delm )
                                                                           if type(x) == list
                                                                           else x )
    if toExcel:
        writer = pd.ExcelWriter( str( outPath ) )
        interactomeCp.to_excel ( writer, sheet_name )
        writer.save()
    else:
        interactomeCp.to_csv( outPath, index=False, sep='\t' )

def read_interface_annotated_interactome ( inPath, delm = None):
    """Read interactome annotated with interface residues from a file using a list of
        specified delimiters.

    Args:
        inPath (str): file directory containing interface-annotated interactome.
    
    Returns:
        DataFrame: interactome with interfaces in list form.

    """
    if delm is None:
        delm = ['|', '+', ',']
    interactome = pd.read_table(inPath, sep='\t')
    interactome["Interfaces"] = interactome["Interfaces"].apply(lambda x: 
                                                                str_to_list_new(x,
                                                                                delm,
                                                                                int))
    if "Chain_pairs" in interactome.columns:
        interactome["Chain_pairs"] = interactome["Chain_pairs"].apply(lambda x: 
                                                                        str_to_list_new(x,
                                                                                        delm[ : -1 ],
                                                                                        str))
    return interactome
    
def write_interactome_sequences (inPath,
                                 sequenceFile,
                                 outPath):
    """Write interactome protein sequences to fasta file.

    Args:
        inPath (str): file directory containing interactome.
        sequenceFile (str): file directory containing protein sequences.
        outPath(str): file directory to write protein sequences to.

    """
    interactome = pd.read_table(inPath, sep='\t')
    sequences = pd.read_table(sequenceFile, sep='\t')
    proteins = list(set(interactome[["Protein_1", "Protein_2"]].values.flatten()))
    interactomeSequences = pd.DataFrame(index=range(len(proteins)),
                                        columns=['ID', 'Sequence'])
    interactomeSequences["ID"] = proteins
    interactomeSequences["Sequence"] = interactomeSequences["ID"].apply(lambda x:
                                                                        sequences.loc[sequences["ID"]==x,
                                                                                      "Sequence"].item())
    write_fasta_file(interactomeSequences,
                     "ID",
                     "Sequence",
                     outPath)

def remove_interactions_reported_once(interactions):
    """Remove protein-protein interactions occuring only one time.

    Args:
        interactions (DataFrame):  protein-protein interactions.
    
    Returns:
        DataFrame: protein-protein interactions with at least one duplicate.
        
    """
    # remove interactions appearing only once
    lessthan = interactions["Protein_2"] < interactions["Protein_1"]
    list1 = list(interactions.loc[lessthan, "Protein_1"].values)
    list2 = list(interactions.loc[lessthan, "Protein_2"].values)
    interactions.loc[lessthan, "Protein_1"] = list2
    interactions.loc[lessthan, "Protein_2"] = list1
    duplicates = interactions.duplicated(subset=["Protein_1", "Protein_2"], keep=False)
    interactions.loc[lessthan, "Protein_1"] = list1
    interactions.loc[lessthan, "Protein_2"] = list2
    duplicateInteractions = interactions[duplicates].reset_index(drop=True)
    return duplicateInteractions
    
def duplicated_PPIs (PPIs, minCount = 2):
    """Return protein-protein interactions occuring at least minOccur times.

    Args:
        PPIs (DataFrame): protein-protein interactions.
        minCount (int): minimum number of occurrences for retained PPIs.
    
    Returns:
        DataFrame: PPIs occuring at least minOccur times.
        
    """
    PPIstrings = PPIs.apply(lambda x: '-'.join( sorted( [ x["Protein_1"], x["Protein_2"] ] ) ), axis=1)
    counts = PPIstrings.apply(lambda x: sum(PPIstrings == x))
    return PPIs[counts >= minCount].reset_index(drop=True)

def remove_duplicate_PPIs (PPIs):
    """Remove duplicate protein-protein interactions (PPIs), keeping first occurence.

    Args:
        PPIs (DataFrame):  protein-protein interactions.
    
    Returns:
        DataFrame: PPIs with no duplicates.
        
    """
    PPIstrings = PPIs.apply(lambda x: '-'.join( sorted( [ x["Protein_1"], x["Protein_2"] ] ) ), axis=1)
    duplicate = PPIstrings.duplicated()
    return PPIs[duplicate == False].reset_index(drop=True)

def sample_noninteracting_pairs(inPath, sampleSize, outPath):
    """Select a sample of non-interacting protein pairs.

    Args:
        inPath (str): file directory containing an interactome.
        sampleSize (int): number of pairs to sample.
        outPath (str): file directory to save non-interacting protein pairs to.
    
    Returns:
        DataFrame: pairs of non-interacting proteins.

    """   
    interactome = pd.read_table(inPath, sep='\t')
    proteins = list(set(interactome[["Protein_1", "Protein_2"]].values.flatten()))
    numProteins = len(proteins)
    PPIindex = pd.DataFrame(columns=("Protein_1","Protein_2"))
    PPIindex["Protein_1"] = interactome["Protein_1"].apply(lambda x: proteins.index(x))
    PPIindex["Protein_2"] = interactome["Protein_2"].apply(lambda x: proteins.index(x))
    PPIs = set([tuple(np.sort([row.Protein_1, row.Protein_2])) for _, row in PPIindex.iterrows()])
    nonPPIs = list()
    for p1 in range(numProteins):
        nonPPIs.extend([(p1,p2) for p2 in range(p1,numProteins) if not ((p1,p2) in PPIs)])
    sample = np.random.choice(len(nonPPIs), size=sampleSize, replace=False)
    nonPPIsample = pd.DataFrame(index=range(len(sample)), columns=("Protein_1","Protein_2"))
    c = -1
    for i in sample:
        pair = nonPPIs[i]
        c += 1
        nonPPIsample.loc[c,:] = proteins[pair[0]], proteins[pair[1]]
    nonPPIsample.to_csv(outPath, index=False, sep='\t')
