#----------------------------------------------------------------------------------------
# Modules that produce mappings of IDs across different databases.
#----------------------------------------------------------------------------------------

import pickle
import pandas as pd
from pathlib import Path
from Bio import SeqIO
from Bio.SubsMat import MatrixInfo
from text_tools import read_list_table

def produce_uniqueGene_swissProtIDs (inPath, geneMapFile, outPath):
    """Remove Swiss-Prot IDs with either no gene name or a gene name shared 
        with other Swiss-Prot IDs.

    Args:
        inPath (Path): path to file containing list of Swiss-Prot IDs.
        geneMapFile (Path): path to file containing dictionary of uniProt ID gene names.
        outPath (Path): path to save list of Swiss-Prot IDs with unique gene names.

    """
    with open(inPath, 'r') as f:
        swissProtIDs = pd.Series(list(set(f.read().split())))
    with open(geneMapFile, 'rb') as f:
        geneMap = pickle.load(f)
    swissProtIDGenes = swissProtIDs.apply(lambda x: geneMap[x] if x in geneMap else 'no_name')
    duplicates = swissProtIDGenes.duplicated(keep=False)
    uniqueGeneSwissProtIDs = swissProtIDs[(duplicates==False) & (swissProtIDGenes!='no_name')]
    uniqueGeneSwissProtIDs.to_csv(outPath, index=False)

def produce_uniqueGene_sequences (inPath,
                                  uniqueGeneSwissProtIDFile,
                                  geneMapFile,
                                  outPath):
    """Remove sequences with either no gene name or a gene name shared with other sequences.

    Args:
        inPath (Path): path to file containing proteome sequences
        uniqueGeneSwissProtIDFile (Path): path to file containing list of Swiss-Prot IDs with unique gene names
        geneMapFile (Path): file path containing dictionary of Swiss-Prot ID gene names
        outPath (Path): file path to save protein sequences with unique gene names

    """
    with open(uniqueGeneSwissProtIDFile, 'r') as f:
        uniqueGeneSwissProtIDs = set(f.read().split())
    with open(geneMapFile, 'rb') as f:
        geneMap = pickle.load(f)
    sequences = pd.read_table(inPath, sep="\t")
    sequencesRefID = sequences["ID"].apply(lambda x: x if x.find('-') == -1 else x[:x.find('-')])
    sequences["Gene_Name"] = sequencesRefID.apply(lambda x: geneMap[x] if x in geneMap else x)
    sequences = sequences[sequencesRefID.isin(uniqueGeneSwissProtIDs)].reset_index(drop=True)
    sequences.to_csv(outPath, index=False, sep='\t')

def produce_proteinSeq_dict (inPath, outPath):
    """Make a dictionary of sequences.

    Args:
        inPath (Path): path to sequence file in fasta format.
        outPath (Path): path to save output pickle file to.

    """
    s = list(SeqIO.parse(str(inPath), 'fasta'))
    proteinSeq = {}
    for _, elm in enumerate(s):
        proteinSeq[elm.id] = str(elm.seq)
    with open(outPath, 'wb') as fOut:
        pickle.dump(proteinSeq, fOut)

def produce_uniprotID_dict(inPath, spList, outPath):
    """Make a map of various gene IDs to UniProt IDs.

    Only use Swiss-Prot proteins.

    Args:
        inPath (Path): path to UniProt ID mapping file.
        spList (Path): path to Swiss-Prot ID list file.
        outPath (Path): path to save output pickle file to.

    """
    with open(spList, 'r') as f:
        swissProtIDs = set(f.read().split())
    with open(inPath, 'r') as fIn:
        idMap = {}
        for line in fIn:
            uniprotID, otherIDtype, otherID = line.strip().split('\t')
            if uniprotID in swissProtIDs:
                if otherIDtype == 'Gene_Name':
                    otherID = otherID.upper()
                idMap[otherID] = uniprotID
    with open(outPath, 'wb') as fOut:
        pickle.dump(idMap, fOut)

def produce_geneName_dict(inPath, spList, outPath):
    """Make a map of various UniProt IDs to gene names.

    Only use Swiss-Prot proteins.

    Args:
        inPath (Path): path to UniProt ID mapping file.
        spList (Path): path to Swiss-Prot ID list file.
        outPath (Path): path to save output pickle file to.

    """
    with open(spList, 'r') as f:
        swissProtIDs = set(f.read().split())
    with open(inPath, 'r') as fIn:
        idMap = {}
        for line in fIn:
            uniprotID, otherIDtype, otherID = line.strip().split('\t')
            if otherIDtype == 'Gene_Name':
                if uniprotID in swissProtIDs:
                    idMap[uniprotID] = otherID.upper()
    with open(outPath, 'wb') as fOut:
        pickle.dump(idMap, fOut)

def produce_isoform_geneName_dict(geneMapFile, isoformFile, outPath):
    """Make a map of isoform UniProt IDs to gene names.

    Args:
        geneMapFile (Path): path to file containing gene name dict for reference UniProt IDs.
        isoformFile (Path): path to file containing list of all isoform Swiss-Prot IDs.
        outPath (Path): path to save gene name dict for all isoforms. 

    """
    isoformData = pd.read_table(isoformFile, sep="\t")
    with open(geneMapFile, 'rb') as f:
        geneMap = pickle.load(f)
    isoformGeneMap = {}
    isoformData["refID"] = isoformData["Isoform"].apply(lambda x: x if x.find('-') == -1 else x[:x.find('-')])
    for _, row in isoformData.iterrows():
        if row.refID in geneMap:
            isoformGeneMap[row.Isoform] = geneMap[row.refID]
    with open(outPath, 'wb') as fOut:
        pickle.dump(isoformGeneMap, fOut)

def produce_protein_interaction_dict (inPath, outPath):
    """Make a dictionary of protein interaction partners.

    Args:
        inPath (Path): path to tab-delimited file containing protein-protein interactions.
        outPath (Path): path to save output pickle file to.

    """ 
    PPIs = pd.read_table(inPath, sep="\t")
    proteins =  set(PPIs[["Protein_1", "Protein_2"]].values.flatten())
    proteinPartners = {}
    for protein in proteins:
        partners = set(PPIs.loc[(PPIs[["Protein_1", "Protein_2"]]==protein).any(1),
                                ["Protein_1", "Protein_2"]].values.flatten()) - {protein}
        if sum((PPIs[["Protein_1", "Protein_2"]]==protein).all(1)) > 0:
            partners.add(protein)
        proteinPartners[protein] = partners
    with open(outPath, 'wb') as fOut:
        pickle.dump(proteinPartners, fOut)

def produce_protein_chain_dict (inPath, outPath):
    """Make a dictionary of protein mapping PDB chains.

    Args:
        inPath (Path): path to tab-delimited file containing protein-chain mapping.
        outPath (Path): path to save output pickle file to.

    """
    chainMap = pd.read_table(inPath, sep="\t")
    proteins =  set(chainMap["Query"])
    proteinChains = {}
    for protein in proteins:
        proteinChains[protein] = set(chainMap.loc[chainMap["Query"]==protein, "Subject"])
    with open(outPath, 'wb') as fOut:
        pickle.dump(proteinChains, fOut)

def produce_protein_chain_alignment_dict (inPath, outPath):
    """Make a dictionary of protein-chain alignments.

    Args:
        inPath (Path): path to file containing protein-chain alignment table.
        outPath (Path): path to save output pickle file to.

    """
    chainMap = read_list_table(inPath, ["Qpos", "Spos"], [int, int], '\t')
    chainMapDict = {}
    for _, row in chainMap.iterrows():
        k = row.Query + '-' + row.Subject
        if k in chainMapDict:
            chainMapDict[k].append((row.Qpos, row.Spos))
        else:
            chainMapDict[k] = [(row.Qpos, row.Spos)]
    with open(outPath, 'wb') as fOut:
        pickle.dump(chainMapDict, fOut)

def produce_chain_dict (inPath, outPath):
    """Produce a dictionary of chain IDs associated with each PDB ID. 

    Args:
        inPath (Path): path to file containing list of chain IDs.
        outPath (Path): path to save pickle file to.

    """
    with open(inPath, 'r') as fin:
        chainIDs = list(fin.read().split())
    chains = {}
    for chainid in chainIDs:
        pdbid = (chainid[ : chainid.find('_') ] if '_' in chainid else chainid)
        if pdbid in chains:
            chains[pdbid].add(chainid)
        else:
            chains[pdbid] = {chainid}
    with open(outPath, 'wb') as fOut:
        pickle.dump(chains, fOut)

def produce_chain_strucRes_dict (inPath, outPath):
    """Produce a dictionary of labels for chain residue positions associated with 3D coordinates. 

    Args:
        inPath (Path): path to file containing PDB chain residue labels.
                        See https://cdn.rcsb.org/etl/kabschSander/ss_dis.txt.gz
        outPath (Path): path to save pickle file to.

    """
    s = list( SeqIO.parse( str(inPath), 'fasta') )
    strucRes = {}
    for row in s:
        if ':disorder' in row.id:  
            pdbid, chainID, _ = list( map( str.strip, row.id.split(':') ) )
            strucRes[ pdbid.lower() + '_' + chainID ] = str( row.seq )
    with open(outPath, 'wb') as fOut:
        pickle.dump(strucRes, fOut)

def produce_chainSeq_dict (inPath, outPath):
    """Make a dictionary of PDB chain sequences.

    Args:
        inPath (Path): path to chain sequence file in fasta format.
        outPath (Path): path to save output pickle file to.

    """
    s = list(SeqIO.parse(str(inPath), 'fasta'))
    chainSeq = {}
    for _, elm in enumerate(s):
        chainSeq[elm.id] = str(elm.seq)
    with open(outPath, 'wb') as fOut:
        pickle.dump(chainSeq, fOut)

def produce_substitution_matrix (name, outPath):
    """Make a dictionary of mutation substitution scores.

    Args:
        name (str): name of substitution table to use.
        outPath (Path): path to save output pickle file to.

    """
    matrices = {'PAM30': MatrixInfo.pam30,
                'PAM60': MatrixInfo.pam60,
                'BLOSUM62' : MatrixInfo.blosum62,
                'BLOSUM95' : MatrixInfo.blosum95}
    matrix = matrices[name]
    scores = {}
    for k, v in matrix.items():
        scores[k] = v
        scores[tuple(reversed(k))] = v
    with open(outPath, 'wb') as fOut:
        pickle.dump(scores, fOut)

def produce_rnaToProtein_refseqID_dict (inPath, outPath):
    """Make a map of RNA RefSeq IDs to protein RefSeq IDs.

    Args:
        inPath (Path): path to RefSeq ID mapping file.
        outPath (Path): path to save output pickle file to.

    """
    idMap = {}
    with open(inPath, 'r') as f:
        next(f)
        for line in f:
            tax_id, gene_id, symbol, rsg, lrg, rna, t, protein, p, category = line.strip().split('\t')
            if (len(rna) > 0) and (len(protein) > 0):
                idMap[rna] = protein
    with open(outPath, 'wb') as fOut:
        pickle.dump(idMap, fOut)
