import io
import csv
import pickle
import pandas as pd
from pathlib import Path
from Bio import SeqIO
from Bio.SubsMat import MatrixInfo
from text_tools import (parse_fasta_file,
                        read_list_table)

def produce_uniqueGene_swissProtIDs(inPath, geneMapFile, outPath):
    """Remove Swiss-Prot IDs with either no gene name or a gene name shared 
    with other Swiss-Prot IDs.

    Args:
        inPath (str): file directory containing list of Swiss-Prot IDs
        geneMapFile (str): file directory containing dictionary of uniProt ID gene names
        outPath (str): file directory to save list of Swiss-Prot IDs with unique gene names

    """
    with open(inPath, 'r') as f:
        swissProtIDs = pd.Series(list(set(f.read().split())))
    with open(geneMapFile, 'rb') as f:
        geneMap = pickle.load(f)
    swissProtIDGenes = swissProtIDs.apply(lambda x: geneMap[x] if x in geneMap else 'no_name')
    duplicates = swissProtIDGenes.duplicated(keep=False)
    uniqueGeneSwissProtIDs = swissProtIDs[(duplicates==False) & (swissProtIDGenes!='no_name')]
    uniqueGeneSwissProtIDs.to_csv(outPath, index=False)

def produce_uniqueGene_sequences(inPath, uniqueGeneSwissProtIDFile, geneMapFile, outPath):
    """Remove sequences with either no gene name or a gene name shared 
    with other sequences.

    Args:
        inPath (str): file directory containing proteome sequences
        uniqueGeneSwissProtIDFile (str): file directory containing list of Swiss-Prot IDs with unique gene names
        geneMapFile (str): file directory containing dictionary of Swiss-Prot ID gene names
        outPath (str): file directory to save protein sequences with unique gene names

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
        inPath (str): human UniProt ID mapping file.
        outPath (str): path to save output pickle file to.

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
        inPath (str): human UniProt ID mapping file.
        outPath (str): path to save output pickle file to.

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
    """Make a map of isoform UniProt IDs to gene names

    Args:
        geneMapFile (str): file directory containing gene name dict for reference proteins
        isoformFile (str): file directory containing list of all isoform Swiss-Prot IDs
        outPath (str): file directory to save gene name dict for all isoforms 

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

def produce_domainName_dict(inPath, outPath):
    """Make a map of various Pfam IDs to domain names.

    Args:
        inPath (str): human Pfam hmm data file mapping file.
        outPath (str): path to save output pickle file to.

    """
    domainName = {}
    with io.open(inPath, "r", encoding="utf-8") as f:
        reader = csv.reader(f)
        for row in reader:
                if row[0][:7] == '#=GF ID':
                    name = row[0][7:].strip()
                elif row[0][:7] == '#=GF AC':
                    if name !='-':
                        acc = row[0][7:].strip()
                        dot = acc.find('.')
                        acc = acc[:dot] if dot!=-1 else acc
                        domainName[acc] = name
                        name = '-'
    with open(outPath, 'wb') as fOut:
        pickle.dump(domainName, fOut)

def produce_domainID_dict(inPath, outPath):
    """Make a map of various domain names to Pfam IDs.

    Args:
        inPath (str): human Pfam hmm data file mapping file.
        outPath (str): path to save output pickle file to.

    """
    domainID = {}
    with io.open(inPath, "r", encoding="utf-8") as f:
        reader = csv.reader(f)
        for row in reader:
                if row[0][:7] == '#=GF ID':
                    name = row[0][7:].strip()
                elif row[0][:7] == '#=GF AC':
                    if name !='-':
                        acc = row[0][7:].strip()
                        dot = acc.find('.')
                        acc = acc[:dot] if dot!=-1 else acc
                        domainID[name] = acc
                        name = '-'
    with open(outPath, 'wb') as fOut:
        pickle.dump(domainID, fOut)

def produce_domain_interaction_dict(inPath, outPath):
    """Make a dictionary of domain interaction partners.

    Args:
        inPath (str): file directory of domain-domain interactions.
        outPath (str): path to save output pickle file to.

    """
    DDIs = pd.read_table(inPath, sep="\t")
    domains =  set(DDIs[["dom1", "dom2"]].values.flatten())
    domainPartners = {}
    for domain in domains:
        partners = set(DDIs.loc[(DDIs[["dom1", "dom2"]]==domain).any(1),
                                ["dom1", "dom2"]].values.flatten()) - {domain}
        if sum((DDIs[["dom1", "dom2"]]==domain).all(1)) > 0:
            partners.add(domain)
        domainPartners[domain] = partners
    with open(outPath, 'wb') as fOut:
        pickle.dump(domainPartners, fOut)

def produce_protein_interaction_dict(inPath, outPath):
    """Make a dictionary of protein interaction partners.

    Args:
        inPath (str): file directory of protein-protein interactions.
        outPath (str): path to save output pickle file to.

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

def produce_protein_domain_dict(inPath, outPath):
    """Make a dictionary of protein mapping domains.

    Args:
        inPath (str): file directory of protein-domain mapping.
        outPath (str): path to save output pickle file to.

    """
    domMap = pd.read_table(inPath, sep="\t")
    proteins =  set(domMap["spid"])
    proteinDomains = {}
    for protein in proteins:
        proteinDomains[protein] = set(domMap.loc[domMap["spid"]==protein, "domain"])
    with open(outPath, 'wb') as fOut:
        pickle.dump(proteinDomains, fOut)

def produce_protein_chain_dict(inPath, outPath):
    """Make a dictionary of protein mapping chains.

    Args:
        inPath (str): file directory of protein-chain mapping.
        outPath (str): path to save output pickle file to.

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
        inPath (str): file directory of protein-chain mapping table.
        outPath (str): path to save output pickle file to.

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
    """Produce dict of chain IDs for each PDB ID 

    Args:
        inPath (str): file directory containing list of chain IDs
        outPath (str): file directory to save dict to

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

def produce_chainSeq_dict (inPath, outPath):
    
    s = list(SeqIO.parse(str(inPath), 'fasta'))
    chainSeq = {}
    for _, elm in enumerate(s):
        chainSeq[elm.id] = str(elm.seq)
    with open(outPath, 'wb') as fOut:
        pickle.dump(chainSeq, fOut)

# def produce_proteinSeq_dict (inPath, outPath):
#     
#     s = pd.read_table(inPath, sep="\t")
#     prSeq = {}
#     for _, row in s.iterrows():
#         prSeq[row.ID] = row.Sequence
#     with open(outPath, 'wb') as fOut:
#         pickle.dump(prSeq, fOut)

def produce_substitution_matrix (name, outPath):
    
    matrices = {'PAM30': MatrixInfo.pam30,
                'PAM60': MatrixInfo.pam60,
                'BLOSUM62' : MatrixInfo.blosum62,
                'BLOSUM95' : MatrixInfo.blosum95}
    matrix = matrices[name]
    scores = {}
    for k in matrix:
        scores[k] = matrix[k]
        scores[tuple(reversed(k))] = matrix[k]
    with open(outPath, 'wb') as fOut:
        pickle.dump(scores, fOut)

def produce_rnaToProtein_refseqID_dict (inPath, outPath):
    """Make a map of RNA RefSeq IDs to protein RefSeq IDs.

    Args:
        inPath (str): RefSeq ID mapping file.
        outPath (str): path to save output pickle file to.

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
