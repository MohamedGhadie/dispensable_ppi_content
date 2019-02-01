import os
import io
import re
import csv
import time
import pickle
from pathlib import Path
from Bio import Seq, SeqIO
import pandas as pd
import numpy as np
from simple_tools import (is_numeric,
                          listStr_to_list,
                          first_mismatch,
                          find_substring,
                          find_masked_substring)

def flanking_sequence(pos, seq, sideLength):
    
    pos = pos - 1
    seqLen = len(seq)
    if 0 <= pos < seqLen:
        s = max(0, pos - sideLength)
        e = min(seqLen, pos + sideLength + 1)
        return pos - s + 1, seq[ s : e ]
    else:
        return '-'

def base_change_pos (codonChange):

    str1, str2 = map(str.strip, codonChange.split('-'))
    return first_mismatch(str1, str2)

def parse_domine_DDI_file(inPath, outPath):
    """Read domain-domain interactions (DDIs) from DOMINE file.

    Args:
        inPath (str): DOMINE DDI file directory.
        outPath (str): file directory to save DDIs to.

    """
    interactions = pd.read_table(inPath, sep="|")
    DDIs = interactions.loc[(interactions["iPfam"]==1)|(interactions["3did"]==1), ["P1","P2"]]
    DDIs.rename(columns={"P1":"dom1","P2":"dom2"}, inplace=True)
    print('%d domain-domain interactions extracted from file ' % len(DDIs) +
          str(inPath))
    DDIs.to_csv(outPath, index=False, sep='\t')

def parse_IntAct_interactions(inPath,
                              spListFile,
                              outPath,
                              geneMapFile = None,
                              selfPPIs = True):
    """Read protein-protein interactions with valid Swiss-Prot IDs from IntAct file.

    Args:
        inPath (str): IntAct interaction file directory.
        spListFile (str): file directory with list of valid Swiss-Prot IDs
        outPath (str): file directory to save interactions to.
        selfPPIs (boolean): False to remove self-PPIs, True to include them.

    """
    with open(spListFile, 'r') as f:
        swissProtIDs = set(f.read().split())
    getGeneNames = False
    if geneMapFile is not None:
        getGeneNames = True
        with open(geneMapFile, 'rb') as f:
            geneMap = pickle.load(f)
    
    with io.open(inPath, "r", encoding="utf-8") as f, io.open(outPath, "w") as fout:
        if getGeneNames:
            geneHeaders = ["Gene_1", "Gene_2"]
        else:
            geneHeaders = []
        fout.write('\t'.join(geneHeaders + 
                             ["Protein_1",
                              "Protein_2",
                              "Protein_1_IntAct_ID",
                              "Protein_2_IntAct_ID",
                              "Interaction_ID"]) + '\n')
        i = c = 0
        next(f)
        for line in f:
            i += 1
            strsplit = line.split('\t')
            if strsplit[35].lower()=='false':
                if (len(strsplit[0])>10) and (len(strsplit[1])>10):
                    if (strsplit[0][:10].lower()=='uniprotkb:') and (strsplit[1][:10].lower()=='uniprotkb:'):
                        protein1 = strsplit[0][10:]
                        protein2 = strsplit[1][10:]
                        if (protein1 in swissProtIDs) and (protein2 in swissProtIDs):
                            if selfPPIs or (protein1 != protein2):
                                IDsplit = list(map(str.strip, strsplit[2].split('|')))
                                hasIntactID = list(map(lambda x: x.find('intact:EBI-'), IDsplit))
                                try:
                                    ind = hasIntactID.index(0)
                                    intactID1 = IDsplit[ind][7:]
                                except ValueError:
                                    intactID1 = '-'
                                IDsplit = list(map(str.strip, strsplit[3].split('|')))
                                hasIntactID = list(map(lambda x: x.find('intact:EBI-'), IDsplit))
                                try:
                                    ind = hasIntactID.index(0)
                                    intactID2 = IDsplit[ind][7:]
                                except ValueError:
                                    intactID2 = '-'
                                IDsplit = list(map(str.strip, strsplit[13].split('|')))
                                hasIntactID = list(map(lambda x: x.find('intact:EBI-'), IDsplit))
                                try:
                                    ind = hasIntactID.index(0)
                                    interactionID = IDsplit[ind][7:]
                                except ValueError:
                                    interactionID = '-'
                                if getGeneNames:
                                    gene1 = geneMap[protein1] if protein1 in geneMap else '-'
                                    gene2 = geneMap[protein2] if protein2 in geneMap else '-'
                                    genes = [gene1, gene2]
                                else:
                                    genes = []
                                c += 1
                                fout.write('\t'.join(genes +
                                                     [protein1,
                                                      protein2,
                                                      intactID1,
                                                      intactID2,
                                                      interactionID]) +  '\n')
    print('%d lines parsed from ' % i + str(inPath) + 
          ', %d PPIs extracted ' % c + 
          'and written to file ' + str(outPath))

def parse_HI_II_14_interactome(inPath,
                               UniProtIDmapFile,
                               outPath,
                               selfPPIs = True):
    """Read protein-protein interactions that map to UniProt IDs from HI-II-14 interactome file.

    Args:
        inPath (str): HI-II-14 interactome file directory.
        UniProtIDmapFile (str): file directory with UniProt ID mapping dictionary.
        outPath (str): file directory to save protein-protein interactions to.
        selfPPIs (boolean): False to remove self-PPIs, True to include them.

    """
    with open(UniProtIDmapFile, 'rb') as f:
        UniProtIDmap = pickle.load(f)
    interactome = pd.read_table(inPath, sep='\t')
    interactome.rename(columns={"Symbol A":"Gene_1",
                                "Symbol B":"Gene_2",
                                "Entrez Gene IDA":"EntrezID_1",
                                "Entrez Gene IDB":"EntrezID_2"},
                       inplace=True)
    for i, row in interactome.iterrows():
        interactome.loc[i, "Protein_1"] = (UniProtIDmap[row.Gene_1] if row.Gene_1 in UniProtIDmap
                                           else (UniProtIDmap[row.EntrezID_1] if row.EntrezID_1 in UniProtIDmap
                                                 else '-'))
        interactome.loc[i, "Protein_2"] = (UniProtIDmap[row.Gene_2] if row.Gene_2 in UniProtIDmap
                                           else (UniProtIDmap[row.EntrezID_2] if row.EntrezID_2 in UniProtIDmap
                                                 else '-'))
    interactome = interactome[(interactome[["Protein_1","Protein_2"]]!='-').all(1)].reset_index(drop=True)
    if not selfPPIs:
        interactome = interactome[interactome["Protein_1"]!=interactome["Protein_2"]].reset_index(drop=True)
    interactome.to_csv(outPath, index=False, sep='\t')

def parse_refSeq_fasta (inPath, outPath):
    
    with io.open(outPath, "w") as fout:
        fout.write('\t'.join(['ID','Length','Sequence']) + '\n')
        for entry in SeqIO.parse(str(inPath), 'fasta'):
            fout.write('\t'.join([entry.id, str(len(entry.seq)), str(entry.seq)]) +  '\n')

def merge_refSeq_sequence_files (inDir, numFiles, outPath):
    
    if not inDir.exists():
        print( '\t' + 'Directory %s does not exist' % str(inDir) )
        return
    
    seq = pd.DataFrame()
    for i in np.array(range(numFiles)) + 1:
        refseqFile = inDir / ('refseq_human_protein_' + str(i) + '.txt')
        if refseqFile.is_file():
            newSeq = pd.read_table(refseqFile, sep='\t')
            seq = seq.append(newSeq, ignore_index=True)
        else:
            print( '\t' + 'File %s not found. Skipping' % refseqFile )
    if not seq.empty:
        seq.drop_duplicates(subset="ID", inplace=True)
        seq.to_csv(outPath, index=False, sep='\t')
    else:
        print( '\t' + 'no RefSeq files found' )

def parse_fasta_file_old(inPath, outPath):
    """Read sequences from fasta file.

    Args:
        inPath (str): fasta file directory containing sequnces.
        outPath (str): path to save tab-delimited sequence file to.

    """
    s = list(SeqIO.parse(inPath, 'fasta'))
    sequences = pd.DataFrame(index=list(range(0,len(s))), columns=('ID','Length','Sequence'))
    c = -1
    for _, row in enumerate(s):
        if row.id[:3]=='sp|':
            strsplit = row.id.split('|')
            c += 1
            sequences.loc[c] = strsplit[1].strip(), len(row.seq), str(row.seq)
    sequences = sequences[:c + 1]
    sequences.to_csv(outPath, index=False, sep='\t')
    print('%d sequences extracted from file ' % (c + 1) + str(inPath) +
                ' and written to file ' + str(outPath))

def parse_fasta_file(inPath, outPath):
    """Read sequences from fasta file.

    Args:
        inPath (str): fasta file directory containing sequnces.
        outPath (str): path to save tab-delimited sequence file to.

    """
    s = list(SeqIO.parse(str(inPath), 'fasta'))
    with io.open(outPath, "w") as fout:
        fout.write('\t'.join(['ID', 'Length', 'Sequence'])  + '\n')
        for i, row in enumerate(s):
            fout.write('\t'.join([row.id, str(len(row.seq)), str(row.seq)])  + '\n')
    print('%d sequences extracted from file ' % (i + 1) + str(inPath) +
                ' and written to file ' + str(outPath))

def parse_hmmscan_file(hmmscanFile, isoformType, coordType, outPath):
    """Read protein domain mappings from HMMER hmmscan output file.

    Args:
        inPath (str): HMMER hmmscan output file directory.
        isoformType (str): type of proteins scanned in file, either 'reference', or
                            'alternative' for both reference and alternative isoforms
        coordType (str): 'alignment' or 'envelope' domain coordinates
        outPath (str): path to save protein domain mapping table to.

    """
    isoformType = isoformType.lower()
    domMap = pd.DataFrame(columns = ('spid', 'domain', 'startPos', 'endPos'))
    if not(os.path.exists(hmmscanFile)):
        print('File %s cannot be found' % hmmscanFile)
    else:
        with io.open(hmmscanFile, 'r', encoding='utf-8') as f:
            for numLines, _ in enumerate(f):
                pass
        domMap["spid"] = [np.nan] * (numLines + 1)
        c = -1
        i = -1
        with io.open(hmmscanFile, "r", encoding="utf-8") as f:
            reader = csv.reader(f)
            for row in reader:
                if row[0][0] != '#':
                    i += 1
                    print('parsed %d of %d lines' % (i,numLines))
                    strsplit = list(map(str.strip,row[0].split()))
                    dot = strsplit[1].find('.')
                    if (dot != -1):
                        if coordType == 'alignment':
                            startPos = strsplit[15]
                            endPos = strsplit[16]
                        elif coordType == 'envelope':
                            startPos = strsplit[17]
                            endPos = strsplit[18]
                        if ((strsplit[1][:2] == 'PF')
                            and strsplit[1][2:dot].isalnum()
                            and strsplit[2].isalnum()
                            and startPos.isalnum()
                            and endPos.isalnum()):
                            c += 1
                            if isoformType == 'reference':
                                domMap.loc[c] = [strsplit[3], strsplit[1][:dot],
                                                int(startPos), int(endPos)]
                            elif isoformType == 'alternative':
                                sp = strsplit[3].split('|')[1].strip()
                                domMap.loc[c] = [sp,
                                                 strsplit[1][:dot],
                                                 int(startPos),
                                                 int(endPos)]
                        else:
                            print('row %d format not readable' % i)
                            print(strsplit)
                    else:
                        print('dot not found in row %d' % i)
                        print(strsplit)
        domMap = domMap[:c + 1]
        domMap.to_csv(outPath, index=False, sep='\t')
        print('%d lines parsed from ' % (i + 1) + str(hmmscanFile) +
              ', %d mappings extracted' % len(domMap) + 
              ' and written to file ' + str(outPath))

def parse_DDI_files(_3didFile, domineFile, _3didOutPath, domineOutPath, outPath):
    """Read domain-domain interactions (DDIs) from DOMINE and 3did files.

    Args:
        _3didFile (str): 3did DDI file directory.
        domineFile (str): DOMINE DDI file directory.
        _3didOutPath (str): file directory to save 3did DDIs to.
        domineOutPath (str): file directory to save DOMINE DDIs to.
        outPath (str): file directory to save all DDIs to, with no repetition.

    """
    parse_3did_DDI_file(_3didFile, _3didOutPath)
    parse_domine_DDI_file(domineFile, domineOutPath)
    _3didDDIs = pd.read_table(_3didOutPath, sep='\t')
    _3didDDIs = _3didDDIs[["dom1","dom2"]]
    domineDDIs = pd.read_table(domineOutPath, sep='\t')
    DDIs = _3didDDIs.append(domineDDIs, ignore_index=True)
    DDIs = DDIs.drop_duplicates(subset=["dom1","dom2"]).reset_index(drop=True)
    # drop duplicates that are in reverse order
    DDIs["duplicate"] = False
    for i, row in DDIs.iterrows():
        isduplicate = list((DDIs["dom1"]==row.dom2) & (DDIs["dom2"]==row.dom1))
        if sum(isduplicate) == 1:
            DDIs.loc[isduplicate,"duplicate"] = (isduplicate.index(True)>i)
    DDIs = DDIs.loc[DDIs['duplicate']==False,["dom1","dom2"]]
    DDIs.to_csv(outPath, index=False, sep='\t')
    print('%d domain-domain interactions written to file ' % len(DDIs) +
          str(outPath))
    
def parse_3did_DDI_file(inPath, outPath):
    """Read domain-domain interactions (DDIs) from 3did file.

    Args:
        inPath (str): 3did DDI file directory.
        outPath (str): file directory to save DDIs to.

    """
    numLines = 10000
    c = -1
    DDIs = pd.DataFrame('NA', index=range(numLines),
                        columns=["dom1", "dom2", "dom1_name", "dom2_name"])
    with io.open(inPath, "r", encoding="utf-8") as f:
        reader = csv.reader(f)
        for row in reader:
            if len(row[0]) > 4:
                if row[0][:4] == '#=ID':
                    ddi = row[0].split('\t')
                    if len(ddi) == 5:
                        if ((ddi[3].find('@Pfam') > -1) and
                                (ddi[4].find('@Pfam') > -1)):
                            s1 = ddi[3].find('PF')
                            e1 = ddi[3].find('.')
                            s2 = ddi[4].find('PF')
                            e2 = ddi[4].find('.')
                            if (-1 < s1 < e1) and (-1 < s2 < e2):
                                c += 1
                                DDIs.loc[c, ["dom1", "dom2", "dom1_name", "dom2_name"]] = ddi[3][s1:e1], ddi[4][s2:e2], ddi[1].strip(), ddi[2].strip()
    print('%d domain-domain interactions extracted from file ' % (c + 1) +
          str(inPath))
    DDIs.to_csv(outPath, index=False, sep='\t')

def parse_3did_DDIinterface_file(inPath, outPath):
    """Read domain-domain interactions (DDIs) with interface residue interactions

    Args:
        inPath (str): 3did interface DDI file directory.
        outPath (str): file directory to save DDIs to.

    """
    numDDIs = 0
    with io.open(inPath, 'r', encoding='utf-8') as f:
        reader = csv.reader(f)
        for row in reader:
            if len(row[0]) > 4:
                if row[0][:4] == '#=ID':
                    numDDIs += 1
    c = -1
    DDIs = pd.DataFrame('NA', index=range(numDDIs),
                        columns=["dom1_name", "dom2_name","interface1_pos","interface2_pos"])
    with io.open(inPath, "r", encoding="utf-8") as f:
        reader = csv.reader(f)
        for row in reader:
            if len(row[0]) > 4:
                if row[0][:4] == '#=ID':
                    ddi = list(map(str.strip, row[0].split('\t')))
                    if len(ddi) == 3:
                        c += 1
                        print(c)
                        DDIs.loc[c, ["dom1_name", "dom2_name"]] = ddi[1], ddi[2]
                        interface1 = []
                        interface2 = []
                elif row[0][:4] == '#=IF':
                    resIx = list(map(str.strip, row[0].split('\t')))
                    if len(resIx) == 2:
                        resIx = list(map(str.strip, resIx[1].split(' ')))
                        resIx = [x.split('-') for x in resIx[1:]]
                        interface1.extend(map(int, [x for [x,_] in resIx]))
                        interface2.extend(map(int, [x for [_,x] in resIx]))
            elif row[0][:2] == '//':
                DDIs.loc[c, ["interface1_pos", "interface2_pos"]] = [sorted(list(set(interface1))),
                                                                           sorted(list(set(interface2)))]
                
    print('%d domain-domain interactions extracted from file ' % (c + 1) +
          str(inPath))
    DDIs.to_csv(outPath, index=False, sep='\t')

def read_DDIinterface_file(inPath):
    """Read DDIs with their interaction interfaces.

    Args:
        inPath (str): file directory containing DDIs with interface residue positions.
    
    Returns:
        DataFrame: DDIs with interface residue positions in list form.

    """
    DDIs = pd.read_table(inPath, sep='\t')
    DDIs["interface1_pos"] = DDIs["interface1_pos"].apply(lambda x: x[1:-1])
    DDIs["interface1_pos"] = DDIs["interface1_pos"].apply(listStr_to_list)
    DDIs["interface1_pos"] = DDIs["interface1_pos"].apply(lambda x: list(map(int,x)))
    DDIs["interface2_pos"] = DDIs["interface2_pos"].apply(lambda x: x[1:-1])
    DDIs["interface2_pos"] = DDIs["interface2_pos"].apply(listStr_to_list)
    DDIs["interface2_pos"] = DDIs["interface2_pos"].apply(lambda x: list(map(int,x)))
    return DDIs

def merge_duplicate_DDI_interfaces(inPath, outPath):
    """Merge interfaces of duplicate DDIs.

    Args:
        inPath (str): file directory containing DDIs with interface residue positions.
        outPath (str): file directory to save unique DDIs with interface residue positions.

    """
    DDIs = read_DDIinterface_file(inPath)
    mergedDDIs = pd.DataFrame('NA', index=list(range(len(DDIs))),
                                columns=("dom1",
                                         "dom2",
                                         "dom1_name",
                                         "dom2_name",
                                         "interface1_pos",
                                         "interface2_pos"))
    c = -1
    checked = pd.Series([False]*len(DDIs))
    for i, ddi in DDIs.iterrows():
        if not checked[i]:
            duplicates = (DDIs["dom1"]==ddi.dom1) & (DDIs["dom2"]==ddi.dom2)
            interface1_array = DDIs.loc[duplicates,"interface1_pos"].values
            interface1 = [pos for interface in interface1_array for pos in interface]
            interface2_array = DDIs.loc[duplicates,"interface2_pos"].values
            interface2 = [pos for interface in interface2_array for pos in interface]
            checked[duplicates] = True
            
            duplicates = (DDIs["dom2"]==ddi.dom1) & (DDIs["dom1"]==ddi.dom2)
            if duplicates.sum() > 0:
                interface1_array = DDIs.loc[duplicates,"interface2_pos"].values
                interface1b = [pos for interface in interface1_array for pos in interface]
                interface2_array = DDIs.loc[duplicates,"interface1_pos"].values
                interface2b = [pos for interface in interface2_array for pos in interface]
                checked[duplicates] = True
                interface1.extend(interface1b)
                interface2.extend(interface2b)
            c += 1
            mergedDDIs.loc[c, ["dom1",
                               "dom2",
                               "dom1_name",
                               "dom2_name"]] = (ddi.dom1,
                                                ddi.dom2,
                                                ddi.dom1_name,
                                                ddi.dom2_name)
            mergedDDIs.loc[c, ["interface1_pos",
                               "interface2_pos"]] = (sorted(list(set(interface1))),
                                                     sorted(list(set(interface2))))
    mergedDDIs = mergedDDIs[:c + 1]
    mergedDDIs.to_csv(outPath, index=False, sep='\t')

def parse_HGMD_disease_mutations (inPath, uniprotIDmapFile, outPath):
    """Parse HGMD file for missense disease substitution mutations whose RefSeq ID maps
        to a UniProt ID

    Args:
        inPath (str): file directory containing HGMD mutations.
        uniprotIDmapFile (str): file directory containing UniProt mapping table.
        outPath (str): file directory to save mutations to.
    
    """    
    mutations = pd.read_table(inPath, sep='\t', encoding='iso8859_11')
    mutations.rename(columns={"gene":"Gene",
                              "sequence_context_hg19":"Context",
                              "codon_change":"Codon_change",
                              "codon_number":"Codon_number",
                              "disease":"Disease"},
                              inplace=True)
    with open(uniprotIDmapFile, 'rb') as f:
        uniprotIDmap = pickle.load(f)
    mutations["Protein"] = mutations["Gene"].apply(lambda x: uniprotIDmap[x] if x in uniprotIDmap else '-')
    mutations = mutations[(mutations["Type"]==1)
                          & mutations["Context"].apply(lambda x: isinstance(x, str))
                          & (mutations["Protein"] != '-')
                          & (mutations["Variant_class"] == 'DM')].reset_index(drop=True)
    mutations.to_csv(outPath, index=False, sep='\t')

def translate_HGMD_flanking_sequences (inPath, outPath):
  
    mutations = pd.read_table(inPath, sep='\t')
    mutations["Base_change_pos"] = mutations["Codon_change"].apply(base_change_pos)
    
    mutations["WT_context"] = mutations.apply(lambda x:
                                               translate_HGMD_wt_flanking_sequence(x["Context"],
                                                                                   x["Base_change_pos"]),
                                               axis=1)
    mutations["Mut_context"] = mutations.apply(lambda x:
                                               translate_HGMD_mut_flanking_sequence(x["Context"],
                                                                                    x["Base_change_pos"]),
                                               axis=1)
    mutations["Mut_context_pos"] = mutations["WT_context"].apply(lambda x: x[0])
    mutations["WT_context"] = mutations["WT_context"].apply(lambda x: x[1])
    mutations["Mut_context"] = mutations["Mut_context"].apply(lambda x: x[1])
    mutations.to_csv(outPath, index=False, sep='\t')

def truncate_HGMD_flanking_sequences (inPath, outPath):

    mutations = pd.read_table(inPath, sep='\t')
    mutations["WT_context"] = mutations.apply(lambda x:
                                               truncate_sequence_left_side(x["Codon_number"],
                                                                           x["Mut_context_pos"],
                                                                           x["WT_context"]),
                                               axis=1)
    mutations["Mut_context"] = mutations.apply(lambda x:
                                               truncate_sequence_left_side(x["Codon_number"],
                                                                           x["Mut_context_pos"],
                                                                           x["Mut_context"]),
                                               axis=1)
    mutations["Mut_context_pos"] = mutations["WT_context"].apply(lambda x: x[0])
    mutations["WT_context"] = mutations["WT_context"].apply(lambda x: x[1])
    mutations["Mut_context"] = mutations["Mut_context"].apply(lambda x: x[1])
    
    mutations["WT_context"] = mutations.apply(lambda x:
                                              truncate_sequence_right_side(x["Mut_context_pos"],
                                                                           x["WT_context"]),
                                              axis=1)
    mutations["Mut_context"] = mutations.apply(lambda x:
                                               truncate_sequence_right_side(x["Mut_context_pos"],
                                                                            x["Mut_context"]),
                                               axis=1)
    mutations.to_csv(outPath, index=False, sep='\t')

def match_flanking_sequences (inPath, sequenceFile, outPath):

    mutations = pd.read_table(inPath, sep='\t')
    sequences = pd.read_table(sequenceFile, sep='\t')
    sequencedProteins = sequences["ID"].values
    mutations["Protein_seq"] = mutations["Protein"].apply(lambda x:
                                                          sequences.loc[sequences["ID"]==x,
                                                                        "Sequence"].item() 
                                                          if x in sequencedProteins
                                                          else '-')
    mutations["Seq_match"] = mutations.apply(lambda x:
                                             find_substring(x["WT_context"],
                                                            x["Protein_seq"]),
                                             axis=1)
    mutations["Seq_match"] = mutations.apply(lambda x:
                                             (x["Codon_number"]-x["Mut_context_pos"])
                                             in x["Seq_match"],
                                             axis=1)
    mutations.drop("Protein_seq", axis=1, inplace=True)
    mutations.to_csv(outPath, index=False, sep='\t')

def remove_hgmd_synon_nonsense_mutations(inPath, outPath):
    
    mutations = pd.read_table(inPath, sep='\t')
    nonsynonymous = mutations.apply(lambda x:
                                    (x["Mut_context"][x["Mut_context_pos"] - 1]
                                     != x["WT_context"][x["Mut_context_pos"] - 1]),
                                    axis=1)
    missense = mutations.apply(lambda x:
                               x["Mut_context"][x["Mut_context_pos"] - 1] != '*',
                               axis=1)
    mutations = mutations[nonsynonymous & missense]
    mutations.to_csv(outPath, index=False, sep='\t')

def translate_HGMD_mut_flanking_sequence (dnaSeq, baseChangePos):

    firstBracket = dnaSeq.find('[')
    mutCodonStart = firstBracket-baseChangePos
    orfStart = mutCodonStart % 3
    mutCodonPos = (mutCodonStart - orfStart) // 3
    mutSeq = dnaSeq[orfStart:].replace(']','')
    mutSeq = mutSeq[:mutSeq.find('[')] + mutSeq[mutSeq.find('/')+1:]
    seqLen = len(mutSeq)
    mutSeq = mutSeq[:seqLen-(seqLen%3)]
    return mutCodonPos + 1, Seq.translate(mutSeq,
                                          table = 'Standard',
                                          stop_symbol = '*',
                                          to_stop = False,
                                          cds = False,
                                          gap = None)

def translate_HGMD_wt_flanking_sequence (dnaSeq, baseChangePos):

    firstBracket = dnaSeq.find('[')
    mutCodonStart = firstBracket-baseChangePos
    orfStart = mutCodonStart % 3
    mutCodonPos = (mutCodonStart - orfStart) // 3
    refSeq = dnaSeq[orfStart:].replace('[','')
    refSeq = refSeq[:refSeq.find('/')] + refSeq[refSeq.find(']')+1:]
    seqLen = len(refSeq)
    refSeq = refSeq[:seqLen-(seqLen%3)]
    return mutCodonPos + 1, Seq.translate(refSeq,
                                          table = 'Standard',
                                          stop_symbol = '*',
                                          to_stop = False,
                                          cds = False,
                                          gap = None)

def truncate_sequence_left_side (mutGlobalPos, mutPos, seq):
    
    if mutGlobalPos < mutPos:
        seq = seq[ mutPos - mutGlobalPos :]
        mutPos = mutGlobalPos
    return mutPos, seq

def truncate_sequence_right_side (mutPos, seq):
    
    stopPos = [m.start() for m in re.finditer('\\*', seq) if m.start() > (mutPos - 1)]
    if  len(stopPos) > 0:
        seq = seq[ : stopPos[0] ]
    return seq

def parse_genome_project_mutations(inPath, outPath):

    with io.open(inPath, "r") as f:
        for headerLines, line in enumerate(f):
            if line[0] != '#':
                break
    
    infoHeaders = set()
    with io.open(inPath, "r") as fin:
        for k in range(headerLines):
            next(fin)
        for numLines, line in enumerate(fin):
            linesplit = list(map(str.strip, line.split('\t')))
            info = list(map(str.strip, linesplit[8].split(';')))
            infosplit = [list(map(str.strip, x.split('='))) for x in info]
            for header, _ in infosplit:
                infoHeaders.add(header)
    
    numLines += 1
    infoHeaders = sorted(infoHeaders)
    allHeaders = ['Chromosome',
                  'Db',
                  'Effect',
                  'Chr_start',
                  'Chr_end',
                  'dot1',
                  'dot2',
                  'dot3'] + infoHeaders
    mut = {k:'-' for k in infoHeaders}
    
    with io.open(inPath, "r") as fin, io.open(outPath, "a") as fout:
        fout.write('\t'.join(allHeaders) + '\n')
        for k in range(headerLines):
            next(fin)
        for i, line in enumerate(fin):
            if (i % 100000) == 0:
                print('\t- %d out of %d lines parsed' % (i, numLines))
            linesplit = list(map(str.strip, line.split('\t')))
            allvalues = linesplit[:8]
            info = list(map(str.strip, linesplit[8].split(';')))
            infosplit = [list(map(str.strip, x.split('='))) for x in info]
            mut.update({k:'-' for k in mut})
            mut.update({k:val for k, val in infosplit})
            allvalues.extend([mut[k] for k in infoHeaders]) 
            fout.write('\t'.join(allvalues) + '\n')
    print('- %d mutations parsed from ' % (i + 1) + str(inPath) +  
          'and written to file ' + str(outPath))

def parse_dbsnp_flatfile_keys (inPath, encod, pausetime, outDir):
    
    RSkeys = set()
    SNPkeys = set()
    VALkeys = set()
    MAFkeys = set()
    CLINkeys = set()
    LOCkeys = set()
    
    with io.open(inPath, 'r', encoding=encod) as fin:
        for numLines, line in enumerate(fin):
            if not (numLines % 1e6):
                print('\t\t%dM lines parsed' % (numLines // 1e6))
            if not ((numLines+1) % 50e6):
                print('\t\tpausing for %d seconds' % pausetime)
                time.sleep(pausetime)
            if line.startswith(('rs','SNP','VAL','GMAF','CLINSIG','LOC')):
                linesplit = list(map(str.strip, line.split('|')))
                entrysplit = [list(map(str.strip, x.split('='))) for x in linesplit]
                entrysplit = [entry for entry in entrysplit if len(entry)==2]
                if linesplit[0][:2] == 'rs':
                    for k, _ in entrysplit:
                        RSkeys.add(k)
                elif linesplit[0] == 'SNP':
                    for k, _ in entrysplit:
                        SNPkeys.add(k)
                elif linesplit[0] == 'VAL':
                    for k, _ in entrysplit:
                        VALkeys.add(k)
                elif linesplit[0] == 'GMAF':
                    for k, _ in entrysplit:
                        MAFkeys.add(k)
                elif linesplit[0] == 'CLINSIG':
                    for k, _ in entrysplit:
                        CLINkeys.add(k)
                elif linesplit[0] == 'LOC':
                    for k, _ in entrysplit:
                        LOCkeys.add(k)
    with open(outDir / 'dbsnp_RSkeys.pkl', 'wb') as fout:
        pickle.dump(RSkeys, fout)
    with open(outDir / 'dbsnp_SNPkeys.pkl', 'wb') as fout:
        pickle.dump(SNPkeys, fout)
    with open(outDir / 'dbsnp_VALkeys.pkl', 'wb') as fout:
        pickle.dump(VALkeys, fout)
    with open(outDir / 'dbsnp_MAFkeys.pkl', 'wb') as fout:
        pickle.dump(MAFkeys, fout)
    with open(outDir / 'dbsnp_CLINSIGkeys.pkl', 'wb') as fout:
        pickle.dump(CLINkeys, fout)
    with open(outDir / 'dbsnp_LOCkeys.pkl', 'wb') as fout:
        pickle.dump(LOCkeys, fout)
        
def parse_dbsnp_flatfile (inPath, keyDir, encod, pausetime, outPath):
    
    with open(keyDir / 'dbsnp_RSkeys.pkl', 'rb') as f:
        RSkeys = pickle.load(f)
    with open(keyDir / 'dbsnp_SNPkeys.pkl', 'rb') as f:
        SNPkeys = pickle.load(f)
    with open(keyDir / 'dbsnp_VALkeys.pkl', 'rb') as f:
        VALkeys = pickle.load(f)
    with open(keyDir / 'dbsnp_MAFkeys.pkl', 'rb') as f:
        MAFkeys = pickle.load(f)
    with open(keyDir / 'dbsnp_CLINSIGkeys.pkl', 'rb') as f:
        CLINkeys = pickle.load(f)
    with open(keyDir / 'dbsnp_LOCkeys.pkl', 'rb') as f:
        LOCkeys = pickle.load(f)
        
    RSkeys = sorted(RSkeys)
    SNPkeys = sorted(SNPkeys)
    VALkeys = sorted(VALkeys)
    MAFkeys = sorted(MAFkeys)
    CLINkeys = sorted(CLINkeys)
    LOCkeys = sorted(LOCkeys)
    
    RSkeys = ['ID',
              'species',
              'class',
              '1000genome',
              'date'] + RSkeys
    VALkeys = VALkeys + ['validation',
                         'method']
    LOCkeys = ['gene'] + LOCkeys
    
    RSinfo = {k:k for k in RSkeys}
    SNPinfo = {k:k for k in SNPkeys}
    VALinfo = {k:k for k in VALkeys}
    MAFinfo = {k:k for k in MAFkeys}
    CLINinfo = {k:k for k in CLINkeys}
    LOCinfo = {k:k for k in LOCkeys}
    
    with io.open(inPath, "r", encoding=encod) as fin:
        for numLines, line in enumerate(fin):
            pass
                
    with io.open(inPath, "r", encoding=encod) as fin, io.open(outPath, "a") as fout:
        for i, line in enumerate(fin):
            if not (i % 1e6):
                print('\t\t%dM out of %gM lines parsed' % (i // 1e6, numLines / 1e6))
            if not ((i+1) % 50e6):
                print('\t\tpausing for %d seconds' % pausetime)
                time.sleep(pausetime)
            if line.startswith(('rs','SNP','VAL','GMAF','CLINSIG','LOC')) or (line.startswith('ss')
                                                                              and (RSinfo["1000genome"] != 'yes')):
                linesplit = list(map(str.strip, line.split('|')))
                entrysplit = [list(map(str.strip, x.split('='))) for x in linesplit]
                if linesplit[0][:2] == 'rs':
                    # save previous mutation
                    allvalues = [RSinfo[k] for k in RSkeys]
                    allvalues.extend([SNPinfo[k] for k in SNPkeys])
                    allvalues.extend([VALinfo[k] for k in VALkeys])
                    allvalues.extend([MAFinfo[k] for k in MAFkeys])
                    allvalues.extend([CLINinfo[k] for k in CLINkeys])
                    allvalues.extend([LOCinfo[k] for k in LOCkeys])
                    fout.write('\t'.join(allvalues) + '\n')
                
                    # start with next mutation
                    RSinfo.update({k:'?' for k in RSkeys})
                    SNPinfo.update({k:'?' for k in SNPkeys})
                    VALinfo.update({k:'?' for k in VALkeys})
                    MAFinfo.update({k:'?' for k in MAFkeys})
                    CLINinfo.update({k:'?' for k in CLINkeys})
                    LOCinfo.update({k:'?' for k in LOCkeys})
                    RSinfo["1000genome"] = 'no'
                    RSinfo["ID"] = linesplit[0]
                    RSinfo["species"] = linesplit[1]
                    RSinfo["class"] = linesplit[3]
                    RSinfo["date"] = linesplit[-1]
                    RSinfo.update({k[0]:k[1] for k in entrysplit if len(k)==2})
                elif linesplit[0][:2] == 'ss':
                    RSinfo["1000genome"] = ('yes' if linesplit[1] == '1000GENOMES' else 'no')
                elif linesplit[0] == 'SNP':
                    SNPinfo.update({k[0]:k[1] for k in entrysplit if len(k)==2})
                elif linesplit[0] == 'VAL':
                    if (len(entrysplit[-1]) == 1) and (len(entrysplit[-2]) == 1):
                        VALinfo["method"] = linesplit[-1]
                        VALinfo["validation"] = linesplit[-2]
                    elif (len(entrysplit[-1]) == 1):
                        VALinfo["validation"] = linesplit[-1]
                    VALinfo.update({k[0]:k[1] for k in entrysplit if len(k)==2})
                elif linesplit[0] == 'GMAF':
                    MAFinfo.update({k[0]:k[1] for k in entrysplit if len(k)==2})
                elif linesplit[0] == 'CLINSIG':
                    CLINinfo.update({k[0]:k[1] for k in entrysplit if len(k)==2})
                elif linesplit[0] == 'LOC':
                    if 'fxn-class=missense' in line:
                        if ((LOCinfo["fxn-class"] != 'missense') or
                            (LOCinfo["aa_position"] == '?') or
                            (LOCinfo["prot_acc"] == '?')):
                            LOCinfo.update({k:'?' for k in LOCkeys})
                            LOCinfo["gene"] = linesplit[1]
                            LOCinfo.update({k[0]:k[1] for k in entrysplit if len(k)==2})           
        allvalues = [RSinfo[k] for k in RSkeys]
        allvalues.extend([SNPinfo[k] for k in SNPkeys])
        allvalues.extend([VALinfo[k] for k in VALkeys])
        allvalues.extend([MAFinfo[k] for k in MAFkeys])
        allvalues.extend([CLINinfo[k] for k in CLINkeys])
        allvalues.extend([LOCinfo[k] for k in LOCkeys])
        fout.write('\t'.join(allvalues) + '\n')
    print('\t%d lines parsed from ' % (i + 1) + str(inPath) +  
          'and written to file ' + str(outPath))

def filter_and_merge_dbsnp_mutations (inDir, uniprotIDmapFile, pausetime, outPath):
    
    with open(uniprotIDmapFile, 'rb') as f:
        uniprotIDmap = pickle.load(f)
    chromosomes = list(map(str, list(np.array(range(22)) + 1) + ['X', 'Y']))
    for i in chromosomes:
        dbsnpChrMutationsFile = inDir / ('dbsnp_chr' + i + '.txt')
        if dbsnpChrMutationsFile.is_file():
            print('\tadding mutations from chromosome ' + i)
            newMut = pd.read_table(dbsnpChrMutationsFile, dtype='str', sep='\t')
            print('\t\tdata read')
            newMut["assertion"] = newMut["assertion"].apply(str.lower)
            newMut["validated"] = newMut["validated"].apply(str.lower)
            newMut["validation"] = newMut["validation"].apply(str.lower)
            print('\t\tfirst conversion done')
            newMut["MAF"] = newMut["MAF"].apply(lambda x:
                                                float(x) if is_numeric(x) else -1)
            newMut["aa_position"] = newMut["aa_position"].apply(lambda x:
                                                                int(x) if is_numeric(x) else -1)
            print('\t\tsecond conversion done')
        
            common = (newMut["MAF"] >= 0.01)
            human = (newMut["species"] == 'human')
            snp = (newMut["class"] == 'snp')
            missense = (newMut["fxn-class"] == 'missense')
            validated = (newMut["validated"] == 'yes')
            notwithdrawn = (newMut["validation"] == 'notwithdrawn')
            positioned = (newMut["aa_position"] > -1)
            sequenced = (newMut["prot_acc"] != '?')
            nonpathogenic = newMut["assertion"].apply(lambda x:
                                                     x not in {'pathogenic',
                                                               'likely pathogenic',
                                                               'uncertain significance',
                                                               'other'})
            print('Common: %d' % sum(common))
            print('Human: %d' % sum(human))
            print('SNP: %d' % sum(snp))
            print('Missense: %d' % sum(missense))
            print('Validated: %d' % sum(validated))
            print('Not withdrawn: %d' % sum(notwithdrawn))
            print('Positioned: %d' % sum(positioned))
            print('Sequenced: %d' % sum(sequenced))
            print('Non-pathogenic: %d' % sum(nonpathogenic))
            print('\t\tfirst selection done')
            selected = (human &
                        snp &
                        missense &
                        nonpathogenic &
                        validated &
                        notwithdrawn &
                        common &
                        positioned &
                        sequenced)
            print('\t\t%d out of %d mutations selected from chromosome ' 
                  % (sum(selected), len(newMut)) + i)
            newMut = newMut[selected]
            newMut["gene"] = newMut["gene"].apply(str.upper)
            newMut["protein"] = newMut["gene"].apply(lambda x: uniprotIDmap[x] if x in uniprotIDmap else '-')
            newMut = newMut[newMut["protein"] != '-']
            newMut.drop_duplicates(subset=["gene", "aa_position"], inplace=True)
            newMut.to_csv(outPath, header=(i=='1'), mode = 'a', index=False, sep='\t')
            print('\t\tpausing for %d seconds' % pausetime)
            time.sleep(pausetime)
    allmut = pd.read_table(outPath, dtype='str', sep='\t')
    allmut.drop_duplicates(subset=["gene", "aa_position"], inplace=True)
    allmut.to_csv(outPath, index=False, sep='\t')

def get_flanking_sequences (inPath, sequenceFile, sideLength, outPath):
    
    mutations = pd.read_table(inPath, sep='\t')
    sequences = pd.read_table(sequenceFile, sep='\t')
    sequencedProteins = sequences["ID"].values
    mutations["protein_seq"] = mutations["prot_acc"].apply(lambda x:
                                                          sequences.loc[sequences["ID"]==x, "Sequence"].item() 
                                                          if x in sequencedProteins
                                                          else '-')
    mutations = mutations[mutations["protein_seq"]!='-'].reset_index(drop=True)
    mutations["context"] = mutations.apply(lambda x:
                                           flanking_sequence(x["aa_position"],
                                                            x["protein_seq"],
                                                            sideLength),
                                           axis=1)
    mutations = mutations[mutations["context"]!='-'].reset_index(drop=True)
    mutations["mut_context_pos"] = mutations["context"].apply(lambda x: x[0])
    mutations["context"] = mutations["context"].apply(lambda x: x[1])
    mutations.drop("protein_seq", axis=1, inplace=True)
    mutations.to_csv(outPath, index=False, sep='\t')

def match_masked_flanking_sequences (inPath, sequenceFile, outPath):

    mutations = pd.read_table(inPath, sep='\t')
    sequences = pd.read_table(sequenceFile, sep='\t')
    sequencedProteins = sequences["ID"].values
    mutations["protein_seq"] = mutations["protein"].apply(lambda x:
                                                          sequences.loc[sequences["ID"]==x, "Sequence"].item() 
                                                          if x in sequencedProteins
                                                          else '-')
    mutations["seq_match"] = mutations.apply(lambda x:
                                             find_masked_substring(x["context"],
                                                                   x["protein_seq"],
                                                                   x["mut_context_pos"]-1),
                                             axis=1)
    mutations["seq_match"] = mutations.apply(lambda x:
                                             (x["aa_position"]-x["mut_context_pos"])
                                             in x["seq_match"],
                                             axis=1)
    mutations.drop("protein_seq", axis=1, inplace=True)
    mutations.to_csv(outPath, index=False, sep='\t')

def remove_dbsnp_synon_nonsense_mutations(inPath, outPath):
    
    mutations = pd.read_table(inPath, sep='\t')
    nonsynonymous = mutations.apply(lambda x: (x["context"][x["mut_context_pos"] - 1]
                                               != x["residue"]), 
                                               axis=1)
    missense = mutations.apply(lambda x: x["residue"] != '*', axis=1)
    mutations = mutations[nonsynonymous & missense]
    mutations.to_csv(outPath, index=False, sep='\t')

def parse_uniprot_missense_mutations (inPath, outPath):
    """Parse UniProt mutation file for missense disease and natural mutations

    Args:
        inPath (str): file directory containing UniProt mutations.
        outPath (str): file directory to save mutations to.
    
    """
    headerLines = 49
    with io.open(inPath, 'r', encoding='utf-8') as f:
        for numLines, _ in enumerate(f):
            pass
    numLines = numLines - headerLines
    mutations = pd.DataFrame(index=list(range(0,numLines+1)),
                                columns=("Gene",
                                         "Protein",
                                         "Mutation_Position",
                                         "Phenotype",
                                         "Disease"),
                                dtype='str')
    
    c = -1
    i = -1
    with io.open(inPath, "r", encoding="utf-8") as f:
        for k in range(headerLines):
            next(f)
        for line in f:
            i += 1
            if (i % 1000) == 0:
                print('%d out of %d lines parsed' % (i, numLines))
                print('%d mutations extracted' % (c + 1))
            strsplit = list(map(str.strip, line.split()))
            mutPos = strsplit[3]
            mutPos = int(mutPos[5:(len(mutPos)-3)])
            disease = ' '.join(strsplit[6:])
            c += 1
            mutations.loc[c] = [strsplit[0],
                                strsplit[1],
                                mutPos,
                                strsplit[4],
                                disease]
    mutations = mutations[:c + 1]
    mutations.to_csv(outPath, index=False, sep='\t')
    print('%d lines parsed from ' % (i + 1) + str(inPath) + 
          ', %d mutations extracted ' % (c + 1) + 
          'and written to file ' + str(outPath))

def parse_blast_file (inPath, encod, pausetime, outPath):
    
    allignKeys = ['Query',
                  'Qlen',
                  'Subject',
                  'Slen',
                  'Score',
                  'Expect',
                  'Identities',
                  'Positives',
                  'Gaps',
                  'Qstart',
                  'Qend',
                  'Sstart',
                  'Send',
                  'Qseq',
                  'Sseq',
                  'Match']
    allignInfo = {k:k for k in allignKeys}
    currentInfo = {k:'-' for k in allignKeys[:4]}
    previous = '-'
    
    with io.open(inPath, "r", encoding=encod) as f:
        for numLines, _ in enumerate(f):
            pass
            
    with io.open(inPath, "r", encoding=encod) as fin, io.open(outPath, "w") as fout:
        for i, line in enumerate(fin):
            if not (i % 1e6):
                print('\t\t%d million out of %g million lines parsed' % (i // 1e6, numLines / 1e6))
            if not ((i+1) % 50e6):
                print('\t\tpausing for %d seconds' % pausetime)
                time.sleep(pausetime)
            line = line.replace('\n','')
            stripline = line.strip()
            if stripline.startswith('Query='):
                linesplit = list(map(str.strip, stripline.split('=')))
                currentInfo["Query"] = linesplit[1]
                previous = 'query'
            elif stripline.startswith('>'):
                currentInfo["Subject"] = stripline[1:].strip()
                previous = 'subject'
            elif stripline.startswith('Length'):
                linesplit = list(map(str.strip, stripline.split('=')))
                if previous=='query':
                    currentInfo["Qlen"] = linesplit[1]
                elif previous=='subject':
                    currentInfo["Slen"] = linesplit[1]
            elif any(x in stripline for x in ['Score =', 'Except =', 'Identities =', 'Positives =', 'Gaps =']):
                linesplit = list(map(str.strip, stripline.split(',')))
                entrysplit = [list(map(str.strip, x.split('='))) for x in linesplit]
                for entry in entrysplit:
                    if len(entry) == 2:
                        if entry[0] == 'Score':
                            allignInfo.update({k:allignInfo[k].replace('\n','')
                                              for k in ['Qseq','Sseq','Match']})
                            info = [allignInfo[k] for k in allignKeys]
                            fout.write('\t'.join(info) + '\n')
                            allignInfo.update({k:'-' for k in allignKeys})
                            allignInfo.update({k:'' for k in ['Qseq','Sseq','Match']})
                            allignInfo.update({k:currentInfo[k] for k in currentInfo})
                            allignInfo["Score"] = entry[1].split()[0]
                        elif entry[0] == 'Expect':
                            allignInfo["Expect"] = entry[1]
                        elif entry[0] in ['Identities', 'Positives', 'Gaps']:
                            allignInfo[entry[0]] = entry[1].split('/')[0].strip()
            elif stripline.startswith('Query'):
                linesplit = stripline.split()
                if allignInfo["Qstart"] == '-':
                    allignInfo["Qstart"] = linesplit[1]
                allignInfo["Qend"] = linesplit[3]
                allignInfo["Qseq"] = allignInfo["Qseq"] + linesplit[2]
                allignLen = len(linesplit[2])
                previous = 'query seq'
            elif stripline.startswith('Sbjct'):
                linesplit = stripline.split()
                if allignInfo["Sstart"] == '-':
                    allignInfo["Sstart"] = linesplit[1]
                allignInfo["Send"] = linesplit[3]
                allignInfo["Sseq"] = allignInfo["Sseq"] + linesplit[2]
                previous = 'subject seq'
            elif (len(line) > 0) and (previous == 'query seq'):
                allignInfo["Match"] = allignInfo["Match"] + line[ - allignLen : ]
        info = [allignInfo[k] for k in allignKeys]
        fout.write('\t'.join(info) + '\n')
    print('\t%d lines parsed from ' % (i + 1) + str(inPath) +  
          'and written to file ' + str(outPath))


def reduce_fasta_headers (inPath, delimiter, first, last, outPath):
    
    s = list(SeqIO.parse(str(inPath), 'fasta'))
    with io.open(outPath, "w") as fout:
        for _, row in enumerate(s):
            idsplit = list(map(str.strip, row.id.split(delimiter)))
            id = '|'.join(idsplit[ first - 1 : last ])
            fout.write('>' + id + '\n')
            fout.write(str(row.seq) + '\n')

def write_fasta_file (data, idCol, seqCol, outPath):
    
    with io.open(outPath, "w") as fout:
        for _, row in data.iterrows():
            fout.write('>' + row[idCol] + '\n')
            fout.write(row[seqCol] + '\n')

def produce_item_list(inPath, col, outPath):
    """Write unique items from a dataframe column to file

    Args:
        inPath (str): file directory containing dataframe
        col (str): name of column to save
        outPath (str): file directory to save list of items to

    """
    df = pd.read_table(inPath, sep='\t')
    items = sorted(set(df[col]))
    with open(outPath, 'w') as fout:
        for item in items:
            fout.write('%s\n' % item)

def read_list_table (inPath, cols, dtyp, delm):
    
    df = pd.read_table(inPath, delm)
    if not isinstance(cols, (list, tuple)):
        cols = [cols]
        dtyp = [dtyp]
    for i, col in enumerate(cols):
        df[col] = df[col].apply(lambda x: list(map(dtyp[i], map(str.strip, x.split(',')))))
    return df

def read_bindprofx_mutations (inPath):
    
    mutations = {}
    ppis = set()
    with io.open(inPath, "r", encoding="utf-8") as f:
        next(f)
        for line in f:
            strsplit = list( map ( str.strip, line.split('\t') ) )
            if len(strsplit) >= 7:
                protein, partner, protein_pos, chain_pos, pdbid, mut, partnerChain = strsplit[:7]
                ppi = '-'.join( [protein, partner, protein_pos, mut[-1]] )
                if ppi not in ppis:
                    if len(strsplit) == 7:
                        ppis.add(ppi)
                        struc = (pdbid,) + tuple(sorted([mut[1], partnerChain]))
                        if struc in mutations:
                            mutations[struc].add(mut)
                        else:
                            mutations[struc] = { mut }
                    elif len(strsplit) == 9:
                        ppis.add(ppi)
    return {k:list(v) for k, v in mutations.items()}

def write_hpc_job (outPath,
                   nodes = 1,
                   ppn = 1,
                   pmem = 7700,
                   walltime = '1:00:00:00',
                   outputfile = 'outputfile',
                   errorfile = 'errorfile',
                   rapid = None,
                   jobid = None,
                   commands = None):
    
    with io.open(outPath, "w") as fout:
        fout.write('#!/bin/bash')
        fout.write('\n' + '#PBS -l nodes=%d:ppn=%d,pmem=%dm,walltime=%s' % (nodes, ppn, pmem, walltime))
        if rapid:
            fout.write('\n' + '#PBS -A %s' % rapid)
        fout.write('\n' + '#PBS -o %s' % outputfile)
        fout.write('\n' + '#PBS -e %s' % errorfile)
        if jobid:
            fout.write('\n' + '#PBS -N %s' % jobid)
        fout.write('\n\n' + 'cd $PBS_O_WORKDIR')
        if commands:
            for cmd in commands:
                fout.write('\n' + cmd)

def read_bindprofx_results (inDir):
    
    processed = {}
    unprocessed = {}
    strucDir = os.listdir(inDir)
    strucDir = [dir for dir in strucDir if os.path.isdir(inDir / dir)]
    for strucID in strucDir:
        struc = tuple(strucID.split('_'))
        if re.match(r'\D\S\d+\D', struc[-1]):
            struc = struc[:-1]
        
        results = []
        resultFile = inDir / strucID / 'result.txt'
        if resultFile.is_file():
            with io.open(resultFile, "r") as f:
                results = list( map(str.strip, f.read().split(';')) )
            results.remove('')
            for result in results:
                ddg, mut = result.split()
                processed[struc + (mut,)] = float(ddg)
        
        if not results:
            mutListFile = inDir / strucID / 'mutList.txt'
            with open(mutListFile, 'r') as f:
                mutList = list( map(str.strip, f.read().split(';')) )
            mutList.remove('')
            if len(mutList) == 1:
                processed[struc + (mutList.pop(),)] = 'X'
            elif len(mutList) > 1:
                errorFile = inDir / strucID / 'errorfile'
                with open(errorFile, 'r') as f:
                    error = f.read()
                if (('FileNotFoundError: [Errno 2] No such file or directory:' in error) and
                    ('align.out' in error)):
                    for mut in mutList:
                        processed[struc + (mut,)] = 'X'
                else:
                    for mut in mutList:
                        unprocessed[struc + (mut,)] = [ mut ]
    
    return processed, unprocessed
