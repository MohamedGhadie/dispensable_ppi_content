import io
import re
import pickle
from Bio import Seq
import pandas as pd
from simple_tools import first_mismatch, find_substring

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

def match_hgmd_flanking_sequences (inPath, sequenceFile, outPath):

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
    sequenceMatch = sum(mutations["Seq_match"])
    print( '\t' + '%d mutations matching to sequence (%.1f %%)'
            % (sequenceMatch, 100 * sequenceMatch / len(mutations)) )
    
    # keep only mutations matching to UniProt sequence
    mutations = mutations[mutations["Seq_match"]]
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

def base_change_pos (codonChange):

    str1, str2 = map(str.strip, codonChange.split('-'))
    return first_mismatch(str1, str2)

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
