import os
import io
import csv
import time
import pickle
from pathlib import Path
from Bio import Seq, SeqIO
import pandas as pd
import numpy as np

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

def produce_item_list(inPath, cols, outPath):
    """Write unique items from a dataframe column to file

    Args:
        inPath (str): file directory containing dataframe
        col (str): name of column to save
        outPath (str): file directory to save list of items to

    """
    df = pd.read_table(inPath, sep='\t')
    items = sorted(set(df[cols].values.flatten()))
    with open(outPath, 'w') as fout:
        for item in items:
            fout.write('%s\n' % item)

def read_list_table (inPath, cols, dtyp, delm = '\t'):
    
    df = pd.read_table(inPath, delm)
    if not isinstance(cols, (list, tuple)):
        cols = [cols]
        dtyp = [dtyp]
    for i, col in enumerate(cols):
        df[col] = df[col].apply(lambda x: list(map(dtyp[i], map(str.strip, x.split(',')))))
    return df

def write_list_table (df, cols, outPath, delm = '\t'):
    
    if isinstance(cols, (list, tuple)):
        for col in cols:
            df[col] = df[col].apply(lambda x: ','.join(map(str, x)))
    else:
        df[cols] = df[cols].apply(lambda x: ','.join(map(str, x)))
    df.to_csv(outPath, index=False, sep=delm)

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
