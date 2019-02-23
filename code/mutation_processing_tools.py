import io
import re
import time
import pickle
from Bio import Seq
import pandas as pd
import numpy as np
from simple_tools import (is_numeric,
                          first_mismatch,
                          find_substring,
                          find_masked_substring)

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
                                                               'drug-response',
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
    allmut.rename(columns={"residue":"mut_res",
                           "aa_position":"mut_position",
                           "allele":"frequency_reporting_allele",
                           "allele.1": "variation_allele"}, inplace=True)
    allmut.to_csv(outPath, index=False, sep='\t')

def filter_clinvar_mutations (inPath,
                              outPath,
                              assembly = None,
                              origin = None,
                              type = None,
                              incClinSig = None,
                              excClinSig = None,
                              status = None,
                              uniprotIDmapFile = None):
    
    mutations = pd.read_table(inPath, sep='\t')
    numMut = len(mutations)
    print('\t' + 'Total number of mutations: %d' % numMut)
    
    if assembly:
        print('\t' + 'Selecting mutations with genome assembly %s' % assembly)
        if isinstance(assembly, str):
            assembly = {assembly}
        sel = mutations["Assembly"].apply(lambda x: x in assembly)
        mutations = mutations[sel]
        print('\t' + 'Number of mutations selected: %d' % len(mutations))
    
    allorigins = set(mutations["OriginSimple"])
    print('\t' + 'Mutation origins (simplified):')
    for org in allorigins:
        sel = mutations["OriginSimple"] == org
        print('\t\t' + '%s: %d (%.2f %%)' % (org, sum(sel), 100 * sum(sel) / len(mutations)))
    
    if origin:
        print('\t' + 'Selecting mutations of origin %s' % origin)
        if isinstance(origin, str):
            origin = {origin}
        sel = mutations["OriginSimple"].apply(lambda x: x in origin)
        mutations = mutations[sel]
    
    alltypes = set(mutations["Type"])
    print('\t' + 'Mutation types:')
    for tp in alltypes:
        sel = mutations["Type"] == tp
        print('\t\t' + '%s: %d (%.2f %%)' % (tp, sum(sel), 100 * sum(sel) / len(mutations)))
    
    if type:
        print('\t' + 'Selecting mutations of type %s' % type)
        if isinstance(type, str):
            type = {type}
        sel = mutations["Type"].apply(lambda x: x in type)
        mutations = mutations[sel]
    
    clinSig = set(mutations["ClinicalSignificance"])
    clinSigTerms = []
    for c in clinSig:
        clinSigTerms.extend(list(map(str.strip, c.split(','))))
    clinSigTerms = set(clinSigTerms)
    print('\t' + 'Mutation clinical significance terms (overlaps possible):')
    for term in clinSigTerms:
        sel = mutations["ClinicalSignificance"].apply(lambda x: term in x)
        print('\t\t' + '%s: %d (%.1f %%)' % (term, sum(sel), 100 * sum(sel) / len(mutations)))
    
    if incClinSig:
        print('\t' + 'Selecting mutations with clinical significance terms: %s' % incClinSig)
        inc = {incClinSig} if isinstance(incClinSig, str) else set(incClinSig)
        terms = mutations["ClinicalSignificance"].apply(lambda x: set(map(str.strip, x.split(','))))
        sel = terms.apply(lambda x: len(inc & x) > 0)
        mutations = mutations[sel]
        print('\t' + 'Number of mutations selected: %d' % len(mutations))
    if excClinSig:
        print('\t' + 'Removing mutations with clinical significance terms: %s' % excClinSig)
        exc = {excClinSig} if isinstance(excClinSig, str) else set(excClinSig)
        terms = mutations["ClinicalSignificance"].apply(lambda x: set(map(str.strip, x.split(','))))
        sel = terms.apply(lambda x: len(exc & x) == 0)
        mutations = mutations[sel]
        print('\t' + 'Number of mutations selected: %d' % len(mutations))
    
    if status:
        print('\t' + 'Selecting mutations with review status terms: %s' % status)
        if isinstance(status, str):
            status = {status}
        sel = mutations["ReviewStatus"].apply(lambda x: x in status)
        mutations = mutations[sel]
        print('\t' + 'Number of mutations selected: %d' % len(mutations))
    
    if uniprotIDmapFile:
        print('\t' + 'Mapping gene symbols to UniProt IDs')
        with open(uniprotIDmapFile, 'rb') as f:
            uniprotIDmap = pickle.load(f)
        geneNames = mutations["GeneSymbol"].apply(str.upper)
        mutations["protein"] = geneNames.apply(lambda x: uniprotIDmap[x] if x in uniprotIDmap else '-')
    
    mutations.to_csv(outPath, index=False, sep='\t')

def decompose_clinvar_snp_mutations (inPath, outPath):
    
    mutations = pd.read_table(inPath, sep='\t')
    decomposed = zip(* decompose_clinvar_snp_names (mutations["Name"].tolist()))
    for col, val in zip(["rna_acc", "cdna_mut", "wt_res", "mut_position", "mut_res"], decomposed):
        mutations[col] = list(map(lambda x: '-' if x is None else x, val))
    mutations = mutations [(mutations["wt_res"] != '-') & 
                           (mutations["mut_res"] != '-') & 
                           (mutations["mut_position"] != '-')]
    mutations.to_csv(outPath, index=False, sep='\t')

def decompose_clinvar_snp_names (names):
    
    # NM_017547.3(FOXRED1):c.694C>T (p.Gln232Ter)
    decomposed = []
    for name in names:
        decomposed.append( decompose_clinvar_snp_name (name) )
    return decomposed

def decompose_clinvar_snp_name (name):        
    
    m = re.match(r"(\w{2}_\d+\.\d+)*(\(\w*\))*(\:)*(c\.[\w\>\-]*)*(\s)*(\(p\.\w*\))*", name.strip())
    if m:
        rna_acc, gene_name, dots, cdna_mut, space, pr_mut = m.groups()
        if pr_mut:
            wt_res, mut_pos, mut_res = decompose_protein_snp (pr_mut[1:-1])
            return rna_acc, cdna_mut, wt_res, mut_pos, mut_res
        else:
            return rna_acc, cdna_mut, None, None, None

def decompose_protein_snp (mut):
    
    snp = mut.strip()
    m = re.match(r"p\.\D{3}\d+\D{3}$", snp)
    if m:
        wt_res, mut_pos, mut_res = snp[2:5], snp[5:-3], snp[-3:]
        return toOneLetterAA(wt_res), int(mut_pos), toOneLetterAA(mut_res)
    else:
        return None, None, None

def map_clinvar_protein_refseq_IDs (inPath, idMapFile, outPath):
    
    mutations = pd.read_table(inPath, sep='\t')
    with open(idMapFile, 'rb') as f:
        protein_acc = pickle.load(f)
    mutations["prot_acc"] = mutations["rna_acc"].apply(lambda x: protein_acc[x] 
                                                                 if x in protein_acc
                                                                 else '-')
    mutations = mutations [mutations["prot_acc"] != '-']
    mutations.to_csv(outPath, index=False, sep='\t')

def get_flanking_sequences (inPath, sequenceFile, sideLength, outPath):
    
    mutations = pd.read_table(inPath, sep='\t')
    sequences = pd.read_table(sequenceFile, sep='\t')
    sequencedProteins = sequences["ID"].values
    mutations["protein_seq"] = mutations["prot_acc"].apply(lambda x:
                                                          sequences.loc[sequences["ID"]==x, "Sequence"].item() 
                                                          if x in sequencedProteins
                                                          else '-')
    mutations = mutations[mutations["protein_seq"]!='-'].reset_index(drop=True)
    mutations["wt_context"] = mutations.apply(lambda x:
                                              flanking_sequence(x["mut_position"],
                                                                x["protein_seq"],
                                                                sideLength),
                                              axis=1)
    mutations = mutations[mutations["wt_context"]!='-'].reset_index(drop=True)
    mutations["mut_context_pos"], mutations["wt_context"] = zip(* mutations["wt_context"].values)
    mutations.drop("protein_seq", axis=1, inplace=True)
    mutations.to_csv(outPath, index=False, sep='\t')

def match_flanking_sequences (inPath, sequenceFile, outPath, mask = False):
    
    mutations = pd.read_table(inPath, sep='\t')
    sequences = pd.read_table(sequenceFile, sep='\t')
    sequencedProteins = sequences["ID"].values
    mutations["protein_seq"] = mutations["protein"].apply(lambda x:
                                                          sequences.loc[sequences["ID"]==x,
                                                                        "Sequence"].item() 
                                                          if x in sequencedProteins
                                                          else '-')
    if mask:
        mutations["seq_match"] = mutations.apply(lambda x:
                                                 find_masked_substring(x["wt_context"],
                                                                       x["protein_seq"],
                                                                       x["mut_context_pos"]-1),
                                                 axis=1)
    else:
        mutations["seq_match"] = mutations.apply(lambda x:
                                                 find_substring(x["wt_context"],
                                                                x["protein_seq"]),
                                                 axis=1)
    seq_match = mutations.apply(lambda x: x["mut_position"] - x["mut_context_pos"] in x["seq_match"], axis=1)
    print( '\t' + '%d mutations matching to UniProt sequence (%.1f %%)'
            % (sum(seq_match), 100 * sum(seq_match) / len(mutations)) )
    
    # keep only mutations matching to UniProt sequence
    mutations = mutations[seq_match]
    mutations["wt_res"] = mutations.apply(lambda x: x["protein_seq"][x["mut_position"]-1], axis=1)
    mutations["wt_context"] = mutations.apply(lambda x: x["wt_context"][:x["mut_context_pos"]-1] +
                                                        x["wt_res"] + 
                                                        x["wt_context"][x["mut_context_pos"]:], axis=1)
    mutations.drop("protein_seq", axis=1, inplace=True)
    mutations.drop("seq_match", axis=1, inplace=True)
    mutations.to_csv(outPath, index=False, sep='\t')

def flanking_sequence(pos, seq, sideLength):
    
    pos = pos - 1
    seqLen = len(seq)
    if 0 <= pos < seqLen:
        s = max(0, pos - sideLength)
        e = min(seqLen, pos + sideLength + 1)
        return pos - s + 1, seq[ s : e ]
    else:
        return '-'

def remove_synon_nonsense_mutations(inPath, outPath):
    
    mutations = pd.read_table(inPath, sep='\t')
    nonsynonymous = mutations.apply(lambda x: x["mut_res"] != x["wt_res"], axis=1)
    missense = mutations["mut_res"] != '*'
    mutations = mutations[nonsynonymous & missense]
    mutations.to_csv(outPath, index=False, sep='\t')

def toOneLetterAA (aminoAcid):
    
    oneLetter = {'Cis': 'C', 'Asp': 'D', 'Ser': 'S', 'Gln': 'Q', 'Lys': 'K', 'Trp': 'W', 
                 'Thr': 'T', 'Asn': 'N', 'Pro': 'P', 'Phe': 'F', 'Ala': 'A', 'Gly': 'G', 
                 'Ile': 'I', 'Leu': 'L', 'His': 'H', 'Arg': 'R', 'Met': 'M', 'Val': 'V', 
                 'Glu': 'E', 'Tyr': 'Y', 'Ter': '*'}
    
    aa = aminoAcid.title()
    if aa in oneLetter:
        return oneLetter[aa] 
    else:
        return '-'

def remove_mutation_overlaps (naturalMutationsFile, diseaseMutationsFile):
    
    naturalMutations = pd.read_table(naturalMutationsFile, sep='\t')
    diseaseMutations = pd.read_table(diseaseMutationsFile, sep='\t')
    
    # identify common mutations overlapping in position with disease mutations
    df = diseaseMutations[ ["protein", "mut_position"] ].copy()
    df2 = naturalMutations[ ["protein", "mut_position"] ].copy()
    df2 = df2.drop_duplicates().reset_index(drop=True)
    df = df.append(df2, ignore_index=True)
    duplicates = df.duplicated(keep='first')[ len(diseaseMutations) : ]
    duplicates = duplicates.reset_index(drop=True)
    print( '\n' + '%d unique non-disease variants overlap in location with disease mutations' % sum(duplicates) )
    
    # remove common mutations overlapping in position with disease mutations
    for i, row in df2[ duplicates ].iterrows():
        todrop = ( ( naturalMutations["protein"] == row.protein ) & 
                   ( naturalMutations["mut_position"] == row.mut_position ) )
        naturalMutations = naturalMutations[todrop == False].reset_index(drop=True)
    print( 'Overlaping non-disease mutations removed' )
    
    numNaturalMutations = len(naturalMutations)
    numDiseaseMutations = len(diseaseMutations)
    
    # remove invalid mutations, i.e., those with mutation residue similar to wild type
    naturalMutations = naturalMutations [naturalMutations["wt_res"] != naturalMutations["mut_res"]]
    diseaseMutations = diseaseMutations [diseaseMutations["wt_res"] != diseaseMutations["mut_res"]]
    
    print( '\n' + 'Number of invalid mutations removed (WT residue = Mut residue)' )
    print( 'non-disease invalid mutations: %d' % (numNaturalMutations - len(naturalMutations)) )
    print( 'disease invalid mutations: %d' % (numDiseaseMutations - len(diseaseMutations)) )
    
    numNaturalMutations = len(naturalMutations)
    numDiseaseMutations = len(diseaseMutations)
    
    print( '\n' + 'Number of mutations after removing invalid mutations:' )
    print( 'non-disease mutations: %d' % numNaturalMutations )
    print( 'disease mutations: %d' % numDiseaseMutations )
    
    # remove duplicate mutations, by position and mutant residue
    naturalMutations = naturalMutations.drop_duplicates(subset=["protein",
                                                                "mut_position",
                                                                "mut_res"]).reset_index(drop=True)
    diseaseMutations = diseaseMutations.drop_duplicates(subset=["protein",
                                                                "mut_position",
                                                                "mut_res"]).reset_index(drop=True)
    
    numNaturalMutations = len(naturalMutations)
    numDiseaseMutations = len(diseaseMutations)
    
    print( '\n' + 'Number of mutations after removing duplicates by position and mutant residue:' )
    print( 'non-disease mutations: %d' % numNaturalMutations )
    print( 'disease mutations: %d' % numDiseaseMutations )
    
    return naturalMutations, diseaseMutations
