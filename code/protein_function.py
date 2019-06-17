#----------------------------------------------------------------------------------------
# Modules for computations on protein function.
#----------------------------------------------------------------------------------------

import os
import io
import pickle
import pandas as pd
import numpy as np
import subprocess
from scipy.stats.stats import pearsonr
from simple_tools import valid_uniprot_id, is_numeric, hamming_dist

def partner_sim (p1, p2, partners):
    """Calculate the fraction of interaction partners shared by two proteins.

    Args:
        p1 (str): protein 1 ID.
        p2 (str): protein 2 ID.
        partners (dict): dictionary containing set of partners for each protein.
    
    Returns:
        float: fraction of shared partners if at least one protein has a partner, otherwise NaN.

    """
    partners_1 = (partners[p1] if p1 in partners else set()) - {p2}
    partners_2 = (partners[p2] if p2 in partners else set()) - {p1}
    total = len(partners_1 | partners_2)
    if total > 0:
        return len(partners_1 & partners_2) / total
    else:
        return np.nan

def go_sim (p1, p2, goAssoc):
    """Calculate the fraction of gene ontology (GO) terms shared by two proteins.

    Args:
        p1 (str): protein 1 ID.
        p2 (str): protein 2 ID.
        goAssoc (dict): dictionary containing list of GO terms for each protein.
    
    Returns:
        float: fraction of shared GO terms if at least one protein is associated with a 
                GO term, otherwise NaN.
    
    """
    go1 = (goAssoc[p1] if p1 in goAssoc else set())
    go2 = (goAssoc[p2] if p2 in goAssoc else set())
    if (len(go1) > 0) or (len(go2) > 0):
        return len(go1 & go2) / len(go1 | go2)
    else:
        return np.nan

def coexpr (p1, p2, expr, minTissues = 3, method = 'pearson_corr'):
    """Calculate tissue co-expression for two proteins.

    Args:
        p1 (str): protein 1 ID.
        p2 (str): protein 2 ID.
        expr (dict): dictionary containing tissue expression values for each protein.
        minTissues (int): minimum required number of tissues with expression levels 
                            defined for both proteins.
        method (str): method used to calculate coexpression. Set to 'pearson_corr' to return 
                        Pearson's correlation coefficient, or 'hamming_dist' to return 
                        1 - hamming_distance / length of valid columns.
    
    Returns:
        float
    
    """
    if (p1 in expr) and (p2 in expr):
        e1, e2 = expr[p1], expr[p2]
        not_nan = (np.isnan(e1) == False) & (np.isnan(e2) == False)
        numTissues = sum(not_nan)
        if numTissues >= minTissues:
            e1, e2 = e1[not_nan], e2[not_nan]
            if method is 'pearson_corr':
                if (len(set(e1)) > 1) and (len(set(e2)) > 1):
                    corr, p = pearsonr(e1, e2)
                    return corr
            elif method is 'hamming_dist':
                return 1 - hamming_dist (e1, e2) / numTissues
    return np.nan

def produce_protein_go_dictionaries (inPath,
                                     GO_outPath,
                                     MF_outPath,
                                     BP_outPath,
                                     CC_outPath):
    """Make dictionaries of protein gene ontology (GO) associations, with each root GO 
        in a separate dictionary.

    Args:
        inPath (Path): path to file containing all protein GO associations.
        GO_outPath (Path): file path to save dict of all GO terms.
        MF_outPath (Path): file path to save dict of molecular function (F) terms.
        BP_outPath (Path): file path to save dict of biological process (P) terms.
        CC_outPath (Path): file path to save dict of cellular component (C) terms.

    """
    produce_protein_go_dict (inPath, GO_outPath)
    produce_protein_go_dict (inPath, MF_outPath, aspect = 'F')
    produce_protein_go_dict (inPath, BP_outPath, aspect = 'P')
    produce_protein_go_dict (inPath, CC_outPath, aspect = 'C')

def produce_protein_go_dict (inPath, outPath, aspect = None):
    """Make a dictionary of protein gene ontology (GO) associations.

    Args:
        inPath (Path): path to file containing protein GO associations.
        outPath (Path): file path to save output dict to.
        aspect (str): GO aspects to select: P, F, C, or all if not provided.

    """
    goa = pd.read_table(inPath, header=None, sep="\t")
    goa.columns = ["db",
                   "id",
                   "gene_symbol",
                   "qualifier",
                   "go_id",
                   "db_reference",
                   "evid_code",
                   "with_or_from",
                   "go_aspect",
                   "gene_name",
                   "synonym",
                   "product_type",
                   "taxon",
                   "date",
                   "source",
                   "extension",
                   "variant"]
    goa = goa.applymap(lambda x: x.strip() if isinstance(x, str) else x)
    goa = goa[goa["qualifier"].apply(lambda x: ('NOT' not in x) if isinstance(x, str) else True) &
              (goa["taxon"] == 'taxon:9606') &
              goa["id"].apply(valid_uniprot_id) &
              (goa["product_type"] == 'protein')]
    if aspect is not None:
        goa = goa[goa["go_aspect"] == aspect]
    go = {}
    for _, row in goa.iterrows():
        if row.id in go:
            go[row.id].add(row.go_id)
        else:
            go[row.id] = {row.go_id}
    with open(outPath, 'wb') as fOut:
        pickle.dump(go, fOut)

def get_all_go_terms (inPath, aspect = None):
    """Return a list of all GO terms.

    Args:
        inPath (Path): path to file containing protein GO associations.
        aspect (str): GO aspects to select: P, F, C, or all if not provided.

    Returns:
        list

    """
    goa = pd.read_table(inPath, header=None, sep="\t")
    goa.columns = ["db",
                   "id",
                   "gene_symbol",
                   "qualifier",
                   "go_id",
                   "db_reference",
                   "evid_code",
                   "with_or_from",
                   "go_aspect",
                   "gene_name",
                   "synonym",
                   "product_type",
                   "taxon",
                   "date",
                   "source",
                   "extension",
                   "variant"]
    goa = goa.applymap(lambda x: x.strip() if isinstance(x, str) else x)
    if aspect is None:
        return list(goa["go_id"])
    else:
        return list(goa.loc[goa["go_aspect"] == aspect, "go_id"])

def produce_fastsemsim_protein_gosim_dict (inPath,
                                           outPath,
                                           task = 'SS',
                                           ont_type = 'GeneOntology',
                                           sim_measure = 'SimGIC',
                                           mix_method = 'BMA',
                                           cutoff = -1,
                                           remove_nan = True,
                                           query_mode = 'pairs',
                                           query_ss_type = 'obj',
                                           ac_species = 'human',
                                           ont_root = 'biological_process',
                                           ont_ignore = None,
                                           ec_ignore = None,
                                           verbosity = '-vv',
                                           ontologyFile = 'go-basic.obo',
                                           annotationFile = None,
                                           paramOutFile = 'fastsemsim_parameters',
                                           fastsemsimOutFile = 'fastsemsim_output'):
    """Calculate gene ontology (GO) similarity using fastsemsim library.
        See https://pythonhosted.org/fastsemsim/
    Args:
        inPath (Path): path to file containing list of proteins (pairs or singles).
        outPath (Path): path to file where protein GO similarity dictionary is saved to.
        task (str): calculation to perform.
        ont_type (str): type of ontology.
        sim_measure (str): measure used to calculate GO similarity.
        mix_method (str): mixing strategy used by pairwise semantic similarity measures.
        cutoff (numeric): filter out GO similarity results below this cutoff.
        remove_nan (bool): remove NaN values from results.
        query_mode (str): input file format; either pairs or list of proteins.
        query_ss_type (str): query type; ex objects annotated with ontology terms.
        ac_species (str): species name.
        ont_root (str): root ontology; biological_process, molecular function or cellular component.
        ont_ignore (list): relationships to ignore.
        ec_ignore (list): evidence codes to ignore.
        verbosity (str): verbosity level.
        ontologyFile (Path): path to gene ontology file.
        annotationFile (Path): path to gene ontology associations.
        paramOutFile (Path): path to save Fastsemsim input parameters.
        fastsemsimOutFile (Path): path to save Fastsemsim output.
    
    """
    if annotationFile:
        cmd = ['fastsemsim', '--ac_file', annotationFile]
    else:
        cmd = ['fastsemsim', '--ac_species', ac_species]
    cmd += [verbosity, '--task', task, '-o', ont_type, '--ontology_file', ontologyFile, 
            '--query_input', 'file', '--query_mode', query_mode, '--query_ss_type', query_ss_type, 
            '--query_file', inPath, '--tss', sim_measure, '--tmix', mix_method, '--root', 
            ont_root, '--cut', str(cutoff), '--remove_nan', '--save_params', paramOutFile, 
            '--output_file', fastsemsimOutFile]
    if ec_ignore:
        for ec in ec_ignore:
            cmd += ['--ignore_EC', ec]
    if ont_ignore:
        for ont in ont_ignore:
            cmd += ['--ontology_ignore', ont]
    if remove_nan:
        cmd.append('--remove_nan')
    print(' '.join(map(str, cmd)))
    subprocess.run(cmd)
    table = pd.read_table(fastsemsimOutFile, sep='\t')
    gosim = {tuple(sorted((p1, p2))):sim for p1, p2, sim in table.values}
    with open(outPath, 'wb') as fOut:
        pickle.dump(gosim, fOut)

def produce_illumina_expr_dict (inPath,
                                uniprotIDmapFile,
                                outPath,
                                headers = None):
    """Make a dictionary of protein tissue expression from Illumina Body Map dataset, log2 transformed.

    Args:
        inPath (Path): path to file containing Illumina Body Map tissue expression data.
        uniprotIDmapFile (Path): path to file containing dict of mappings to UniProt IDs.
        outPath (Path): file path to save output dict to.
        headers (list): list of expression column indices starting with gene name column index
                        followed by indices for expression data. If None, column 1 is used as
                        gene name and columns 2 to 17 are used as expression data.

    """
    with open(uniprotIDmapFile, 'rb') as f:
        uniprotID = pickle.load(f)
    expr = pd.read_table(inPath, comment='#', sep="\t")
    if not headers:
        headers = list(range(1, 18))
    e = {}
    for _, row in expr.iterrows():
       if row[headers[0]] in uniprotID:
            id = uniprotID[row[headers[0]]]
            e[id] = np.log2(np.array([row[k] for k in headers[1:]]))
    with open(outPath, 'wb') as fOut:
        pickle.dump(e, fOut)

def produce_gtex_expr_dict (inDir,
                            uniprotIDmapFile,
                            outPath,
                            uniprotIDlistFile = None):
    """Make a dictionary of protein tissue expression from the GTEx dataset.

    Args:
        inDir (Path): file directory containing GTEx tissue expression data files.
        uniprotIDmapFile (Path): path to file containing dict of mappings to UniProt IDs.
        outPath (Path): file path to save output dict to.
        uniprotIDlistFile (Path): path to file containing list of UniProt IDs.

    """
    with open(uniprotIDmapFile, 'rb') as f:
        uniprotID = pickle.load(f)
    if uniprotIDlistFile:
        with open(uniprotIDlistFile, 'r') as f:
            uniprotIDlist = list(set(f.read().split()))
    else:
        uniprotIDlist = list(uniprotID.values())
    tissueExpr = {k:[] for k in uniprotIDlist}
    filenames = [f for f in os.listdir(inDir) if f.endswith('.bed') and not f.startswith('.')]
    for i, filename in enumerate(filenames):
        print('processing file %d of %d: %s' % (i + 1, len(filenames), filename))
        expr = {}
        inPath = inDir / filename
        with io.open(inPath, "r", errors='replace') as f:
            next(f)
            for j, line in enumerate(f):
                linesplit = list(map(str.strip, line.strip().split('\t')))
                if len(linesplit) > 4:
                    chr, start, end, gene_id = linesplit[:4]
                    gene_id = gene_id.split('.')[0]
                    if gene_id in uniprotID:
                        id = uniprotID[gene_id]
                        expr[id] = np.nanmean([float(e) for e in linesplit[4:] if is_numeric(e)])
        for id in uniprotIDlist:
            tissueExpr[id].append(expr[id] if id in expr else np.nan)
        print('%d lines processed' % (j + 1))
    for k, v in tissueExpr.items():
        tissueExpr[k] = np.array(v)
    with open(outPath, 'wb') as fOut:
        pickle.dump(tissueExpr, fOut)

def produce_hpa_expr_dict (inPath, uniprotIDmapFile, outPath):
    """Make a dictionary of protein tissue expression from the HPA dataset.

    Args:
        inPath (Path): path to file containing HPA tissue expression.
        uniprotIDmapFile (Path): path to file containing dict of mappings to UniProt IDs.
        outPath (Path): file path to save output dict to.

    """
    with open(uniprotIDmapFile, 'rb') as f:
        uniprotID = pickle.load(f)
    matrix = pd.read_table(inPath, sep='\t')
    matrix = matrix [matrix['Reliability'] != 'Uncertain']
    matrix["unique_tissue"] = matrix["Tissue"] + '_' + matrix["Cell type"]
    matrix.rename(columns={"Gene name":"Gene_name"}, inplace=True)
    
    tissue_expr = {}
    IDs = set()
    for _, row in matrix.iterrows():
        if row.Gene in uniprotID:
            id = uniprotID[row.Gene]
        elif row.Gene_name in uniprotID:
            id = uniprotID[row.Gene_name]
        else:
            id = row.Gene_name
        IDs.add(id)
        tissue_expr[(id, row.unique_tissue)] = row.Level
    
    unique_tissues = list(set(matrix["unique_tissue"]))
    expr_vectors = {}
    for id in IDs:
        expr = []
        for tissue in unique_tissues:
            k = id, tissue
            expr.append(tissue_expr[k] if k in tissue_expr else '-')
        expr_vectors[id] = expr
        
    with open(outPath, 'wb') as fOut:
        pickle.dump(expr_vectors, fOut)

def produce_fantom5_expr_dict (inPath,
                               uniprotIDmapFile,
                               outPath,
                               sampleTypes = None,
                               sampleTypeFile = None,
                               uniprotIDlistFile = None):
    """Make a dictionary of protein tissue expression from the Fantom5 dataset.

    Args:
        inDir (Path): file directory containing Fantom5 tissue expression data files.
        uniprotIDmapFile (Path): path to file containing dict of mappings to UniProt IDs.
        outPath (Path): file path to save output dict to.
        sampleTypes (str): type of sample; ex tissues.
        sampleTypeFile (Path): path to Fantom5 sample type spreadsheet.
        uniprotIDlistFile (Path): path to file containing list of UniProt IDs for output.

    """
    with open(uniprotIDmapFile, 'rb') as f:
        uniprotID = pickle.load(f)
    if uniprotIDlistFile:
        with open(uniprotIDlistFile, 'r') as f:
            uniprotIDlist = list(set(f.read().split()))
    else:
        uniprotIDlist = list(uniprotID.values())
    
    if sampleTypes:
        sampleCategory = pd.read_excel(sampleTypeFile)
        if isinstance(sampleTypes, str):
            sampleTypes = [sampleTypes]
        selcols = sampleCategory.loc[sampleCategory["Sample category"].apply(lambda x: x in sampleTypes), "FF ontology id"]
        selcols = list(set(selcols.values))
    else:
        selcols = None
    
    with io.open(inPath, "r") as f:
        while True:
            line = f.readline()
            if line.startswith('##ParemeterValue[genome_assembly]='):
                geneAssembly = line.strip().split('=')[1]
                break
    
    with io.open(inPath, "r") as f:
        i, hgncIndex, tpmIndex, tissues = -1, -1, [], []
        while True:
            line = f.readline()
            if (line == '') or (line[:2] != '##'):
                break
            i += line.startswith('##ColumnVariables')
            if line.startswith('##ColumnVariables[tpm'):
                if selcols:
                    for col in selcols:
                        if '.%s.%s' % (col, geneAssembly) in line:
                            tpmIndex.append(i)
                            tissues.append(col)
                            break
                else:
                    tpmIndex.append(i)
                    tissues.append(line.strip()[22:].split(']', maxsplit=1)[0])
            elif line.startswith('##ColumnVariables[hgnc_id]'):
                hgncIndex = i
    print('%d columns selected: %s' % (len(tissues), tissues))
    
    print('processing expression values')
    with io.open(inPath, "r") as f:
        expr, line = {}, '##'
        while line.startswith('##'):
            line = f.readline()
        while True:
            if line is '':
                break
            linesplit = list(map(str.strip, line.strip().split('\t')))
            if len(linesplit) > 7:
                hgncID = linesplit[hgncIndex]
                if hgncID in uniprotID:
                    id = uniprotID[hgncID]
                    if id in uniprotIDlist:
                        tpms = [tpm for i, tpm in enumerate(linesplit) if i in tpmIndex]
                        tpms = list(map(float, tpms))
                        if id in expr:
                            expr[id].append(tpms)
                        else:
                            expr[id] = [tpms]
            line = f.readline()
    
    for k, v in expr.items():
        expr[k] = np.nanmean(v, axis=0)
    with open(outPath, 'wb') as fOut:
        pickle.dump(expr, fOut)
