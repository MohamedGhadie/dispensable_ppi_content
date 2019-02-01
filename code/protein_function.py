import pickle
import pandas as pd
import numpy as np
from scipy.stats.stats import pearsonr
from simple_tools import valid_uniprot_id

def partner_sim (p1, p2, partners):
    """Calculate the fraction of interaction partners shared by two proteins.

    Args:
        p1 (str): protein 1 ID.
        p2 (str): protein 2 ID.
        partners (dict): dictionary containing list of partners for each protein.
    
    Returns:
        float: fraction of shared partners if at least one protein has a partner, 
                otherwise NaN.

    """
    partners_1 = (partners[p1] if p1 in partners else set()) - {p2}
    partners_2 = (partners[p2] if p2 in partners else set()) - {p1}
    directPartners_1 = set(partners_1)
    for p in directPartners_1:
        if p in partners:
            partners_1 = partners_1 | partners[p]
    directPartners_2 = set(partners_2)
    for p in directPartners_2:
        if p in partners:
            partners_2 = partners_2 | partners[p]
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

def coexpr (p1, p2, expr, minTissues = 3):
    """Calculate the tissue co-expression for two proteins using Pearson's correlation
        coefficient.

    Args:
        p1 (str): protein 1 ID.
        p2 (str): protein 2 ID.
        expr (dict): dictionary containing tissue expression values for each protein.
        minTissues (int): minimum required number of tissues with expression levels 
                            defined for both proteins.
    
    Returns:
        float: tissue co-expression if both proteins have expression values 
                defined in at least 'minTissues' tissues together, otherwise NaN.
    
    """
    if (p1 in expr) and (p2 in expr):
        e1 = expr[p1]
        e2 = expr[p2]
        not_nan = (np.isnan(e1) | np.isnan(e2)) == False
        if sum(not_nan) >= minTissues:
            corr, p = pearsonr(e1[not_nan], e2[not_nan])
            return corr
    return np.nan

def produce_protein_go_dictionaries (inPath,
                                     GO_outPath,
                                     MF_outPath,
                                     BP_outPath,
                                     CC_outPath):
    """Make dictionaries of protein gene ontology (GO) associations, with each GO 
        aspect in a separate dictionary.

    Args:
        inPath (str): file directory of all protein GO associations.
        GO_outPath (str): file directory to save single dict for all GO aspects.
        MF_outPath (str): file directory to save single dict for Molecular Function (F).
        BP_outPath (str): file directory to save single dict for Biological Process (P).
        CC_outPath (str): file directory to save single dict for Cellular Component (C).

    """
    produce_protein_go_dict(inPath,
                            GO_outPath)
    
    produce_protein_go_dict(inPath,
                            MF_outPath,
                            aspect = 'F')
    
    produce_protein_go_dict(inPath,
                            BP_outPath,
                            aspect = 'P')
    
    produce_protein_go_dict(inPath,
                            CC_outPath,
                            aspect = 'C')

def produce_protein_go_dict(inPath, outPath, aspect = None):
    """Make a dictionary of protein gene ontology (GO) associations.

    Args:
        inPath (str): file directory of protein GO associations.
        outPath (str): file directory to save output dict to.
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

def get_all_go_terms(inPath, aspect = None):
    """Return a list of all GO terms.

    Args:
        inPath (str): file directory of protein GO associations.
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
    if aspect is None:
        return list(goa["go_id"])
    else:
        return list(goa.loc[goa["go_aspect"] == aspect, "go_id"])

def produce_protein_expr_dict (inPath,
                               uniprotIDmapFile,
                               outPath,
                               headers = None):
    """Make a dictionary of protein tissue expression, log2 transformed.

    Args:
        inPath (str): file directory of gene tissue expression data.
        uniprotIDmapFile (Path): file directory containing dict of mappings to UniProt IDs.
        outPath (str): file directory to save output dict to.
        headers (list): list of expression column indices starting with gene name column index
                        followed by indices for expression data. If None, column 1 is used as
                        gene name and columns 2 to 17 are used as expression data.

    """
    with open(uniprotIDmapFile, 'rb') as f:
        uniprotID = pickle.load(f)
    expr = pd.read_table(inPath, comment='#', sep="\t")
    if headers is None:
        headers = list(range(1, 18))
    e = {}
    for _, row in expr.iterrows():
       if row[headers[0]] in uniprotID:
            id = uniprotID[row[headers[0]]]
            e[id] = np.log2(np.array([row[k] for k in headers[1:]]))
    with open(outPath, 'wb') as fOut:
        pickle.dump(e, fOut)
