
import os
import pickle
import pandas as pd
from pathlib import Path
from text_tools import read_list_table

def main():
    
    # reference interactome name
    interactome_name = 'experiment'
    
    # parent directory of all data files
    dataDir = Path('../data')
    
    # parent directory of all processed data files
    procDir = dataDir / 'processed'
    
    # directory of processed data files specific to interactome
    interactomeDir = procDir / interactome_name
    
    # input data files
    nondiseaseEdgotypeFile = interactomeDir / 'nondisease_mutation_edgotype_experiment.txt'
    diseaseEdgotypeFile = interactomeDir / 'disease_mutation_edgotype_experiment.txt'
    uniprotMapFile = procDir / 'to_human_uniprotID_map.pkl'
    geneMapFile = procDir / 'to_human_geneName_map.pkl'
    
    # output data files
    neutralPPIFile = interactomeDir / 'neutral_ppis.txt'
    deleteriousPPIFile = interactomeDir / 'deleterious_ppis.txt'
    
    # create output directories if not existing
    if not interactomeDir.exists():
        os.makedirs(interactomeDir)
    
    neutralPerturbs = read_list_table (nondiseaseEdgotypeFile,
                                       cols=["partners", "perturbations"],
                                       dtyp=[str, float])
    neutralPerturbs["Entrez_Gene_ID"] = neutralPerturbs["Entrez_Gene_ID"].apply(str)
    
    neutralPPIs = set()
    for gene, partners, perturbations, edgotype in neutralPerturbs[["Entrez_Gene_ID",
                                                                    "partners",
                                                                    "perturbations",
                                                                    "Edgotype_class"]].values:
        if edgotype == 'Edgetic':
            for p, pert in zip(partners, perturbations):
                if pert == 1:
                    ppi = tuple(sorted([gene, p]))
                    neutralPPIs.add(ppi)
    
    deleteriousPerturbs = read_list_table (diseaseEdgotypeFile,
                                           cols=["partners", "perturbations"],
                                           dtyp=[str, float])
    deleteriousPerturbs["Entrez_Gene_ID"] = deleteriousPerturbs["Entrez_Gene_ID"].apply(str)
    
    deleteriousPPIs = set()
    for gene, partners, perturbations, edgotype in deleteriousPerturbs[["Entrez_Gene_ID",
                                                                        "partners",
                                                                        "perturbations",
                                                                        "Edgotype_class"]].values:
        if edgotype == 'Edgetic':
            for p, pert in zip(partners, perturbations):
                if pert == 1:
                    ppi = tuple(sorted([gene, p]))
                    deleteriousPPIs.add(ppi)
    
    #neutralPPIs = neutralPPIs - (neutralPPIs & deleteriousPPIs)
    
    with open(geneMapFile, 'rb') as f:
        geneMap = pickle.load(f)
    with open(uniprotMapFile, 'rb') as f:
        toUniProtID = pickle.load(f)
    
    genes_1, genes_2 = zip(* neutralPPIs)
    neutralDisruptions = pd.DataFrame(data={"Entrez_Gene_ID_1":genes_1, "Entrez_Gene_ID_2":genes_2})
    neutralDisruptions["Protein_1"] = neutralDisruptions["Entrez_Gene_ID_1"].apply(
                                            lambda x: toUniProtID[x] if x in toUniProtID else '-')
    neutralDisruptions["Protein_2"] = neutralDisruptions["Entrez_Gene_ID_2"].apply(
                                            lambda x: toUniProtID[x] if x in toUniProtID else '-')
    neutralDisruptions["Gene_1"] = neutralDisruptions["Protein_1"].apply(
                                            lambda x: geneMap[x] if x in geneMap else '-')
    neutralDisruptions["Gene_2"] = neutralDisruptions["Protein_2"].apply(
                                            lambda x: geneMap[x] if x in geneMap else '-')
    neutralDisruptions.to_csv(neutralPPIFile, index=False, sep='\t')
    
    genes_1, genes_2 = zip(* deleteriousPPIs)
    deleteriousDisruptions = pd.DataFrame(data={"Entrez_Gene_ID_1":genes_1, "Entrez_Gene_ID_2":genes_2})
    deleteriousDisruptions["Protein_1"] = deleteriousDisruptions["Entrez_Gene_ID_1"].apply(
                                            lambda x: toUniProtID[x] if x in toUniProtID else '-')
    deleteriousDisruptions["Protein_2"] = deleteriousDisruptions["Entrez_Gene_ID_2"].apply(
                                            lambda x: toUniProtID[x] if x in toUniProtID else '-')
    deleteriousDisruptions["Gene_1"] = deleteriousDisruptions["Protein_1"].apply(
                                            lambda x: geneMap[x] if x in geneMap else '-')
    deleteriousDisruptions["Gene_2"] = deleteriousDisruptions["Protein_2"].apply(
                                            lambda x: geneMap[x] if x in geneMap else '-')
    deleteriousDisruptions.to_csv(deleteriousPPIFile, index=False, sep='\t')

if __name__ == "__main__":
    main()
