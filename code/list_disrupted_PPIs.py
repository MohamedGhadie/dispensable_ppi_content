
import os
import pandas as pd
from pathlib import Path
from text_tools import read_list_table

def main():
    
    # reference interactome name
    # options: HI-II-14, IntAct
    interactome_name = 'HI-II-14'
    
    # method of calculating mutation ∆∆G for which results will be used
    # options: bindprofx, foldx
    ddg_method = 'bindprofx'
    
    # parent directory of all data files
    dataDir = Path('../data')
    
    # parent directory of all processed data files
    procDir = dataDir / 'processed'
    
    # directory of processed data files specific to interactome
    interactomeDir = procDir / interactome_name
    
    # input data files
    interactomeFile = interactomeDir / 'human_interactome.txt'
    nondiseaseEdgotypeFile = interactomeDir / ('nondisease_mutation_edgotype_physics_%s.txt' % ddg_method)
    diseaseEdgotypeFile = interactomeDir / ('disease_mutation_edgotype_physics_%s.txt' % ddg_method)
    
    # output data files
    neutralPPIFile = interactomeDir / ('neutral_ppis_%s.txt' % ddg_method)
    deleteriousPPIFile = interactomeDir / ('deleterious_ppis_%s.txt' % ddg_method)
    
    # create output directories if not existing
    if not interactomeDir.exists():
        os.makedirs(interactomeDir)
    
    neutralPerturbs = read_list_table (nondiseaseEdgotypeFile,
                                       cols=["partners", "perturbations"],
                                       dtyp=[str, float])
    neutralPPIs = set()
    for protein, partners, perturbations, edgotype in neutralPerturbs[["protein",
                                                                       "partners",
                                                                       "perturbations",
                                                                       "edgotype"]].values:
        if edgotype == 'edgetic':
            for p, pert in zip(partners, perturbations):
                if pert == 1:
                    ppi = tuple(sorted([protein, p]))
                    neutralPPIs.add(ppi)
    
    deleteriousPerturbs = read_list_table (diseaseEdgotypeFile,
                                           cols=["partners", "perturbations"],
                                           dtyp=[str, float])
    deleteriousPPIs = set()
    for protein, partners, perturbations, edgotype in deleteriousPerturbs[["protein",
                                                                           "partners",
                                                                           "perturbations",
                                                                           "edgotype"]].values:
        if edgotype == 'edgetic':
            for p, pert in zip(partners, perturbations):
                if pert == 1:
                    ppi = tuple(sorted([protein, p]))
                    deleteriousPPIs.add(ppi)
    
    neutralPPIs = neutralPPIs - (neutralPPIs & deleteriousPPIs)
    
    interactome = pd.read_table(interactomeFile, sep='\t')
    sel = interactome.apply(lambda x: tuple(sorted([x["Protein_1"], x["Protein_2"]])) in neutralPPIs, axis=1)
    neutralDisruptions = interactome[sel]
    neutralDisruptions.to_csv(neutralPPIFile, index=False, sep='\t')
    
    sel = interactome.apply(lambda x: tuple(sorted([x["Protein_1"], x["Protein_2"]])) in deleteriousPPIs, axis=1)
    deleteriousDisruptions = interactome[sel]
    deleteriousDisruptions.to_csv(deleteriousPPIFile, index=False, sep='\t')

if __name__ == "__main__":
    main()
