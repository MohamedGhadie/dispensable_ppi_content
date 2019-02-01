import pandas as pd
import numpy as np
from random import seed, random
from simple_tools import reverseTuples

def single_mutation_PPI_perturb (protein, interfaces, pos, dist):
    """Calculate the fraction of a protein's interaction interfaces in whose neighborhood 
        a mutation occurs. Neighborhood is defined within a specified number of residues. 

    Args:
        protein (str): protein ID.
        interfaces (list): protein interfaces in tuple form, with second tuple element
                            being the partner's interface.
        pos (int): mutation position on protein sequence.
        dist (int): max distance measured by residue number for a mutation to be considered
                    in the neighborhood of an interface.
    
    Returns:
        float: fraction of interfaces perturbed by mutation.
    
    """
    leftside = [i for i, j in interfaces]
    mutationNeighbor = range(pos - dist, pos + dist + 1)
    onInterface = []
    for interface in leftside:
        onInterface.append(any([i in mutationNeighbor for i in interface]))
    return sum(onInterface) / len(onInterface)
    
def single_mutation_PPI_perturbs (protein, interactome, maxInterfaces, pos, dist):
    """Calculate the fraction of a protein's interaction interfaces with multiple partners 
        in whose neighborhood a mutation occurs. Neighborhood is defined within a specified  
        number of residues. 

    Args:
        protein (str): protein ID.
        interactome (dataframe): protein interactions annotated with interaction interfaces.
        maxInterfaces (int): select only interaction partners with which the protein has 
                                this number of interfaces or less.
        pos (int): mutation position on protein sequence.
        dist (int): max distance measured by residue number for a mutation to be considered
                    in the neighborhood of an interface.
    
    Returns:
        list: fraction of interfaces perturbed by mutation per interaction partner.
    
    """
    PPIs = interactome[(interactome[["Protein_1","Protein_2"]]==protein).any(1)
                       & interactome["Interfaces"].apply(lambda x:
                                                         0 < len(x) <= maxInterfaces)
                      ].reset_index(drop=True)
    if not PPIs.empty:
        inCol2 = (PPIs["Protein_2"] == protein)
        if sum(inCol2) > 0:
            PPIs.loc[inCol2, "Protein_2"] = PPIs.loc[inCol2, "Protein_1"].values
            PPIs.loc[inCol2, "Protein_1"] = protein
            PPIs.loc[inCol2, "Interfaces"] = PPIs.loc[inCol2, "Interfaces"].apply(reverseTuples)
        
        PPIs_Perturbed = PPIs.apply(lambda x:
                                    single_mutation_PPI_perturb(x["Protein_1"],
                                                                x["Interfaces"],
                                                                pos,
                                                                dist),
                                    axis=1)
        return (list(PPIs["Protein_2"]), PPIs_Perturbed.values)
    else:
        return ([], np.array([]))
    
def mutation_PPI_interface_perturbations (mutations, interactome, maxInterfaces, dist):
    """Calculate the fraction of protein interaction interfaces in whose neighborhood 
        specified mutations occur. Neighborhood is defined within a specified number of 
        residues.

    Args:
        mutations (dataframe): table of mutations with their associated protein IDs.
        interactome (dataframe): protein interactions annotated with interaction interfaces.
        maxInterfaces (int): select only interaction partners with which the protein has 
                                this number of interfaces or less.
        dist (int): max distance measured by residue number for a mutation to be considered
                    in the neighborhood of an interface.
    
    Returns:
        Series: list of fraction of interfaces perturbed per interaction partner for each
                mutation.
    
    """
    return mutations.apply(lambda x:
                           single_mutation_PPI_perturbs(x["Protein"],
                                                        interactome,
                                                        maxInterfaces,
                                                        x["Mutation_Position"],
                                                        dist),
                           axis=1)

def energy_based_perturb_prediction ( perturbs, ddg, cutoff ):

    pertProb = sum( ddg["DDG"] > cutoff ) / len( ddg )
    knownDDG = unknownDDG = 0
    allPred = []
    seed()
    for i, row in perturbs.iterrows():
        pred = []
        for p, pert in zip( row.partners, row.perturbations ):
            if pert == 0:
                pred.append( 0 )
            else:
                ddgVal = ddg.loc[ ( ddg["protein"] == row.Protein ) & 
                                  ( ddg["partner"] == p ) &
                                  ( ddg["protein_pos"] == row.Mutation_Position ),
                                  "DDG" ]
                if not ddgVal.empty:
                    knownDDG += 1
                    pred.append( 1 if ddgVal.item() > cutoff else 0 )
                else:
                    unknownDDG += 1
                    if random() >= pertProb:
                        pred.append( 0 )
                    else:
                        pred.append( 1 )
        allPred.append( np.array( pred ) )
    return pd.Series( allPred ), knownDDG, unknownDDG

def energy_based_perturb_prediction2 ( perturbs, ddg, cutoff ):

    knownDDG = unknownDDG = 0
    allPred = []
    for i, row in perturbs.iterrows():
        pred = []
        for p, pert in zip( row.partners, row.perturbations ):
            if pert == 0:
                pred.append( 0 )
            else:
                ddgVal = ddg.loc[ ( ddg["protein"] == row.Protein ) & 
                                  ( ddg["partner"] == p ) &
                                  ( ddg["protein_pos"] == row.Mutation_Position ),
                                  "DDG" ]
                if not ddgVal.empty:
                    knownDDG += 1
                    pred.append( 1 if ddgVal.item() > cutoff else 0 )
                else:
                    pred.append( np.nan )
        allPred.append( np.array( pred ) )
    return pd.Series(allPred), knownDDG, unknownDDG
