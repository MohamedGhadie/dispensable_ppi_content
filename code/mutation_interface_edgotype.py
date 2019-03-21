import io
import numpy as np
import pandas as pd
from random import seed, random
from simple_tools import reverseTuples
    
def mutation_PPI_interface_perturbations (mutations,
                                          interactome,
                                          maxInterfaces = np.inf,
                                          dist = 0):
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
    perturbations = []
    for protein, mut_pos in mutations[["protein", "mut_position"]].values:
        perturbs = single_mutation_PPI_perturbs(protein,
                                                mut_pos,
                                                interactome,
                                                maxInterfaces = maxInterfaces,
                                                dist = dist)
        perturbations.append(perturbs)
    return perturbations
#     return mutations.apply(lambda x:
#                            single_mutation_PPI_perturbs(x["Protein"],
#                                                         interactome,
#                                                         maxInterfaces,
#                                                         x["Mutation_Position"],
#                                                         dist = dist),
#                            axis=1)

def single_mutation_PPI_perturbs (protein,
                                  pos,
                                  interactome,
                                  maxInterfaces = np.inf,
                                  dist = 0):
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
                                                                pos,
                                                                x["Interfaces"],
                                                                dist = dist),
                                    axis=1)
        return PPIs["Protein_2"].tolist(), PPIs_Perturbed.values
    else:
        return [], np.array([])

def single_mutation_PPI_perturb (protein, pos, interfaces, dist = 0):
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
    return sum(onInterface)

def energy_based_perturbation (perturbs, ddg, cutoff, probabilistic = False):

    knownDDG = unknownDDG = 0
    pertProb = sum([d > cutoff for _, _, _, _, d in ddg.values()]) / len(ddg)
    seed()
    allPred = []
    for protein, pos, mut_res, partners, perturbations in perturbs[["protein",
                                                                    "mut_position",
                                                                    "mut_res",
                                                                    "partners",
                                                                    "perturbations"]].values:
        pred = []
        for p, pert in zip(partners, perturbations):
            if pert == 0:
                pred.append(0)
            else:
                k = protein, p, pos, mut_res
                if k in ddg:
                    knownDDG += 1
                    pdb_id, chainID, ch_partner, chain_mut, ddgVal = ddg[k]
                    pred.append(1 if ddgVal > cutoff else 0)
                else:
                    unknownDDG += 1
                    if probabilistic:
                        pred.append(1 if random() < pertProb else 0)
                    else:
                        pred.append(np.nan)
        allPred.append(pred)
    return allPred, knownDDG, unknownDDG

def create_perturbed_network (interactome,
                              perturbations,
                              network_outPath,
                              nodeColor_outPath):
    """Create network from perturbed interactome for plotting by Cytoscape software.

    Args:
        interactome (DataFrame): interactome with all edges.
        perturbations (DataFrame): interactome perturbations.
        network_outPath (str): file directory to save edges for Cytoscape plotting.
        nodeColor_outPath (str): file directory to save node colors for Cytoscape plotting.
    
    """
    network = pd.DataFrame()
    edges = [tuple(e) for e in interactome[["Protein_1", "Protein_2"]].values]
    edgeColors = {e:'black' for e in edges}
    
    network["Node_1"], network["Node_2"] = zip(*edges)
    nodes = list(set(network.values.flatten()))
    numMut = {p:0 for p in nodes}
    network["Edge"] = network["Node_1"] + ' (pp) ' + network["Node_2"]
    
    # assign colors to perturbed edges
    for _, row in perturbations.iterrows():
        if row.protein in numMut:
            numMut[row.protein] += 1 
        for partner, perturb in zip(row.partners, row.perturbations):
            if perturb > 0:
                if (row.protein, partner) in edgeColors:
                    edgeColors[(row.protein, partner)] = 'red'
                elif (partner, row.protein) in edgeColors:
                    edgeColors[(partner, row.protein)] = 'red'
    network["Edge_color"] = [edgeColors[e] for e in edges]
    network.to_csv(network_outPath, index=False, sep='\t')

    # assign colors to nodes
    with io.open(nodeColor_outPath, "w") as fout:
        fout.write('Node' + '\t' + 'Color' + '\n')
        nodeColors = []
        for n in nodes:
            if numMut[n] > 9:
                nodeColors.append('blue')
            elif numMut[n] > 2:
                nodeColors.append('mediumpurple')
            elif numMut[n] > 0:
                nodeColors.append('magenta')
            else:
                nodeColors.append('grey')
            fout.write(n + '\t' + nodeColors[-1] + '\n')
    
    return nodes, edges, nodeColors, network["Edge_color"].tolist()

def assign_edgotypes (perturbs, mono_edgetic = False):
    
    edgotypes = []
    for pert in perturbs:
        etype = assign_edgotype (pert, mono_edgetic = mono_edgetic)
        edgotypes.append(etype)
    return edgotypes

def assign_edgotype (perturbs, mono_edgetic = False):
    
    pert = []
    for p in perturbs:
        pert.append(p > 0 if not np.isnan(p) else np.nan)
    numPerturbs = pert.count(True)
    if (numPerturbs >= 2) or ((numPerturbs == 1) and not mono_edgetic):
        return 'edgetic'
    elif pert.count(np.nan) > 0:
        return '-'
    elif numPerturbs == 1:
        return 'mono-edgetic'
    else:
        return 'non-edgetic'
