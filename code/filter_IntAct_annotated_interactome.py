import os
import pandas as pd
from pathlib import Path
from interactome_tools import (read_single_interface_annotated_interactome,
                               write_single_interface_annotated_interactome)

def main():
    
    # directory where interactome processed data is saved
    interactomeDir = Path('../data/processed/IntAct')
    
    interactomeDupFile = interactomeDir / 'human_interactome_all_withtype.txt'
    interactomeDup = pd.read_table(interactomeDupFile)
    print('\n' + 'Size of reference interactome with duplicates: %d' % len(interactomeDup) )
    
    interactomeFile = interactomeDir / 'human_interactome.txt'
    interactome = pd.read_table(interactomeFile)
    print('\n' + 'Size of reference interactome with no duplicates: %d' % len(interactome) )
    
    interfaceAnnotatedInteractomeFile = interactomeDir / 'human_interface_annotated_interactome.txt'
    annotatedInteractome = read_single_interface_annotated_interactome(interfaceAnnotatedInteractomeFile)
    print('\n' + 'Size of structural interactome: %d' % len( annotatedInteractome ) )
    
#     typedict = {}
#     for i, row in interactome.iterrows():
#         types = set( interactomeDup.loc[ (interactomeDup[ ["Protein_1", "Protein_2"] ] == row.Protein_1).any(1)
#                                          & (interactomeDup[ ["Protein_1", "Protein_2"] ] == row.Protein_2).any(1),
#                                          "Interaction_type" ].values )
#         if 'psi-mi:"MI:0407"(direct interaction)' in types:
#             typ = 'psi-mi:"MI:0407"(direct interaction)'
#         elif 'psi-mi:"MI:0915"(physical association)' in types:
#             typ = 'psi-mi:"MI:0915"(physical association)'
#         elif 'psi-mi:"MI:0914"(association)' in types:
#             typ = 'psi-mi:"MI:0914"(association)'
#         elif not types.issubset( {'psi-mi:"MI:0403"(colocalization)', 'psi-mi:"MI:0208"(genetic interaction)'} ):
#             types = types - {'psi-mi:"MI:0403"(colocalization)', 'psi-mi:"MI:0208"(genetic interaction)'}
#             typ = types.pop()
#         elif 'psi-mi:"MI:0403"(colocalization)' in types:
#             typ = 'psi-mi:"MI:0403"(colocalization)'
#         elif 'psi-mi:"MI:0208"(genetic interaction)' in types:
#             typ = 'psi-mi:"MI:0208"(genetic interaction)'
#         else:
#             typ = '-'
#         interactome.loc[i, "type"] = typ
#         
#         if typ in typedict:
#             typedict[typ] += 1
#         else:
#             typedict[typ] = 1
#     
#     print('\n' + 'Number of PPIs per interaction type in reference interactome')
#     for k in typedict:
#         print( '%s: %d' % (k, typedict[k]) )
    
    typedict = {}
    for i, row in annotatedInteractome.iterrows():
        if i % 100 == 0:
            print( '%d PPIs processed for type' % i )
        types = set( interactomeDup.loc[ (interactomeDup[ ["Protein_1", "Protein_2"] ] == row.Protein_1).any(1)
                                         & (interactomeDup[ ["Protein_1", "Protein_2"] ] == row.Protein_2).any(1),
                                         "Interaction_type" ].values )
        if 'psi-mi:"MI:0407"(direct interaction)' in types:
            typ = 'psi-mi:"MI:0407"(direct interaction)'
        elif 'psi-mi:"MI:0915"(physical association)' in types:
            typ = 'psi-mi:"MI:0915"(physical association)'
        elif 'psi-mi:"MI:0914"(association)' in types:
            typ = 'psi-mi:"MI:0914"(association)'
        elif not types.issubset( {'psi-mi:"MI:0403"(colocalization)', 'psi-mi:"MI:0208"(genetic interaction)'} ):
            types = types - {'psi-mi:"MI:0403"(colocalization)', 'psi-mi:"MI:0208"(genetic interaction)'}
            typ = types.pop()
        elif 'psi-mi:"MI:0403"(colocalization)' in types:
            typ = 'psi-mi:"MI:0403"(colocalization)'
        elif 'psi-mi:"MI:0208"(genetic interaction)' in types:
            typ = 'psi-mi:"MI:0208"(genetic interaction)'
        else:
            typ = '-'
        annotatedInteractome.loc[i, "type"] = typ
        
        if typ in typedict:
            typedict[typ] += 1
        else:
            typedict[typ] = 1
    
    print('\n' + 'Number of PPIs per interaction type in the interface-annotated interactome')
    for k in typedict:
        print( '%s: %d' % (k, typedict[k]) )
    
    outPath = interactomeDir / 'human_interface_annotated_interactome_highconf.txt'
    tothrow = ['psi-mi:"MI:0914"(association)',
               'psi-mi:"MI:0403"(colocalization)',
               'psi-mi:"MI:0208"(genetic interaction)']
    filteredInteractome = annotatedInteractome[ annotatedInteractome["type"].apply(lambda x: x not in tothrow ) 
                                               ].reset_index(drop=True)
    write_single_interface_annotated_interactome ( filteredInteractome, outPath )

if __name__ == "__main__":
    main()
