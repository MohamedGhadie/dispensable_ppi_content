#----------------------------------------------------------------------------------------
# Process SKEMPI file.
#
# Run the following scripts before running this script:
# - produce_data_mappings.py
#----------------------------------------------------------------------------------------

import io
import os
import pandas as pd
from pathlib import Path
from Bio import SeqIO
from text_tools import parse_skempi_interactome
from structural_annotation import process_skempi_mutations

def main():
    
    # allow downloading of PDB structures while constructing the structural interactome
    allow_pdb_downloads = True
    
    # suppress PDB warnings when constructing the structural interactome
    suppress_pdb_warnings = True
    
    # parent directory of all data files
    dataDir = Path('../data')
    
    # directory of data files from external sources
    extDir = dataDir / 'external'
    
    # parent directory of all processed data files
    procDir = dataDir / 'processed'
    
    skempiDir = procDir / 'skempi'
    
    # directory for PDB structure files
    pdbDir = Path('../pdb_files')
    
    # input data files
    skempiFile = extDir / 'skempi_v2.csv'
    pdbSeqresFile = extDir / 'pdb_seqres_reduced.fasta'
    chainSeqFile = procDir / 'chain_sequences.pkl'
    chainStrucResFile = procDir / 'chain_strucRes.pkl'
    
    # output data files
    interactomeFile = skempiDir / 'skempi_interactome.txt'
    skempiSeqFile = skempiDir / 'skempi_sequences.fasta'
    processedSkempiFile = skempiDir / 'processed_skempi_mutations.txt'
        
    # create directories if not existing
    if not skempiDir.exists():
        os.makedirs(skempiDir)
    if not pdbDir.exists():
        os.makedirs(pdbDir)
    
    if not interactomeFile.is_file():
        parse_skempi_interactome (skempiFile, interactomeFile)    
    
    if not skempiSeqFile.is_file():
        seqres = {}
        s = list(SeqIO.parse(str(pdbSeqresFile), 'fasta'))
        for _, row in enumerate(s):
            seqres[row.id] = str(row.seq)
    
        with io.open(skempiSeqFile, "w") as fout:
            for id in chainIDs:
                if id in seqres:
                    fout.write('>' + id + '\n')
                    fout.write(seqres[id] + '\n')
    
    if not processedSkempiFile.is_file():
        process_skempi_mutations (skempiFile,
                                  chainSeqFile,
                                  chainStrucResFile,
                                  pdbDir,
                                  processedSkempiFile,
                                  downloadPDB = allow_pdb_downloads,
                                  suppressWarnings = suppress_pdb_warnings)

if __name__ == '__main__':
    main()
