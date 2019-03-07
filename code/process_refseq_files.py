#----------------------------------------------------------------------------------------
# Process RefSeq data files
#----------------------------------------------------------------------------------------

import os
import numpy as np
from pathlib import Path
from text_tools import parse_refSeq_fasta, merge_refSeq_sequence_files

def main():
    
    # number of files to process in RefSeq directory
    numfiles = 39
    
    # parent directory of all data files
    dataDir = Path('../data')
    
    # directory of data files from external sources
    extDir = dataDir / 'external'
    
    # parent directory of all processed data files
    procDir = dataDir / 'processed'
    
    # directory for RefSeq raw data files
    refseqInDir = extDir / 'refSeq'
    
    # directory to save RefSeq processed data files
    refseqOutDir = procDir / 'refSeq_intermediate'
    
    # output data files
    refseqFile = procDir / 'refseq_human_protein_sequences.txt'
    
    if not procDir.exists():
        os.makedirs(procDir)
    if not refseqOutDir.exists():
        os.makedirs(refseqOutDir)
    
    for i in np.arange(1, numfiles + 1):
        print('parsing RefSeq sequence file %d of %d' % (i, numfiles))
        parse_refSeq_fasta (refseqInDir / ('human.%d.protein.faa' % i),
                            refseqOutDir / ('refseq_human_protein_%d.txt' % i))
    
    if not refseqFile.is_file():
        print('merging RefSeq sequences from all files')
        merge_refSeq_sequence_files (refseqOutDir, numfiles, refseqFile)

if __name__ == '__main__':
    main()