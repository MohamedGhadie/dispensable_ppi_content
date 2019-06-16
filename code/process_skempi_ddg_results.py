#----------------------------------------------------------------------------------------
# This script processes results from foldx calculations for SKEMPI mutations. Since each 
# foldx result file contains only one mutation, no second-round jobs will be produced.
#
# Requirements:
# Files must be in format produced by script produce_mutation_structure_maps.py
#----------------------------------------------------------------------------------------

import os
from pathlib import Path
from ddg_tools import (read_foldx_results,
                       read_bindprofx_results,
                       produce_foldx_and_beluga_jobs,
                       produce_bindprofx_and_guillimin_jobs,
                       write_mutation_ddg_tofile)

def main():
    
    # use crystal or homology model structures 
    # options: model, crystal
    structure = 'model'
    
    # method of calculating mutation ∆∆G
    # options: bindprofx, foldx
    ddg_method = 'bindprofx'
    
    # parent directory of all data files
    dataDir = Path('../data')
    
    # parent directory of all processed data files
    procDir = dataDir / 'processed'
    
    # directory of SKEMPI data files
    skempiDir = procDir / 'skempi'
    
    # directory of foldx results
    inDir = skempiDir / ddg_method / 'results'
    
    # directory of foldx output jobs
    outDir = skempiDir / ddg_method
    
    # directory of PDB structure files
    pdbDir = Path('../pdb_files')
    
    # output file whose ∆∆G values wil be updated
    mutationsFile = skempiDir / ('skempi_%s_mutations_%s_ddg.txt' % (structure, ddg_method))
    
    tempFile = skempiDir / ('skempi_%s_mutations_%s_ddg_2.txt' % (structure, ddg_method))
    
    # create output directories if not existing
    if not outDir.exists():
        os.makedirs(outDir)
    
    if ddg_method is 'foldx':
        processed, unprocessed = read_foldx_results (inDir)
        produce_foldx_and_beluga_jobs (unprocessed,
                                       pdbDir,
                                       outDir,
                                       account = 'ctb-yxia',
                                       walltime = '1-00',
                                       ntasks = 1,
                                       nodes = 1,
                                       ntasks_per_node = 1,
                                       cpus_per_task = 1,
                                       mem_per_cpu = '8000M',
                                       outputfile = '/project/def-yxia/ghadie84/foldx/data/%x-%j.out',
                                       username = 'ghadie84',
                                       serverDataDir = '/project/def-yxia/ghadie84/foldx/data')
    elif ddg_method is 'bindprofx':
        processed, unprocessed = read_bindprofx_results (inDir)
        produce_bindprofx_and_guillimin_jobs (unprocessed,
                                              pdbDir,
                                              outDir,
                                              nodes = 1,
                                              ppn = 1,
                                              pmem = 7700,
                                              walltime = '1:00:00:00',
                                              rapid = 'evf-115-aa',
                                              username = 'ghadie84',
                                              extraCommands = ['source /home/ghadie84/venv/bin/activate'],
                                              serverDataDir = '/sf1/project/evf-115-aa/ghadie84/bindprofx/data')
    
    write_mutation_ddg_tofile (processed, mutationsFile, tempFile, type = 'binding')
    os.remove (mutationsFile)
    os.rename (tempFile, mutationsFile)

if __name__ == "__main__":
    main()
