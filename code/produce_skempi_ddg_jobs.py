#----------------------------------------------------------------------------------------
# This script produces jobs for SKEMPI mutations mapped onto structures to be submitted 
# to FoldX or BindProfX for ∆∆G calculations. Mutations with existing ∆∆G values in the 
# input file are skipped.
#
# Run the following script at least once before running this script:
# - produce_skempi_mutation_structure_maps.py
#----------------------------------------------------------------------------------------

import os
from pathlib import Path
from ddg_tools import (read_unprocessed_ddg_mutations,
                       produce_foldx_and_beluga_jobs,
                       produce_bindprofx_and_guillimin_jobs)

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
    
    # directory of foldx output jobs
    outDir = skempiDir / ddg_method
    
    # directory of PDB structure files
    pdbDir = Path('../pdb_files')
    
    # input file containing mutations to submit to bindprofx
    mutationsFile = skempiDir / ('skempi_%s_mutations_%s_ddg.txt' % (structure, ddg_method))
    
    # create output directories if not existing
    if not outDir.exists():
        os.makedirs(outDir)
    
    mutations = read_unprocessed_ddg_mutations (mutationsFile, 'binding')
    
    if ddg_method is 'foldx':
        produce_foldx_and_beluga_jobs (mutations,
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
        produce_bindprofx_and_guillimin_jobs (mutations,
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

if __name__ == "__main__":
    main()
