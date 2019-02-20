import os
from pathlib import Path
from ddg_tools import read_bindprofx_results, produce_bindprofx_jobs, write_mutation_ddg_tofile

def main():
    
    # valid reference interactome names
    interactome_names = ['HI-II-14', 'IntAct']
    
    # choose interactome (index in interactome_names)
    interactome_choise = 1
    
    # select reference interactome
    interactome_name = interactome_names[ interactome_choise ]
    
    # directory for input data files
    dataDir = Path('/Volumes/MG_Samsung/junk_ppi_content/data/processed') / interactome_name
    
    # directory for input data files
    inDir = dataDir / 'bindprofx' / 'results_all_12'
    
    # directory for output data files
    outDir = dataDir / 'bindprofx'
    
    # directory for PDB structure files
    pdbDir = Path('/Volumes/MG_Samsung/pdb_files')
    
    
    # create output directories if not existing
    if not outDir.exists():
        os.makedirs(outDir)
    
    processed, unprocessed = read_bindprofx_results (inDir)
    
#     write_mutation_ddg_tofile (processed,
#                                dataDir / 'nondisease_mutations_bindprofx_ddg.txt',
#                                dataDir / 'nondisease_mutations_bindprofx_ddg_2.txt',
#                                'binding')
#     write_mutation_ddg_tofile (processed,
#                                dataDir / 'disease_mutations_bindprofx_ddg.txt',
#                                dataDir / 'disease_mutations_bindprofx_ddg_2.txt',
#                                'binding')
    write_mutation_ddg_tofile (processed,
                               dataDir / 'all_mutations_bindprofx_ddg.txt',
                               dataDir / 'all_mutations_bindprofx_ddg_2.txt',
                               'binding')
    
    produce_bindprofx_jobs (unprocessed,
                            pdbDir,
                            outDir,
                            write_hpc_jobfiles = True,
                            nodes = 1,
                            ppn = 1,
                            pmem = 7700,
                            walltime = '1:00:00:00',
                            rapid = 'evf-115-aa',
                            username = 'ghadie84',
                            hpcCommands = ['source /home/ghadie84/venv/bin/activate'],
                            serverDataDir = '/sf1/project/evf-115-aa/ghadie84/bindprofx/data')

if __name__ == "__main__":
    main()
