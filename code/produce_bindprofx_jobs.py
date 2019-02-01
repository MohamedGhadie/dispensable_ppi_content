import os
from pathlib import Path
from text_tools import read_bindprofx_mutations
from structural_annotation import produce_bindprofx_jobs

def main():
    
    # valid reference interactome names
    interactome_names = ['HI-II-14', 'IntAct']
    
    # choose interactome (index in interactome_names)
    interactome_choise = 1
    
    # select reference interactome
    interactome_name = interactome_names[ interactome_choise ]
    
    # directory for input data files
    inDir = Path('../data/processed') / interactome_name
    
    # directory for output data files
    outDir = inDir / 'bindprofx'
    
    # directory for PDB structure files
    pdbDir = Path('/Volumes/MG_Samsung/pdb_files')
    
    # input file containing mutations to submit to bindprofx
    mutationsFile = inDir / 'all_mutations_onchains_ddg_7.txt'
    
    # create output directories if not existing
    if not outDir.exists():
        os.makedirs( outDir )
    
    mutations = read_bindprofx_mutations (mutationsFile)
    
    produce_bindprofx_jobs (mutations,
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
