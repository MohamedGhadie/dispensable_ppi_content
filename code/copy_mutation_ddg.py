from pathlib import Path
from ddg_tools import copy_mutation_ddg

def main():
    
    dataDir = Path('../data/processed')
    
    fromPath = dataDir / 'IntAct' / 'disease_mutations_foldx_ddg.txt'
    toPath = dataDir / 'IntAct' / 'disease_mutations_onchains.txt'
    outPath = dataDir / 'IntAct' / 'disease_mutations_foldx_ddg_2.txt'
    copy_mutation_ddg (fromPath, toPath, outPath, 'binding')

if __name__ == "__main__":
    main()
