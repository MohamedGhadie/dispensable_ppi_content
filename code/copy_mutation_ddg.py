from pathlib import Path
from ddg_tools import copy_mutation_ddg

def main():
    
    dataDir = Path('/Volumes/MG_Samsung/junk_ppi_content/data/processed')
    
    inPath1 = dataDir / 'HI-II-14' / 'disease_mutations_onchains_ddg.txt'
    inPath2 = dataDir / 'IntAct' / 'nondisease_mutations_onchains_ddg.txt'
    outPath = dataDir / 'IntAct' / 'nondisease_mutations_onchains_ddg_2.txt'
    copy_mutation_ddg (inPath1, inPath2, outPath, 'binding')

if __name__ == "__main__":
    main()
