from pathlib import Path
from text_tools import copy_mutation_ddg

def main():
    
    dataDir = Path('/Volumes/MG_Samsung/junk_ppi_content/data/processed')
    
    inPath1 = dataDir / 'allinteractome_mutations_onchains_ddg.txt'
    inPath2 = dataDir / 'IntAct' / 'disease_mutations_onchains.txt'
    outPath = dataDir / 'IntAct' / 'disease_mutations_onchains_ddg.txt'
    copy_mutation_ddg (inPath1, inPath2, outPath, 'binding')

if __name__ == "__main__":
    main()
