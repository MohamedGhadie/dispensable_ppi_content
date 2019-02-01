from pathlib import Path
from structural_annotation import copy_mutation_ddg

def main():
    
    # directory for input data files
    inDir = Path('../data/processed')
    
    inPath1 = inDir / 'HI-II-14' / 'all_mutations_onchains_ddg_3.txt'
    inPath2 = inDir / 'IntAct' / 'nondisease_mutations_onchains_ddg_6.txt'
    outPath = inDir / 'IntAct' / 'nondisease_mutations_onchains_ddg_7.txt'
    copy_mutation_ddg (inPath1, inPath2, outPath)
    
    inPath1 = inDir / 'HI-II-14' / 'all_mutations_onchains_ddg_3.txt'
    inPath2 = inDir / 'IntAct' / 'disease_mutations_onchains_ddg_6.txt'
    outPath = inDir / 'IntAct' / 'disease_mutations_onchains_ddg_7.txt'
    copy_mutation_ddg (inPath1, inPath2, outPath)

if __name__ == "__main__":
    main()
