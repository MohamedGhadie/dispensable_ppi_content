
from pathlib import Path
from text_tools import sample_files

def main():
    
    sampleSize = 200
    
    # parent directory of all data files
    dataDir = Path('../data')
    
    # parent directory of all processed data files
    procDir = dataDir / 'processed'
    
    # directory of SKEMPI data files
    skempiDir = procDir / 'skempi'
    
    # directory of foldx output jobs
    inDir = skempiDir / 'foldx_crystals_all'
    
    # directory of foldx output jobs
    outDir = skempiDir / 'foldx_sample'
    
    sample_files (inDir, outDir, sampleSize)

if __name__ == "__main__":
    main()
