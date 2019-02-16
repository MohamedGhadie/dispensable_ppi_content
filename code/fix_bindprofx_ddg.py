import io
from pathlib import Path

def main():
    
    dataDir = Path('/Volumes/MG_Samsung/junk_ppi_content/data/processed')
    inPath = dataDir / 'allinteractome_mutations_onchains_ddg.txt'
    outPath = dataDir / 'allinteractome_mutations_onchains_ddg_new.txt'
    
    with io.open(inPath, "r", encoding="utf-8") as f, io.open(outPath, "w") as fout:
        fout.write( '\t'.join(['protein',
                               'partner',
                               'protein_pos',
                               'pdb_id',
                               'chain_id',
                               'chain_pos',
                               'chain_mutation',
                               'partner_chain']) + '\n' )
        next(f)
        for line in f:
            strsplit = list(map(str.strip, line.split('\t')))
            if len(strsplit) == 7:
                protein, partner, pr_pos, ch_pos, pdbid, ch_mut, ch_partner = strsplit
                chainID = ch_mut[1]
                fout.write( '\t'.join([protein,
                                       partner,
                                       pr_pos,
                                       pdbid,
                                       chainID,
                                       ch_pos,
                                       ch_mut,
                                       ch_partner]) + '\n' )
            elif len(strsplit) == 8:
                protein, partner, pr_pos, ch_pos, pdbid, ch_mut, ch_partner, status = strsplit
                chainID = ch_mut[1]
                if status is 'S':
                    fout.write( '\t'.join([protein,
                                           partner,
                                           pr_pos,
                                           pdbid,
                                           chainID,
                                           ch_pos,
                                           ch_mut,
                                           ch_partner]) + '\n' )
                else:
                    fout.write( '\t'.join([protein,
                                           partner,
                                           pr_pos,
                                           pdbid,
                                           chainID,
                                           ch_pos,
                                           ch_mut,
                                           ch_partner,
                                           status]) + '\n' )
            elif len(strsplit) == 9:
                protein, partner, pr_pos, ch_pos, pdbid, ch_mut, ch_partner, status, ddg = strsplit
                chainID = ch_mut[1]
                if status is 'S':
                    fout.write( '\t'.join([protein,
                                           partner,
                                           pr_pos,
                                           pdbid,
                                           chainID,
                                           ch_pos,
                                           ch_mut,
                                           ch_partner,
                                           ddg]) + '\n' )
            else:    
                fout.write('\n')

if __name__ == "__main__":
    main()