import os
import re
import time
import pickle
import warnings
from Bio import Seq, SeqIO
from Bio.PDB import *
from Bio.SeqUtils import seq1
from urllib import error
from itertools import compress

downloadStructures = True
suppressWarnings = False
MAX_STRUC = 100
strucList = []
structures = {}
chainSeq = {}
chainStrucResLabel = {}
chainResOrder = {}
resPosToID = {}

class ChainSelect( Select ):
    
    def __init__(self, chain_letters):
        self.chain_letters = chain_letters
    
    def accept_chain(self, chain):
        return (chain.get_id() in self.chain_letters)

def allow_pdb_downloads (download):
    
    global downloadStructures
    downloadStructures = download

def suppress_pdb_warnings (suppress):
    
    global suppressWarnings
    suppressWarnings = suppress

def load_pdbtools_chain_sequences (inPath):
    
    global chainSeq
    with open( inPath, 'rb' ) as f:
        chainSeq = pickle.load(f)

def load_pdbtools_chain_strucRes_labels (inPath):
    
    global chainStrucResLabel
    with open( inPath, 'rb' ) as f:
        chainStrucResLabel = pickle.load(f)

def load_chainResOrder ( inPath ):
    """Load dict of sequence positions for structured residues in each PDB chain.
    
    """
    global chainResOrder
    with open(inPath, 'rb') as f:
        chainResOrder = pickle.load(f)

def clear_structures ():
    """Clear dict of PDB structures
    
    """
    global structures, strucList
    structures.clear()
    strucList = []

def load_structure ( pdbid, structureFile):
    
    if structureFile.is_file():
        try:
            return PDBParser( QUIET=suppressWarnings ).get_structure( pdbid, str(structureFile) )
        except:
            return None
    else:
        print( '\t' + 'Structure file not found' )
        return None

def add_structure ( pdbid, struc ):
    """Add structure to dict.
    
    """
    global structures, strucList
    if len( strucList ) == MAX_STRUC:
        del structures[ strucList.pop(0) ]
    structures[ pdbid ] = struc
    strucList.append( pdbid )

def retrieve_structure ( pdbid, pdbDir):
    
    if pdbid not in structures:
        structureFile = pdbDir / ('pdb' + pdbid + '.ent')
        if downloadStructures and not structureFile.is_file():
            try:
                download_structure ( pdbid, pdbDir )
            except error.URLError:
                pass
        if structureFile.is_file():
            struc = load_structure (pdbid, structureFile)
            if struc:
                add_structure (pdbid, struc)

def return_structure (pdbid, pdbDir):
    
    retrieve_structure (pdbid, pdbDir)
    if pdbid in structures:
        return structures[ pdbid ]
    else:
        return None

def download_structures ( inPath, outDir ):
    """Download PDB structure files for a list of pdb IDs.

    Args:
        inPath (str): file directory containing list of PDB IDs.
        outDir (str): directory to save PDB files to.

    """
    with open(inPath, 'r') as fin:
        pdbIDs = list( fin.read().split() )
    print( '\t%d PDB IDs to check' % len( pdbIDs ) )
    failed = 0
    for i, id in enumerate( pdbIDs ):
        print( '\t' + '%d downloads attempted, %d downloads failed' % (i, failed) )
        filename = outDir / ('pdb' + id + '.ent')
        if not filename.is_file():
            try:
                download_structure ( id, outDir )
            except error.URLError:
                failed += 1
    print( '\t' + '%d downloads attempted, %d downloads failed' % (i + 1, failed) )

def download_structure ( pdbid, outDir, fileFormat = 'pdb' ):
    """Download PDB structure file.

    Args:
        pdbid (str): PDB ID of structure file to download.
        outDir (str): directory to save structure file to.

    """
    try:
        PDBList().retrieve_pdb_file( pdbid,
                                     pdir = outDir,
                                     file_format = fileFormat )
        filename = outDir / ('pdb' + pdbid + '.pdb')
        if ( fileFormat == 'pdb' ) and filename.is_file():
            base = os.path.splitext( filename )[ 0 ]
            os.rename( filename, base + '.ent' )
    except error.URLError:
        raise

def return_chain_sequence (chainID):
    
    if len(chainSeq):
        if chainID in chainSeq:
            return chainSeq[chainID]
    else:
        warnings.warn( "Chain sequence dictionary is empty" )
    return None

def return_chain_strucRes_label (chainID):
    
    if len(chainStrucResLabel):
        if chainID in chainStrucResLabel:
            return chainStrucResLabel[chainID]
    else:
        warnings.warn( "Chain structured residue label dictionary is empty" )
    return None

def residueID ( residue ):
    """Return residue ID

    Args:
        residue (Bio.PDB.Residue): residue whose ID is to be returned
    
    Returns:
        str: residue ID
    
    """
    het, pos, inCode = residue.get_id()
    return (het + '_' + str(pos) + '_' + inCode).strip()

def all_atoms ( chain, residues ):
    """Return all atoms combined for all residues in a chain

    Args:
        chain (Bio.PDB.Chain): chain whose atoms are to be returned
        residues (dict): residues per chain
    
    Returns:
        list: list of atoms
    
    """
    chainID = chain.get_id()
    allatoms = []
    atoms = atoms_per_residue( residues[ chainID ] )
    for residue in residues[ chainID ]:
        resID = residueID( residue )
        allatoms.extend( atoms[ resID ] )
    return allatoms

def atoms_per_residue ( residues ):
    """Return all atoms per residue for a list of PDB residues 

    Args:
        residues (list): list of residues
    
    Returns:
        dict: list of atoms for each residue
    
    """
    allatoms = {}
    for residue in residues:
        resID = residueID( residue )
        atoms = residue.get_unpacked_list()
        if len( atoms ) > 0:
            allatoms[ resID ] = atoms
    return allatoms
    
def residues_per_chain (model):
    """Return all residues per chain for all chains in a PDB structure 

    Args:
        model (Bio.PDB.Model): PDB model
    
    Returns:
        dict: list of all residues for each chain, including hetero-residues.
    
    """
    allresidues = {}
    for chain in model.get_chains():
        residues = [res for res in chain.get_residues()]
        if len(residues) > 0:
            allresidues[chain.get_id()] = residues
    return allresidues

def ordered_residues_per_chain (pdbid, model, pdbDir):
    """Return all ordered residues per chain for all chains in a PDB structure 

    Args:
        pdbid (str): PDB ID for model parent structure.
        model (Bio.PDB.Model): PDB model
    
    Returns:
        dict: list of all ordered residues for each chain, including hetero-residues.
    
    """
    resPerChain = {}
    for chain in model.get_chains():
        residues = ordered_chain_residues (pdbid, model, chain.get_id(), pdbDir)
        if residues:
            resPerChain[chain.get_id()] = residues
    return resPerChain

def nonhet_residues_per_chain (model):
    """Return all non-hetero residues per chain for all chains in a PDB structure 

    Args:
        model (Bio.PDB.Model): PDB model
    
    Returns:
        dict: list of all non-hetero residues for each chain.
    
    """
    allresidues = {}
    for chain in model.get_chains():
        residues = []
        for res in chain.get_residues():
            het, resnum, incode = res.get_id()
            if het == ' ':
                residues.append(res)
        if len(residues) > 0:
            allresidues[chain.get_id()] = residues
    return allresidues

def ordered_chain_residues (pdbid, model, chainID, pdbDir):
    """Return all residues of a chain that have a structure 

    Args:
        pdbid (str): PDB ID for model parent structure.
        model (Bio.PDB.Model): PDB model
    
    Returns:
        dict: list of all structured residues, including hetero-residues.
    
    """
    residues = chain_residues (model, chainID)
    pos = return_chain_res_IDsToPos (pdbid, chainID, pdbDir)
    return [res for res in residues if res.get_id() in pos]

def chain_residues (model, chainID):
    """Return all residues for a chain, including those with no structure 

    Args:
        pdbid (str): PDB ID for model parent structure.
        model (Bio.PDB.Model): PDB model
    
    Returns:
        dict: list of all residues, including hetero-residues.
    
    """
    chain = model[chainID]
    return [res for res in chain.get_residues()]

def is_contact (res, other_residues, maxDist):
    
    for r in other_residues:
        if get_distance(res, r) < maxDist:
            return True
    return False

def get_distance (res1, res2):
    
    dist = [a1 - a2 for a1 in res1.get_unpacked_list() for a2 in res2.get_unpacked_list()]
    return min(dist)

def get_interface_by_chainIDs (inDir,
                               chain_id1,
                               chain_id2,
                               maxDist = 5,
                               suppressWarnings = False):
    """Identify residues in one PDB chain that are binding to another PDB chain.

    Args:
        indDir (Path): directory containing PDB structure files.
        chain_id1 (str): ID of chain whose binding residues are to be determined.
        chain_id2 (str): ID of chain to which interface residues bind.
        maxDist (float): max distance for binding.
        suppressWarnings (boolean): True to suppress warnings while reading PDB structures.
    
    Returns:
        list: binding residues of the first chain

    """
    global structures
    pdbID1, _ = chain_id1.split('_')
    pdbID2, _ = chain_id2.split('_')
    if pdbID1 == pdbID2:
        return get_chain_interfaces (chain_id1, chain_id2, inDir, maxDist = maxDist)
    else:
        return [], []        

def get_chain_interfaces (chain_id1, chain_id2, pdbDir, maxDist = 5):
    """Identify residues in one PDB chain that are binding to another PDB chain.

    Args:
        residues1 (list): list of residues from which binding residues are identified.
        residues2 (list): list of residues to whom binding is checked.
        maxDist (float): max distance for binding. 
    
    Returns:
        list: binding residues.

    """
    pdbID, id1 = chain_id1.split('_')
    _, id2 = chain_id2.split('_')
    struc = return_structure (pdbID, pdbDir)
    if struc:
        model = struc[0]
        residues1 = ordered_chain_residues (pdbID, model, id1, pdbDir)
        residues2 = ordered_chain_residues (pdbID, model, id2, pdbDir)
        interfaceIndex1, interfaceIndex2 = get_interface_indices (residues1,
                                                                  residues2,
                                                                  maxDist = maxDist)
        interfaceResID1 = [residues1[i].get_id() for i in interfaceIndex1]
        interfaceResID2 = [residues2[i].get_id() for i in interfaceIndex2]
        chain1ResIDsToPos = return_chain_res_IDsToPos (pdbID, id1, pdbDir)
        chain2ResIDsToPos = return_chain_res_IDsToPos (pdbID, id2, pdbDir)
        interface1 = [chain1ResIDsToPos[id] for id in interfaceResID1 if id in chain1ResIDsToPos]
        interface2 = [chain2ResIDsToPos[id] for id in interfaceResID2 if id in chain2ResIDsToPos]
        return sorted(interface1), sorted(interface2)
    else:
        warnings.warn('Warning: structure file for chains %s and %s not found' % (chain_id1, chain_id2))
        return [], []

def get_interface_indices (residues1, residues2, maxDist = 5):
    """Identify residues at the interface between two groups of residues.

    Args:
        residues1 (list): first list of residues from which first binding interface is computed.
        residues2 (list): second list of residues from which second binding interface is computed.
        maxDist (float): max distance for binding. 
    
    Returns:
        list: binding residues.

    """
    numRes1 = len(residues1)
    numRes2 = len(residues2)
    resPairs = [(r1, r2) for r1 in range(numRes1) for r2 in range(numRes2)]
    distances = [get_distance(res1, res2) for res1 in residues1 for res2 in residues2]
    interfaceIndices = [pair for pair, d in zip(resPairs, distances) if d < maxDist]
    if len(interfaceIndices) > 0:
        interfaceIndex1, interfaceIndex2 = zip(*interfaceIndices)
        return sorted(set(interfaceIndex1)), sorted(set(interfaceIndex2))
    else:
        return [], []

def return_multichain_res_posToIDs (pdbid, chainResPos, pdbDir):
    
    IDs = []
    for chainID, pos in chainResPos:
        resID = return_chain_res_posToID (pdbid, chainID, pos, pdbDir)
        IDs.append( (chainID, pos, resID) )
    return IDs

def return_chain_res_posToID (pdbid, chainID, resPos, pdbDir):
    
    IDmaps = return_chain_res_posToIDs (pdbid, chainID, pdbDir)
    if resPos in IDmaps:
        return IDmaps[ resPos ]
    else:
        return None

def return_chain_res_IDToPos (pdbid, chainID, resID, pdbDir):
    
    posmaps = return_chain_res_IDsToPos (pdbid, chainID, pdbDir)
    if resID in posmaps:
        return posmaps[ resID ]
    else:
        return None
    
def return_chain_res_IDsToPos (pdbid, chainID, pdbDir):
    
    idmaps = return_chain_res_posToIDs (pdbid, chainID, pdbDir)
    posmaps = {}
    for pos, id in idmaps.items():
        posmaps[ id ] = pos
    return posmaps

def return_chain_res_posToIDs (pdbid, chainID, pdbDir):
    
    get_chain_res_posToIDs (pdbid, chainID, pdbDir)
    fullID = pdbid + '_' + chainID
    if fullID in resPosToID:
        return resPosToID[ fullID ]
    else:
        return {}

def get_chain_res_posToIDs (pdbid, chainID, pdbDir):
    
    global resPosToID
    fullID = pdbid + '_' + chainID
    if fullID not in resPosToID:
        struc = return_structure (pdbid, pdbDir)
        if struc:
            maps = allchain_res_pos_to_IDs (pdbid, struc)
            for id, posmaps in maps.items():
                resPosToID[ id ] = posmaps

def allchain_res_pos_to_IDs (pdbid, struc):
    
    allmaps = {}
    model = struc[ 0 ]
    for chain in model:
        idmaps = chain_res_pos_to_IDs (pdbid, chain)
        if idmaps:
            fullID = pdbid + '_' + chain.get_id()
            allmaps[ fullID ] = idmaps
    return allmaps

def chain_res_pos_to_IDs (pdbid, chain, minFracValidRes = 0.8):
    
    IDs = {}
    fullID = pdbid + '_' + chain.get_id()
    seq = return_chain_sequence (fullID)
    strucRes = return_chain_strucRes_label (fullID)
    if seq and strucRes:
        if len( seq ) == len( strucRes ):
            # construct chain sequence of structured residues based on SEQRES record 
            seqres = ''.join( compress(seq, [c == '-' for c in strucRes]) )
            seqresLen = len( seqres )
        
            # construct chain sequence of structured residues as read by BioPython
            coordresIDs = [ res.get_id() for res in chain.get_residues() ]
            coordres = [ seq1(res.get_resname(), undef_code='.') for res in chain.get_residues() ]
            coordresIDs = coordresIDs[ : seqresLen ]
            coordres = ''.join( coordres[ : seqresLen ] )
            
            # check if fraction of valid structured residue names read by Biopython 
            # meets cutoff
            if ( 1 - coordres.count('.') / len(coordres) ) >= minFracValidRes:
                # check if the sequence of structured residues read by Biopython 
                # is a subset of SEQRES 
                match = re.search(coordres, seqres)
                if match:
                    mStart = match.start()
                    for seqPos, label in enumerate(strucRes):
                        if (seqPos >= mStart) and (label == '-'):
                            if coordres:
                                IDs[seqPos + 1] = coordresIDs.pop(0)
    return IDs

def produce_chain_list (inPath, outPath):    
    """Produce list of PDB chain IDs from a fasta file of chain sequences.

    Args:
        inPath (str): file directory containing PDB chain fasta sequences.
        outPath (str): file directory to save chain IDs to.

    """
    chainList = []
    s = list(SeqIO.parse(str(inPath), 'fasta'))
    chainList = map(str.strip, [row.id for row in s])
    with open(outPath, 'w') as fout:
        for id in sorted(set(chainList)):
            fout.write(id + '\n')

def produce_pdb_ids (inPath, outPath):
    """Produce list of PDB IDs from a list of chain IDs.

    Args:
        inPath (str): file directory containing list of chain IDs.
        outPath (str): file directory to save PDB IDs to.

    """
    with open(inPath, 'r') as fin:
        chainIDs = list(fin.read().split())
    pdbIDs = [id[ : id.find('_')] for id in chainIDs if '_' in id]
    with open(outPath, 'w') as fout:
        for id in sorted(set(pdbIDs)):
            fout.write(id + '\n')

def produce_chainRes_posToID_dict (inPath,
                                   pdbDir,
                                   chainSeqFile,
                                   outPath,
                                   printInvalids = False,
                                   download = True,
                                   suppressWarnings = False):
    """Produce dict of resseq positions of residues with known coordinates for PDB chains
        with defined sequences and available structure files.

    Args:
        inPath (Path): file directory containing PDB's ss_dis.txt file.
        pdbDir (Path): file directory containing PDB structure files.
        chainSeqFile (Path): file directory containing chain sequence dict.
        outPath (Path): file directory to save dict to.

    """
    allow_pdb_downloads (download)
    suppress_pdb_warnings (suppressWarnings)
    clear_structures()
    
    fullcount = count = processed = matches = mismatches = lenMismatches = 0
    
    with open(chainSeqFile, 'rb') as f:
        chainSeq = pickle.load(f)
    
    if outPath.is_file():
        with open(geneMapFile, 'rb') as f:
            posToID = pickle.load(f)
    else:
        posToID = {}
    
    s = list( SeqIO.parse( str(inPath), 'fasta') )
    for i, row in enumerate(s):
        if ':disorder' in row.id:
            fullcount += 1
    print( '\t' + '%d chain sequence orders in total' % fullcount )
    
    s = list( SeqIO.parse( str(inPath), 'fasta') )
    for i, row in enumerate( s ):
        if ':disorder' in row.id:
            count += 1
            if not (count % 1000):
                print( '\t' + 'saving to file' )
                with open(outPath, 'wb') as fOut:
                    pickle.dump( posToID, fOut )
                print( '\t' + '%d chain sequence orders encountered' % count )
                print( '\t' + '%d chain sequence orders processed' % processed )
                print( '\t' + '%d chain length mismatches with ordered sequence' % lenMismatches )
                print( '\t' + 'Among those matching in length with structure file available:' )
                print( '\t' + '%d chain sequences match with sequence of structured residues' % matches )
                print( '\t' + '%d chain sequences mismatch with sequence of tructured residues' % mismatches )
            pdbid, chainID, _ = list( map( str.strip, row.id.split(':') ) )
            pdbid = pdbid.lower()
            fullID = pdbid + '_' + chainID
            seq = return_chain_sequence (fullID)
            if seq and (fullID not in posToID):
                if len(seq) == len(row.seq):
                    struc = return_structure (pdbid, pdbDir)
                    if struc:
                        processed += 1
                        chain = struc[ 0 ][ chainID ]
                        
                        # construct chain sequence of structured residues based on SEQRES record 
                        seqres = ''.join( [ r for k, r in enumerate(seq) if row.seq[ k ] == '-' ] )
                    
                        # construct chain sequence of structured residues as read by BioPython
                        coordres = []
                        coordresIDs = []
                        for r in chain.get_residues():
                            coordresIDs.append( r.get_id() )
                            coordres.append( seq1(r.get_resname(), undef_code='.') )
                        coordres = ''.join( coordres[ : len( seqres ) ] )
                    
                        # check if the sequence of structured residues read by Biopython 
                        # matches with SEQRES 
                        match = [ m.start() for m in re.finditer('(?=%s)' % coordres, seqres) ]
                        if 0 in match:
                            matches += 1
                            IDs = []
                            for c in row.seq:
                                IDs.append( coordresIDs.pop(0) if c == '-' else () )
                            posToID[ fullID ] = IDs
                        else:
                            mismatches += 1
                            if printInvalids:
                                print( '\t' + 'Sequence of structured residues does not match SEQRES:' )
                                print(fullID)
                                print( '\t' + 'Full chain sequence (len = %d):' % len(seq) )
                                print(seq)
                                print( '\t' + 'Structured residues (len = %d):' % len(seqres) )
                                print(seqres)
                                print( '\t' + 'structured residues as read (len = %d)' % len(coordres) )
                                print(coordres)
                else:
                    lenMismatches += 1
    
    with open(outPath, 'wb') as fOut:
        pickle.dump( posToID, fOut )

def write_partial_structure (pdbid, chainIDs, inDir, outPath):
    
    struc = return_structure (pdbid, inDir)
    
    if struc:
        pdbio = PDBIO()
        pdbio.set_structure(struc)
        pdbio.save(str(outPath), select = ChainSelect(chainIDs))

def produce_seqres_order_dict (inPath, pdbDir, chainSeqFile, outPath):
    """Produce dict of resseq positions of residues with known coordinates for PDB chains
        with defined sequences and available structure files.

    Args:
        inPath (Path): file directory containing PDB's ss_dis.txt file.
        pdbDir (Path): file directory containing PDB structure files.
        chainSeqFile (Path): file directory containing chain sequence dict.
        outPath (Path): file directory to save dict to.

    """
    count = matches = mismatches = 0
    with open(chainSeqFile, 'rb') as f:
        chainSeq = pickle.load(f)
    seqres_order = {}
    s = list(SeqIO.parse(str(inPath), 'fasta'))
    for i, row in enumerate(s):
        if ':disorder' in row.id:
            count += 1
            if not (count % 1000):
                print(count)
            if not (count % 100000):
                time.sleep(180)
            pdbid, chainID, _ = list( map( str.strip, row.id.split(':') ) )
            pdbid = pdbid.lower()
            fullID = pdbid + '_' + chainID
            if fullID in chainSeq:
                seq = chainSeq[fullID]
                if len(seq) == len(row.seq):
                    structureFile = pdbDir / ('pdb' + pdbid + '.ent')
                    if structureFile.is_file():
                        struc = PDBParser(QUIET=True).get_structure(pdbid, str(structureFile))
                        model = struc[0]
                        # construct chain sequence of structured residues based on SEQRES record 
                        seqres = ''.join( [ r for k, r in enumerate(seq) if row.seq[ k ] == '-' ] )
                        # construct chain sequence of structured residues as read by BioPython
                        coordres = []
                        for r in model[ chainID ].get_residues():
                            coordres.append( seq1(r.get_resname(), undef_code='.') )
                        coordres = ''.join( coordres[ : len( seqres ) ] )
                        # check if sequence of structured residues matches between SEQRES and BioPython
                        # or at least SEQRES sequence starts with sequence read by BioPython 
                        match = [ m.start() for m in re.finditer('(?=%s)' % coordres, seqres) ]
                        if 0 in match:
                            matches += 1
                            seqres_order[fullID] = [k + 1 for k, ch in enumerate(row.seq) if ch == '-']
                        else:
                            mismatches += 1
                            print('%d mismatches' % mismatches)
                            print(fullID)
                            print('Full chain sequence (len = %d):' % len(seq))
                            print(seq)
                            print('Structured residues (len = %d):' % len(seqres))
                            print(seqres)
                            print('structured residues as read (len = %d)' % len(coordres))
                            print(coordres)
    with open(outPath, 'wb') as fOut:
        pickle.dump(seqres_order, fOut)
