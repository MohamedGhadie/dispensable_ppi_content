import re
import time
import pickle
import pandas as pd
import numpy as np
from Bio import Seq, SeqIO
from Bio.PDB import *
from Bio.SeqUtils import seq1
from urllib import error

downloadStructures = True
suppressWarnings = False
MAX_STRUC = 100
strucList = []
structures = {}
chainResOrder = {}
validpos = seqmatch = numStruct = 0

class ChainSelect( Select ):
    
    def __init__(self, chain_letters):
        self.chain_letters = chain_letters
    
    def accept_chain(self, chain):
        return (chain.get_id() in self.chain_letters)

def clear_structures ():
    """Clear dict of PDB structures
    
    """
    global structures, strucList, numStruct
    structures.clear()
    strucList = []
    numStruct = 0

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
    global structures, strucList, numStruct
    if len( strucList ) == MAX_STRUC:
        del structures[ strucList.pop(0) ]
    structures[ pdbid ] = struc
    strucList.append( pdbid )
    numStruct += 1

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

def load_chainResOrder (inPath):
    """Load dict of sequence positions for structured residues in each PDB chain.
    
    """
    global chainResOrder
    print('\tloading chain residue order dictionary')
    with open(inPath, 'rb') as f:
        chainResOrder = pickle.load(f)

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

def download_pdb_files (inPath, outDir):
    """Old replicate function: Download PDB structure files for a list of pdb IDs.

    Args:
        inPath (str): file directory containing list of PDB IDs.
        outDir (str): directory to save PDB files to.

    """
    with open(inPath, 'r') as fin:
        pdbIDs = list(fin.read().split())
    print('\t%d PDB IDs to check' % len(pdbIDs))
    failed = 0
    pdbl = PDBList()
    for i, id in enumerate(pdbIDs):
        print('\t%d downloads attempted, %d downloads failed' % (i, failed))
        filename = outDir / ('pdb' + id + '.ent')
        if not filename.is_file():
            try:
                pdbl.retrieve_pdb_file(id, pdir = outDir)
            except error.URLError:
                print('\tPDB ID ' + id + ' file failed to download')
                failed += 1
    print('\t%d downloads attempted, %d downloads failed' % (i + 1, failed))

def residueID (residue):
    """Return residue ID

    Args:
        residue (Bio.PDB.Residue): residue whose ID is to be returned
    
    Returns:
        str: residue ID
    
    """
    het, pos, inCode = residue.get_id()
    return (het + '_' + str(pos) + '_' + inCode).strip()

def all_atoms (chain, residues):
    """Return all atoms combined for all residues in a chain

    Args:
        chain (Bio.PDB.Chain): chain whose atoms are to be returned
        residues (dict): residues per chain
    
    Returns:
        list: list of atoms
    
    """
    chainID = chain.get_id()
    allatoms = []
    atoms = atoms_per_residue(residues[chainID])
    for residue in residues[chainID]:
        resID = residueID(residue)
        allatoms.extend(atoms[resID])
    return allatoms

def atoms_per_residue (residues):
    """Return all atoms per residue for a list of PDB residues 

    Args:
        residues (list): list of residues
    
    Returns:
        dict: list of atoms for each residue
    
    """
    allatoms = dict()
    for residue in residues:
        resID = residueID(residue)
        atoms = residue.get_unpacked_list()
        if len(atoms) > 0:
            allatoms[resID] = atoms
    return allatoms

def residues_per_chain (model):
    """Return all residues per chain for all chains in a PDB structure 

    Args:
        model (Bio.PDB.Model): PDB model
    
    Returns:
        dict: list of all residues for each chain, including hetero-residues.
    
    """
    allresidues = dict()
    for chain in model.get_chains():
        residues = [res for res in chain.get_residues()]
        if len(residues) > 0:
            allresidues[chain.get_id()] = residues
    return allresidues

def ordered_residues_per_chain (pdbid, model):
    """Return all ordered residues per chain for all chains in a PDB structure 

    Args:
        pdbid (str): PDB ID for model parent structure.
        model (Bio.PDB.Model): PDB model
    
    Returns:
        dict: list of all ordered residues for each chain, including hetero-residues.
    
    """
    residues = residues_per_chain(model)
    for id in residues:
        fullid = pdbid + '_' + id
        if fullid in chainResOrder:
            residues[id] = residues[id][:len(chainResOrder[fullid])]
    return residues

def nonhet_residues_per_chain (model):
    """Return all non-hetero residues per chain for all chains in a PDB structure 

    Args:
        model (Bio.PDB.Model): PDB model
    
    Returns:
        dict: list of all non-hetero residues for each chain.
    
    """
    allresidues = dict()
    for chain in model.get_chains():
        residues = []
        for res in chain.get_residues():
            het, resnum, incode = res.get_id()
            if het == ' ':
                residues.append(res)
        if len(residues) > 0:
            allresidues[chain.get_id()] = residues
    return allresidues

def is_contact (res,
                other_residues,
                maxDist):
    
    for r in other_residues:
        if get_distance(res, r) < maxDist:
            return True
    return False

def get_distance (res1, res2):
    
    dist = [a1 - a2 for a1 in res1.get_unpacked_list() for a2 in res2.get_unpacked_list()]
    return min(dist)

def count_neighbors(res,
                    otherRes,
                    nbDist):
    
    nb = [get_distance(res, res2) < nbDist for res2 in otherRes]
    return sum(nb)

def get_num_neighbors (pdbDir,
                       pdbid,
                       resPos,
                       hostchainID,
                       chainIDs,
                       nbDist,
                       suppressWarnings):
    
    global structures
    if hostchainID in chainResOrder:
        try:
            resInd = chainResOrder[hostchainID].index(resPos)
        except ValueError:
            return np.nan
        hostchainID = hostchainID[ hostchainID.find('_') + 1 :]
        chainIDs = [id[ id.find('_') + 1 :] for id in chainIDs]
        if pdbid not in structures:
            structureFile = pdbDir / ('pdb' + pdbid + '.ent')
            if not structureFile.is_file():
                try:
                    PDBList().retrieve_pdb_file(pdbid, pdir = pdbDir)
                    print('Structure file for ' + pdbid + ' downloaded')
                except error.URLError:
                    print('\tPDB file for ' + pdbid + ' could not be downloaded')
                    return np.nan
            structures[pdbid] = PDBParser(QUIET=suppressWarnings).get_structure(pdbid, str(structureFile))
        if pdbid in structures:
            model = structures[pdbid][0]
            modelChains = [chain.get_id() for chain in model.get_chains()]
            matchedChains = sorted(set(chainIDs) & set(modelChains))
            if hostchainID in matchedChains:
                residues = ordered_residues_per_chain(pdbid, model)
                hostResidues = residues[hostchainID]
                res = hostResidues[resInd]
                otherRes = hostResidues.copy()
                otherRes.pop(resInd)
                for id in chainIDs:
                    if id != hostchainID:
                        otherRes.extend(residues[id])
                return count_neighbors(res, otherRes, nbDist)
    return np.nan

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
        resOrder1 (list): resseq positions of residues with known coordinates for
                            first chain.
        resOrder2 (list): resseq positions of residues with known coordinates for
                            second chain.
        maxDist (float): max distance for binding.
        suppressWarnings (boolean): True to suppress warnings while reading PDB structures.
    
    Returns:
        list: binding residues of the first chain

    """
    global structures, numStruct
    pdbID1, _ = chain_id1.split('_')
    pdbID2, _ = chain_id2.split('_')
    if pdbID1 == pdbID2:
        if numStruct > 99:
            clear_structures()
        if pdbID1 not in structures:
            structureFile = inDir / ('pdb' + pdbID1 + '.ent')
            if structureFile.is_file():
                structures[pdbID1] = PDBParser(QUIET=suppressWarnings).get_structure(pdbID1, str(structureFile))
                numStruct += 1
            else:
                print('Warning: structure file for chains %s and %s not found' % (chain_id1, chain_id2))
                return [], []
        return get_chain_interfaces(chain_id1,
                                    chain_id2,
                                    maxDist = maxDist)
    else:
        return [], []

def get_chain_interfaces (chain_id1,
                          chain_id2,
                          maxDist = 5):
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
    model = structures[pdbID][0]
    chain1 = model[id1]
    chain2 = model[id2]
    resOrder1 = chainResOrder[chain_id1]
    resOrder2 = chainResOrder[chain_id2]
    numRes1 = len(resOrder1)
    numRes2 = len(resOrder2)
    residues1 = [res for res in chain1.get_residues()][ : numRes1 ]
    residues2 = [res for res in chain2.get_residues()][ : numRes2 ]
    interfaceIndex1, interfaceIndex2 = get_interface_indices (residues1,
                                                              residues2,
                                                              maxDist = maxDist)
    interface1 = [resOrder1[i] for i in interfaceIndex1]
    interface2 = [resOrder2[i] for i in interfaceIndex2]
    return interface1, interface2

def get_interface_indices (residues1,
                           residues2,
                           maxDist = 5):
    """Identify residues in one PDB chain that are binding to another PDB chain.

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
        return sorted( set( interfaceIndex1 ) ), sorted( set( interfaceIndex2 ) )
    else:
        return [], []

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

def write_partial_structure (pdbid, chainIDs, inDir, outPath):
    
    struc = return_structure (pdbid, inDir)
    
    if struc:
        pdbio = PDBIO()
        pdbio.set_structure( structures[ pdbid ] )
        pdbio.save(str(outPath), select = ChainSelect(chainIDs))
