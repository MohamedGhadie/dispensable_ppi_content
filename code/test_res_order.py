import io
import pickle
from pathlib import Path
from text_tools import read_list_table

def main():
    
    inDir = Path('../data/processed')
    
    inPath = inDir / 'chain_seqres_order_old.pkl'
    with open(inPath, 'rb') as f:
        chainResOrder_old = pickle.load(f)
    
    inPath = inDir / 'chain_seqres_order.pkl'
    with open(inPath, 'rb') as f:
        chainResOrder = pickle.load(f)
    
    pdbIDs = []
    for k in chainResOrder:
        pdbIDs.append( k.split('_')[ 0 ] )
    pdbIDs = sorted( set( pdbIDs ) )
    
    found = notfound = match = mismatch = low = newpdb = 0
    chainsToRemove = []
    for k in chainResOrder_old:
        if k in chainResOrder:
            found += 1
            if chainResOrder_old[k] == chainResOrder[k]:
                match += 1
            else:
                chainsToRemove.append( k )
                mismatch += 1
        else:
            chainsToRemove.append( k )
            pdbid, id = k.split('_')
            if id.islower():
                low += 1
            if pdbid not in pdbIDs:
                newpdb += 1
            notfound += 1
    
    print( '\t' + 'Found = %d' % found )
    print( '\t' + 'Not found = %d' % notfound )
    print( '\t' + 'Match = %d' % match )
    print( '\t' + 'Mismatch = %d' % mismatch )
    print( '\t' + 'Lower case among not found = %d' % low )
    print( '\t' + 'New PDBs among those not found = %d' % newpdb )
    print( '\t' + 'Chains to be removed: %s' % ','.join( chainsToRemove ) )
    
    return
    
    pdb_interfaces_old = read_chain_interfaces ( inDir / 'pdb_interfaces_oldResOrder.txt' )
    pdb_interfaces = dict()
    for k in pdb_interfaces_old:
        chain1, chain2 = k.split('-')
        if (chain1 not in chainsToRemove) and (chain2 not in chainsToRemove):
            pdb_interfaces[ k ] = pdb_interfaces_old[ k ]
    write_chain_interfaces (pdb_interfaces, inDir / 'pdb_interfaces.txt' )
    
    print( '\t' 'Checking new interface file' )
    pdb_interfaces = read_chain_interfaces ( inDir / 'pdb_interfaces.txt' )
    errorCount = 0
    for k in pdb_interfaces:
        chain1, chain2 = k.split('-')
        if (chain1 in chainsToRemove) or (chain2 in chainsToRemove):
            print( '\t' + 'Chain error: %s' % k )
            errorCount += 1
    if errorCount == 0:
        print( '\t' + 'All chains valid' )

def read_chain_interfaces (inPath):
    
    interfaces = dict()
    if inPath.is_file():
        interfaces_df = read_list_table(inPath, "Chain1_interface", int, '\t')
        for _, row in interfaces_df.iterrows():
            chainKey = row.Chain_1 + '-' + row.Chain_2
            if -1 in row.Chain1_interface:
                interfaces[ chainKey ] = []
            else:
                interfaces[ chainKey ] = row.Chain1_interface
    return interfaces

def write_chain_interfaces (interfaces, outPath):
    
    with io.open(outPath, "w") as fout:
        fout.write('\t'.join(["Chain_1", "Chain_2", "Chain1_interface"]) + '\n')
        for chainKey in sorted(interfaces.keys()):
            chain1, chain2 = chainKey.split('-')
            if len(interfaces[chainKey]) > 0:
                interface = ','.join([str(elm) for elm in interfaces[chainKey]])
            else:
                interface = '-1'
            fout.write('\t'.join([chain1, chain2, interface]) + '\n')

if __name__ == "__main__":
    main()