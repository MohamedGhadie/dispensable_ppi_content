
def read_unprocessed_ddg_mutations (inPath):
    
    mutations = {}
    done = set()
    with io.open(inPath, "r", encoding="utf-8") as f:
        next(f)
        for line in f:
            strsplit = list( map ( str.strip, line.split('\t') ) )
            if len(strsplit) >= 8:
                protein, partner, pr_pos, pdbid, chainID, ch_pos, mut, ch_partner = strsplit[:8]
                if type is 'binding':
                    mutation = '-'.join( [protein, partner, pr_pos, mut[-1]] )
                elif type is 'folding':
                    mutation = '-'.join( [protein, pr_pos, mut[-1]] )
                if mutation not in done:
                    if len(strsplit) == 8:
                        done.add(mutation)
                        if type is 'binding':
                            struc = (pdbid,) + tuple(sorted([chainID, ch_partner]))
                        elif type is 'folding':
                            struc = pdbid, chainID
                        if struc in mutations:
                            mutations[struc].add(mut)
                        else:
                            mutations[struc] = {mut}
                    elif strsplit[8] is not 'X':
                        done.add(mutation)
    return {k:list(v) for k, v in mutations.items()}

def produce_bindprofx_jobs (mutations,
                            pdbDir,
                            outDir,
                            write_hpc_jobfiles = True,
                            nodes = 1,
                            ppn = 1,
                            pmem = 7700,
                            walltime = '1:00:00:00',
                            rapid = None,
                            username = '',
                            hpcCommands = None,
                            serverDataDir = '../data'):
    
    clear_structures()
    
    dataDir = outDir / 'data'
    jobDir = outDir / 'jobs'
    if not dataDir.exists():
        os.makedirs(dataDir)
    if not jobDir.exists():
        os.makedirs(jobDir)
    
    for struc, mutList in mutations.items():
        strucid = '_'.join(struc)
        strucDir = dataDir / strucid
        if not strucDir.exists():
            os.makedirs(strucDir)
        mutListFile = strucDir / 'mutList.txt'
        mutList = [ '%s;' % mutList.pop(0) ] + ['\n%s;' % mut for mut in mutList]
        with io.open(mutListFile, "w") as fout:
            for mut in mutList:
                fout.write(mut)
        pdbid, chainID1, chainID2 = struc[:3]
        write_partial_structure (pdbid,
                                 [chainID1, chainID2],
                                 pdbDir,
                                 strucDir / 'complex.pdb')
        
        if write_hpc_jobfiles:
            commands = [ '../bin/get_final_score.py %s/%s' % (serverDataDir, strucid) ]
            if hpcCommands:
                commands = hpcCommands + commands
            write_hpc_job (jobDir / (strucid + '_job.txt'),
                           nodes = nodes,
                           ppn = ppn,
                           pmem = pmem,
                           walltime = walltime,
                           outputfile = '%s/%s/outputfile' % (serverDataDir, strucid),
                           errorfile = '%s/%s/errorfile' % (serverDataDir, strucid),
                           rapid = rapid,
                           jobid = username + '_bindprofx_' + strucid,
                           commands = commands)

def read_bindprofx_results (inDir):
    
    processed = {}
    unprocessed = {}
    strucDir = os.listdir(inDir)
    strucDir = [dir for dir in strucDir if os.path.isdir(inDir / dir)]
    for strucID in strucDir:
        struc = tuple(strucID.split('_'))
        if re.match(r'\D\S\d+\D', struc[-1]):
            struc = struc[:-1]
        
        results = []
        resultFile = inDir / strucID / 'result.txt'
        if resultFile.is_file():
            with io.open(resultFile, "r") as f:
                results = list( map(str.strip, f.read().split(';')) )
            results.remove('')
            for result in results:
                ddg, mut = result.split()
                processed[struc + (mut,)] = float(ddg)
        
        if not results:
            mutListFile = inDir / strucID / 'mutList.txt'
            with open(mutListFile, 'r') as f:
                mutList = list( map(str.strip, f.read().split(';')) )
            mutList.remove('')
            if len(mutList) == 1:
                processed[struc + (mutList.pop(),)] = 'X'
            elif len(mutList) > 1:
                errorFile = inDir / strucID / 'errorfile'
                with open(errorFile, 'r') as f:
                    error = f.read()
                if (('FileNotFoundError: [Errno 2] No such file or directory:' in error) and
                    ('align.out' in error)):
                    for mut in mutList:
                        processed[struc + (mut,)] = 'X'
                else:
                    for mut in mutList:
                        unprocessed[struc + (mut,)] = [mut]
    
    return processed, unprocessed

def write_foldx_config (outPath,
                        command,
                        pdb_dir = None,
                        output_dir = None,
                        pdb_file = None,
                        mutant_file = None,
                        temp = 298,
                        ph = 7,
                        ionStrength = 0.05,
                        water = '-IGNORE',
                        vdwDesign = 2):
    
    with io.open(outPath, "w") as fout:
        fout.write('command=%s' % command)
        if pdb_dir:
            fout.write('\n' + 'pdb-dir=%s' % pdb_dir)
        if output_dir:
            fout.write('\n' + 'output-dir=%s' % output_dir)
        if pdb_file:
            fout.write('\n' + 'pdb=%s' % pdb_file)
        if mutant_file:
            fout.write('\n' + 'mutant-file=%s' % mutant_file)
        fout.write('\n' + 'temperature=%.1f' % temp)
        fout.write('\n' + 'pH=%.1f' % ph)
        fout.write('\n' + 'ionStrength=%f' % ionStrength)
        fout.write('\n' + 'water=%s' % water)
        fout.write('\n' + 'vdwDesign=%d' % vdwDesign)

def read_foldx_results (inDir):
    
    processed = {}
    unprocessed = {}
    strucDir = os.listdir(inDir)
    strucDir = [dir for dir in strucDir if os.path.isdir(inDir / dir)]
    for strucID in strucDir:
        struc = tuple(strucID.split('_'))
        if re.match(r'\D\S\d+\D', struc[-1]):
            struc = struc[:-1]
        mutListFile = inDir / strucID / 'individual_list.txt'
        with io.open(mutListFile, "r") as f:
            mutList = list( map(str.strip, f.read().split(';')) )
        mutList.remove('')
        
        resultFile = None
        for filename in os.listdir(inDir / strucID):
            if filename.startswith('Average'):
                resultFile = filename
        
        results = []
        if resultFile:
            resultFile = inDir / strucID / resultFile
            headers = ["Pdb", "SD", "total energy", "Backbone Hbond", "Sidechain Hbond",
                       "Van der Waals", "Electrostatics", "Solvation Polar", "Solvation Hydrophobic", 
                       "Van der Waals clashes", "entropy sidechain", "entropy mainchain", 
                       "sloop_entropy", "mloop_entropy", "cis_bond", "torsional clash", "backbone clash", 
                       "helix dipole", "water bridge", "disulfide", "electrostatic kon", 
                       "partial covalent bonds", "energy Ionisation", "Entropy Complex"]
            with io.open(resultFile, "r") as f:
                for line in f:
                    linesplit = list( map ( str.strip, line.split('\t') ) )
                    if linesplit == headers:
                        break
                for line in f:
                    linesplit = list( map ( str.strip, line.split('\t') ) )
                    if len(linesplit) > 3:
                        results.append(linesplit[2])
                if len(results) == len(mutList):
                    for mut, ddg in zip(mutList, results):
                        processed[struc + (mut,)] = float(ddg)
        
        if not results:
            if len(mutList) == 1:
                processed[struc + (mutList.pop(),)] = 'X'
            elif len(mutList) > 1:
                for mut in mutList:
                    unprocessed[struc + (mut,)] = [mut]
    
    return processed, unprocessed

def write_hpc_job (outPath,
                   nodes = 1,
                   ppn = 1,
                   pmem = 7700,
                   walltime = '1:00:00:00',
                   outputfile = 'outputfile',
                   errorfile = 'errorfile',
                   rapid = None,
                   jobid = None,
                   commands = None):
    
    with io.open(outPath, "w") as fout:
        fout.write('#!/bin/bash')
        fout.write('\n' + '#PBS -l nodes=%d:ppn=%d,pmem=%dm,walltime=%s' % (nodes, ppn, pmem, walltime))
        if rapid:
            fout.write('\n' + '#PBS -A %s' % rapid)
        fout.write('\n' + '#PBS -o %s' % outputfile)
        fout.write('\n' + '#PBS -e %s' % errorfile)
        if jobid:
            fout.write('\n' + '#PBS -N %s' % jobid)
        fout.write('\n\n' + 'cd $PBS_O_WORKDIR')
        if commands:
            for cmd in commands:
                fout.write('\n' + cmd)

# def read_interface_mutation_ddg (inPath):
#     
#     """Read change in binding free energy upon mutation.
# 
#     Args:
#         inPath (str): file directory containing free energy change. 
# 
#     """
#     ddgDict = {}
#     with io.open(inPath, "r", encoding="utf-8") as f:
#         next(f)
#         for line in f:
#             linesplit = list( map ( str.strip, line.split('\t') ) )
#             if len(linesplit) == 8:
#                 (protein, partner, protein_pos, chain_pos, pdb_id, 
#                  chain_mut, partner_chain, ddg) = linesplit 
#                 if ddg != 'X':
#                     k = protein, partner, int(protein_pos)
#                     if k not in ddgDict:
#                         ddgDict[k] = pdb_id, chain_mut, partner_chain, float(ddg)
#     return ddgDict
# 
# def read_mutation_ddg (inPath):
#     
#     """Read change in free energy of structures upon mutation.
# 
#     Args:
#         inPath (str): file directory containing free energy change. 
# 
#     """
#     ddgDict = {}
#     with io.open(inPath, "r", encoding="utf-8") as f:
#         next(f)
#         for line in f:
#             linesplit = list( map ( str.strip, line.split('\t') ) )
#             if len(linesplit) == 8:
#                 gene, protein, protein_pos, pdbid, chainID, chain_pos, mut, ddg = linesplit
#                 if ddg != 'X':
#                     k = protein, int(protein_pos), '_'.join([pdbid, chainID]), mut[-1] 
#                     ddgDict[k] = float(ddg)
#     return ddgDict

def read_protein_mutation_ddg (inPath):
    
    """Read change in free energy of structures upon mutation.

    Args:
        inPath (str): file directory containing free energy change. 

    """
    ddgDict = {}
    with io.open(inPath, "r", encoding="utf-8") as f:
        next(f)
        for line in f:
            linesplit = list( map ( str.strip, line.split('\t') ) )
            if len(linesplit) > 8:
                if linesplit[8] is not 'X':
                    protein, partner, pr_pos, pdbid, chainID, ch_pos, ch_mut, ch_partner, ddg = strsplit
                    if type is 'binding':
                        k = protein, partner, int(pr_pos), ch_mut[-1]
                        val = pdb_id, chainID, ch_partner, float(ddg)
                    elif type is 'folding':
                        k = protein, int(pr_pos), ch_mut[-1]
                        val = pdbid, chainID, float(ddg) 
                    if k not in ddgDict:
                        ddgDict[k] = val
    return ddgDict

def read_chain_mutation_ddg (inPath, type):
    
    """Read change in free energy of structures upon mutation.

    Args:
        inPath (str): file directory containing free energy change. 

    """
    ddgDict = {}
    with io.open(inPath, "r", encoding="utf-8") as f:
        next(f)
        for line in f:
            strsplit = list( map ( str.strip, line.split('\t') ) )
            if len(strsplit) > 8:
                protein, partner, pr_pos, pdbid, chainID, ch_pos, ch_mut, ch_partner, ddg = strsplit[:9]
                if type is 'binding':
                    k = (pdbid,) + tuple(sorted([chainID, ch_partner])) + (ch_mut,)
                elif type is 'folding':
                    k = pdbid, chainID, ch_mut
                ddgDict[k] = ddg
    return ddgDict

def copy_mutation_ddg (inPath1, inPath2, outPath, type):
    
    """Copy mutation change in free energy from one file to another.

    Args:
        inPath1 (str): file directory containing change in binding free energy.
        inPath2 (str): file directory containing mutations with unknown change in binding free energy.
        outPath (str): file directory to output mutations with change in binding free energy.

    """
    ddg = read_chain_mutation_ddg (inPath1, type)
    write_mutation_ddg_tofile (ddg, inPath2, outPath, type)

def write_mutation_ddg_tofile (ddg, inPath, outPath, type):
    
    with io.open(inPath, "r", encoding="utf-8") as f, io.open(outPath, "w") as fout:
        for line in f:
            strsplit = list( map ( str.strip, line.split('\t') ) )
            if len(strsplit) == 8:
                protein, partner, pr_pos, pdbid, chainID, ch_pos, ch_mut, ch_partner = strsplit
                if type is 'binding':
                    k = (pdbid,) + tuple(sorted([chainID, ch_partner])) + (ch_mut,)
                elif type is 'folding':
                    k = pdbid, chainID, ch_mut
                if k in ddg:
                    strsplit.append(str(ddg[k]))
            fout.write( '\t'.join(map(str, strsplit)) +  '\n' )

