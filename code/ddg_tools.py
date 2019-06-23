#----------------------------------------------------------------------------------------
# Modules for calculating and processing change in free energy (∆∆G) data.
#----------------------------------------------------------------------------------------

import os
import io
from pathlib import Path
from text_tools import write_guillimin_job, write_beluga_job
from pdb_tools import clear_structures, write_partial_structure

def read_unprocessed_ddg_mutations (inPath, type = 'binding'):
    """Read PDB chain mutations with missing ∆∆G values from file.

    Args:
        inPath (Path): path to file containing mutations.
        type (str): type of ∆∆G values, 'binding' for interface ∆∆G, 'folding' for protein folding ∆∆G.

    Returns:
        dict

    """
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

def produce_bindprofx_and_guillimin_jobs (mutations,
                                          pdbDir,
                                          outDir,
                                          nodes = 1,
                                          ppn = 1,
                                          pmem = 7700,
                                          walltime = '1:00:00:00',
                                          rapid = None,
                                          username = '',
                                          extraCommands = None,
                                          serverDataDir = '../data'):
    """Produce data files for BindProfX ∆∆G calculations as well as job files for Guillimin 
        server. See http://www.hpc.mcgill.ca

    Args:
        mutations (dict): mutations associated with each structural model.
        pdbDir (Path): file directory containing PDB structures.
        outDir (Path): file directory to save BindProfX data files and Guillimin job files to.
        nodes (numeric): number of server nodes to be allocated.
        ppn (numeric): total number of CPU cores to be allocated.
        pmem (numeric): default random access memory (RAM) in MB to be reserved per core.
        walltime (str): maximum time allowed for job to run on server in format days:hr:min:sec.
        rapid (str): resource allocation project identifier (RAPid) for HPC user.
        username (str): user name to associate with HPC server job.
        extraCommands (list): additional non-bindprofx commands to be written to job file.
        serverDataDir (Path): data directory used by HPC server relative to job directory.

    """
    produce_bindprofx_jobs (mutations, pdbDir, outDir / 'data')
    produce_guillimin_bindprofx_jobs (mutations,
                                      outDir / 'jobs',
                                      nodes = nodes,
                                      ppn = ppn,
                                      pmem = pmem,
                                      walltime = walltime,
                                      rapid = rapid,
                                      username = username,
                                      extraCommands = extraCommands,
                                      serverDataDir = serverDataDir,)

def produce_bindprofx_and_beluga_jobs (mutations,
                                       pdbDir,
                                       outDir,
                                       account = 'ctb-yxia',
                                       walltime = '1-00',
                                       ntasks = 1,
                                       nodes = 1,
                                       ntasks_per_node = 1,
                                       cpus_per_task = 1,
                                       mem = None,
                                       mem_per_cpu = '4000M',
                                       outputfile = '%x-%j.out',
                                       errorfile = None,
                                       username = '',
                                       extraCommands = None,
                                       serverDataDir = '../data'):
    """Produce data files for BindProfX ∆∆G calculations as well as job files for Beluga server.
        See https://docs.computecanada.ca/wiki/Béluga/en

    Args:
        mutations (dict): mutations associated with each structural model.
        pdbDir (Path): file directory containing PDB structures.
        outDir (Path): file directory to save jobs to.
        account (str): project account name.
        walltime (str): maximum time allowed for job to run.
        ntasks (numeric): number of processes to be allocated.
        nodes (numeric): number of server nodes to be allocated.
        ntasks_per_node (numeric): number of processes to be allocated per node.
        cpus_per_task (numeric): number of nodes to be allocated per process.
        mem (str): memory per node.
        mem_per_cpu (str): memory per core.
        outputfile (Path): path to file where standard output is written.
        errorfile (Path): path to file where runtime error is written.
        username (str): user name to associate with Beluga server job.
        extraCommands (list): additional non-bindprofx commands to be written to job file.
        serverDataDir (Path): data directory used by Beluga server relative to job directory.

    """
    produce_bindprofx_jobs (mutations, pdbDir, outDir / 'data')
    produce_beluga_bindprofx_jobs (mutations,
                                   outDir / 'jobs',
                                   account = account,
                                   walltime = walltime,
                                   ntasks = ntasks,
                                   nodes = nodes,
                                   ntasks_per_node = ntasks_per_node,
                                   cpus_per_task = cpus_per_task,
                                   mem = mem,
                                   mem_per_cpu = mem_per_cpu,
                                   outputfile = outputfile,
                                   errorfile = errorfile,
                                   username = username,
                                   extraCommands = extraCommands,
                                   serverDataDir = serverDataDir)

def produce_bindprofx_jobs (mutations, pdbDir, outDir):
    """Produce jobs to be submitted to BindProfX server for ∆∆G calculations.

    Args:
        mutations (dict): mutations associated with each structural model.
        pdbDir (Path): file directory containing PDB structures.
        outDir (Path): file directory to save BindProfX jobs to.

    """
    clear_structures()
    if not outDir.exists():
        os.makedirs(outDir)
    
    for struc, mutList in mutations.items():
        strucid = '_'.join(struc)
        strucDir = outDir / strucid
        if not strucDir.exists():
            os.makedirs(strucDir)
        mutListFile = strucDir / 'mutList.txt'
        mutList = ['%s;' % mutList.pop(0)] + ['\n%s;' % mut for mut in mutList]
        with io.open(mutListFile, "w") as fout:
            for mut in mutList:
                fout.write(mut)
        pdbid, chainID1, chainID2 = struc[:3]
        write_partial_structure (pdbid,
                                 [chainID1, chainID2],
                                 pdbDir,
                                 strucDir / 'complex.pdb')

def produce_guillimin_bindprofx_jobs (mutations,
                                      outDir,
                                      nodes = 1,
                                      ppn = 1,
                                      pmem = 7700,
                                      walltime = '1:00:00:00',
                                      rapid = None,
                                      username = '',
                                      extraCommands = None,
                                      serverDataDir = '../data'):
    """Produce Guillimin server job files specific to BindProfX ∆∆G calculations.
        See http://www.hpc.mcgill.ca

    Args:
        mutations (dict): mutations associated with each structural model.
        outDir (Path): file directory to save jobs to.
        nodes (numeric): number of server nodes to be allocated.
        ppn (numeric): total number of CPU cores to be allocated.
        pmem (numeric): default random access memory (RAM) in MB to be reserved per core.
        walltime (str): maximum time allowed for job to run on server in format days:hr:min:sec.
        rapid (str): resource allocation project identifier (RAPid) for HPC user.
        username (str): user name to associate with HPC server job.
        extraCommands (list): additional non-bindprofx commands to be written to job file.
        serverDataDir (Path): data directory used by HPC server relative to job directory.

    """
    if not outDir.exists():
        os.makedirs(outDir)
    
    for struc, mutList in mutations.items():
        strucid = '_'.join(struc)
        commands = ['../bin/get_final_score.py %s/%s' % (serverDataDir, strucid)]
        if extraCommands:
            commands = extraCommands + commands
        write_guillimin_job (outDir / (strucid + '_job.txt'),
                             nodes = nodes,
                             ppn = ppn,
                             pmem = pmem,
                             walltime = walltime,
                             outputfile = '%s/%s/outputfile' % (serverDataDir, strucid),
                             errorfile = '%s/%s/errorfile' % (serverDataDir, strucid),
                             rapid = rapid,
                             jobid = username + '_bindprofx_' + strucid,
                             commands = commands)

def produce_beluga_bindprofx_jobs (mutations,
                                   outDir,
                                   account = 'ctb-yxia',
                                   walltime = '1-00',
                                   ntasks = 1,
                                   nodes = 1,
                                   ntasks_per_node = 1,
                                   cpus_per_task = 1,
                                   mem = None,
                                   mem_per_cpu = '4000M',
                                   outputfile = '%x-%j.out',
                                   errorfile = None,
                                   username = '',
                                   extraCommands = None,
                                   serverDataDir = '../data'):
    """Produce Beluga server job files specific to BindProfX ∆∆G calculations.
        See https://docs.computecanada.ca/wiki/Béluga/en

    Args:
        mutations (dict): mutations associated with each structural model.
        outDir (Path): file directory to save jobs to.
        account (str): project account name.
        walltime (str): maximum time allowed for job to run.
        ntasks (numeric): number of processes to be allocated.
        nodes (numeric): number of server nodes to be allocated.
        ntasks_per_node (numeric): number of processes to be allocated per node.
        cpus_per_task (numeric): number of nodes to be allocated per process.
        mem (str): memory per node.
        mem_per_cpu (str): memory per core.
        outputfile (Path): path to file where standard output is written.
        errorfile (Path): path to file where runtime error is written.
        username (str): user name to associate with Beluga server job.
        extraCommands (list): additional non-bindprofx commands to be written to job file.
        serverDataDir (Path): data directory used by Beluga server relative to job directory.

    """
    if not outDir.exists():
        os.makedirs(outDir)
    
    for struc, mutList in mutations.items():
        strucid = '_'.join(struc)
        commands = ['../bin/get_final_score.py %s/%s' % (serverDataDir, strucid)]
        if extraCommands:
            commands = extraCommands + commands
        write_beluga_job (outDir / (strucid + '_job.txt'),
                          account = account,
                          walltime = walltime,
                          ntasks = ntasks,
                          nodes = nodes,
                          ntasks_per_node = ntasks_per_node,
                          cpus_per_task = cpus_per_task,
                          mem = mem,
                          mem_per_cpu = mem_per_cpu,
                          outputfile = outputfile,
                          errorfile = errorfile,
                          jobname = username + '_bindprofx_' + strucid,
                          commands = commands)

def read_bindprofx_results (inDir):
    """Read mutation ∆∆G results produced by BindProfX.

    Args:
        inDir (Path): file directory containing BindProfX results.

    Returns:
        dict, dict: processed and unprocessed mutations.

    """
    processed, unprocessed = {}, {}
    strucDir = os.listdir(inDir)
    strucDir = [dir for dir in strucDir if os.path.isdir(inDir / dir)]
    for strucID in strucDir:
        struc = tuple(strucID.split('_'))
        if len(struc) > 3:
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

def produce_foldx_and_guillimin_jobs (mutations,
                                      pdbDir,
                                      outDir,
                                      foldxParam = None,
                                      nodes = 1,
                                      ppn = 1,
                                      pmem = 7700,
                                      walltime = '1:00:00:00',
                                      rapid = None,
                                      username = '',
                                      extraCommands = None,
                                      serverDataDir = '../data'):
    """Produce data files for FoldX ∆∆G calculations as well as job files for Guillimin 
    server. See http://www.hpc.mcgill.ca

    Args:
        mutations (dict): mutations associated with each structural model.
        pdbDir (Path): file directory containing PDB structures.
        outDir (Path): file directory to save FoldX data files and Guillimin job files to.
        foldxParam (dict): FoldX parameters used for each mutation job, otherwise default.
        nodes (numeric): number of server nodes to be allocated.
        ppn (numeric): total number of CPU cores to be allocated.
        pmem (numeric): default random access memory (RAM) in MB to be reserved per core.
        walltime (str): maximum time allowed for job to run on server in format days:hr:min:sec.
        rapid (str): resource allocation project identifier (RAPid) for HPC user.
        username (str): user name to associate with HPC server job.
        extraCommands (list): additional non-foldx commands to be written to job file.
        serverDataDir (Path): data directory used by HPC server relative to job directory.

    """
    produce_foldx_jobs (mutations, pdbDir, outDir / 'data', parameters = foldxParam)
    produce_guillimin_foldx_jobs (mutations,
                                  outDir / 'jobs',
                                  nodes = nodes,
                                  ppn = ppn,
                                  pmem = pmem,
                                  walltime = walltime,
                                  rapid = rapid,
                                  username = username,
                                  extraCommands = extraCommands,
                                  serverDataDir = serverDataDir,)

def produce_foldx_and_beluga_jobs (mutations,
                                   pdbDir,
                                   outDir,
                                   foldxParam = None,
                                   account = 'ctb-yxia',
                                   walltime = '1-00',
                                   ntasks = 1,
                                   nodes = 1,
                                   ntasks_per_node = 1,
                                   cpus_per_task = 1,
                                   mem = None,
                                   mem_per_cpu = '4000M',
                                   outputfile = '%x-%j.out',
                                   errorfile = None,
                                   username = '',
                                   extraCommands = None,
                                   serverDataDir = '../data'):
    """Produce data files for FoldX ∆∆G calculations as well as job files for Beluga server.
        See https://docs.computecanada.ca/wiki/Béluga/en

    Args:
        mutations (dict): mutations associated with each structural model.
        pdbDir (Path): file directory containing PDB structures.
        outDir (Path): file directory to save jobs to.
        foldxParam (dict): FoldX parameters used for each mutation job, otherwise default.
        account (str): project account name.
        walltime (str): maximum time allowed for job to run.
        ntasks (numeric): number of processes to be allocated.
        nodes (numeric): number of server nodes to be allocated.
        ntasks_per_node (numeric): number of processes to be allocated per node.
        cpus_per_task (numeric): number of nodes to be allocated per process.
        mem (str): memory per node.
        mem_per_cpu (str): memory per core.
        outputfile (Path): path to file where standard output is written.
        errorfile (Path): path to file where runtime error is written.
        username (str): user name to associate with Beluga server job.
        extraCommands (list): additional non-foldx commands to be written to job file.
        serverDataDir (Path): data directory used by Beluga server relative to job directory.

    """
    produce_foldx_jobs (mutations, pdbDir, outDir / 'data', parameters = foldxParam)
    produce_beluga_foldx_jobs (mutations,
                               outDir / 'jobs',
                               account = account,
                               walltime = walltime,
                               ntasks = ntasks,
                               nodes = nodes,
                               ntasks_per_node = ntasks_per_node,
                               cpus_per_task = cpus_per_task,
                               mem = mem,
                               mem_per_cpu = mem_per_cpu,
                               outputfile = outputfile,
                               errorfile = errorfile,
                               username = username,
                               extraCommands = extraCommands,
                               serverDataDir = serverDataDir)

def produce_foldx_jobs (mutations, pdbDir, outDir, parameters = None):
    """Produce jobs to be submitted to FoldX server for ∆∆G calculations.

    Args:
        mutations (dict): mutations associated with each structural model.
        pdbDir (Path): file directory containing PDB structures.
        outDir (Path): file directory to save FoldX jobs to.
        parameters (dict): FoldX parameters used for each mutation job, otherwise default.

    """
    clear_structures()
    if not outDir.exists():
        os.makedirs(outDir)
    
    default_param = {'temp':298,
                     'ph':7,
                     'ionStrength':0.05,
                     'water':'-IGNORE',
                     'vdwDesign':2}
    if not parameters:
        parameters = {}
    for struc, mutList in mutations.items():
        strucid = '_'.join(struc)
        pdbid, chainID1, chainID2 = struc[:3]
        for mut in mutList:
            mutID = '_'.join([strucid, mut])
            mutDir = outDir / mutID
            if not mutDir.exists():
                os.makedirs(mutDir)
            write_partial_structure (pdbid,
                                     [chainID1, chainID2],
                                     pdbDir,
                                     mutDir / (strucid + '.pdb'))
            write_foldx_config (mutDir / 'config_repairPDB.cfg',
                                'RepairPDB',
                                pdb_dir = '../data/%s' % mutID,
                                output_dir = '../data/%s' % mutID,
                                pdb_file = '%s.pdb' % strucid)
            mutParam = {k:v for k, v in default_param.items()}
            if (struc, mut) in parameters:
                for k, v in parameters[(struc, mut)].items():
                    mutParam[k] = v
            write_foldx_config (mutDir / 'config_pssm.cfg',
                                'Pssm',
                                other_cmd = ['analyseComplexChains=%s,%s' % (chainID1, chainID2),
                                             'aminoacids=%s' % mut[-1],
                                             'positions=%sa' % mut[:-1]],
                                pdb_dir = '../data/%s' % mutID,
                                output_dir = '../data/%s' % mutID,
                                pdb_file = '%s_Repair.pdb' % strucid,
                                temp = mutParam['temp'],
                                ph = mutParam['ph'],
                                ionStrength = mutParam['ionStrength'],
                                water = mutParam['water'],
                                vdwDesign = mutParam['vdwDesign'])

def write_foldx_config (outPath,
                        command,
                        other_cmd = None,
                        pdb_dir = None,
                        output_dir = None,
                        pdb_file = None,
                        mutant_file = None,
                        temp = 298,
                        ph = 7,
                        ionStrength = 0.05,
                        water = '-IGNORE',
                        vdwDesign = 2):
    """Produce FoldX configuration file.

    Args:
        outPath (Path): file path to save FoldX configurations.
        command (str): main command to be run by FoldX.
        other_cmd (list): additional commands in string format to be run by FoldX.
        pdb_dir (Path): file directory containing PDB structures.
        output_dir (Path): file directory where FoldX results are saved.
        pdb_file (Path): path to file containing PDB structure to be processed.
        mutant_file (Path): path to file containing list of mutations if applicable.
        temp (numeric): temperature in Kelvins.
        ph (numeric): PH value.
        ionStrength (numeric): ionic strength of the solution in Moles.
        water (str): how to handle water molecules, '-CRYSTAL', '-PREDICT', '-IGNORE' or '-COMPARE'.
        vdwDesign (numeric): VDW design of the experiment, 0 very soft, 1 medium soft, 2 strong.

    """
    with io.open(outPath, "w") as fout:
        fout.write('command=%s' % command)
        if other_cmd:
            fout.write('\n' + '\n'.join(other_cmd))
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

def produce_guillimin_foldx_jobs (mutations,
                                  outDir,
                                  nodes = 1,
                                  ppn = 1,
                                  pmem = 7700,
                                  walltime = '1:00:00:00',
                                  rapid = None,
                                  username = '',
                                  extraCommands = None,
                                  serverDataDir = '../data'):
    """Produce Guillimin server job files specific to FoldX ∆∆G calculations.
        See http://www.hpc.mcgill.ca

    Args:
        mutations (dict): mutations associated with each structural model.
        outDir (Path): file directory to save jobs to.
        nodes (numeric): number of server nodes to be allocated.
        ppn (numeric): total number of CPU cores to be allocated.
        pmem (numeric): default random access memory (RAM) in MB to be reserved per core.
        walltime (str): maximum time allowed for job to run on server in format days:hr:min:sec.
        rapid (str): resource allocation project identifier (RAPid) for HPC user.
        username (str): user name to associate with HPC server job.
        extraCommands (list): additional non-foldx commands to be written to job file.
        serverDataDir (Path): data directory used by HPC server relative to job directory.

    """
    if not outDir.exists():
        os.makedirs(outDir)
    
    for struc, mutList in mutations.items():
        strucid = '_'.join(struc)
        for mut in mutList:
            mutID = '_'.join([strucid, mut])
            commands = ['../foldx -f %s/%s/config_repairPDB.cfg' % (serverDataDir, mutID),
                        '../foldx -f %s/%s/config_pssm.cfg' % (serverDataDir, mutID)]
            if extraCommands:
                    commands = extraCommands + commands
            write_guillimin_job (outDir / (mutID + '_job.txt'),
                                 nodes = nodes,
                                 ppn = ppn,
                                 pmem = pmem,
                                 walltime = walltime,
                                 outputfile = '%s/%s/outputfile' % (serverDataDir, mutID),
                                 errorfile = '%s/%s/errorfile' % (serverDataDir, mutID),
                                 rapid = rapid,
                                 jobid = username + '_foldx_' + mutID,
                                 commands = commands)

def produce_beluga_foldx_jobs (mutations,
                               outDir,
                               account = 'ctb-yxia',
                               walltime = '1-00',
                               ntasks = 1,
                               nodes = 1,
                               ntasks_per_node = 1,
                               cpus_per_task = 1,
                               mem = None,
                               mem_per_cpu = '4000M',
                               outputfile = '%x-%j.out',
                               errorfile = None,
                               username = '',
                               extraCommands = None,
                               serverDataDir = '../data'):
    """Produce Beluga server job files specific to FoldX ∆∆G calculations.
        See https://docs.computecanada.ca/wiki/Béluga/en

    Args:
        mutations (dict): mutations associated with each structural model.
        outDir (Path): file directory to save jobs to.
        account (str): project account name.
        walltime (str): maximum time allowed for job to run.
        ntasks (numeric): number of processes to be allocated.
        nodes (numeric): number of server nodes to be allocated.
        ntasks_per_node (numeric): number of processes to be allocated per node.
        cpus_per_task (numeric): number of nodes to be allocated per process.
        mem (str): memory per node.
        mem_per_cpu (str): memory per core.
        outputfile (Path): path to file where standard output is written.
        errorfile (Path): path to file where runtime error is written.
        username (str): user name to associate with Beluga server job.
        extraCommands (list): additional non-foldx commands to be written to job file.
        serverDataDir (Path): data directory used by Beluga server relative to job directory.

    """
    if not outDir.exists():
        os.makedirs(outDir)
    
    for struc, mutList in mutations.items():
        strucid = '_'.join(struc)
        for mut in mutList:
            mutID = '_'.join([strucid, mut])
            commands = ['../foldx -f %s/%s/config_repairPDB.cfg' % (serverDataDir, mutID),
                        '../foldx -f %s/%s/config_pssm.cfg' % (serverDataDir, mutID)]
            if extraCommands:
                commands = extraCommands + commands
            write_beluga_job (outDir / (mutID + '_job.sh'),
                              account = account,
                              walltime = walltime,
                              ntasks = ntasks,
                              nodes = nodes,
                              ntasks_per_node = ntasks_per_node,
                              cpus_per_task = cpus_per_task,
                              mem = mem,
                              mem_per_cpu = mem_per_cpu,
                              outputfile = outputfile,
                              errorfile = errorfile,
                              jobname = username + '_foldx_' + mutID,
                              commands = commands)

def read_foldx_results (inDir):
    """Read mutation ∆∆G results produced by FoldX.

    Args:
        inDir (Path): file directory containing FoldX results.

    Returns:
        dict, dict: processed and unprocessed mutations.

    """
    processed, unprocessed = {}, {}
    strucDirs = os.listdir(inDir)
    strucDirs = [dir for dir in strucDirs if os.path.isdir(inDir / dir)]
    for strucID in strucDirs:
        strucDir = inDir / strucID
        struc = tuple(strucID.split('_'))
        if len(struc) > 3:
            struc = struc[:-1]
        
        mutListFile = resultFile = None
        for filename in os.listdir(strucDir):
            if filename.startswith('individual_list'):
                mutListFile = strucDir / filename
                with io.open(mutListFile, "r") as f:
                    mutList = list(map(str.strip, f.read().split(';')))
                    mutList.remove('')
            elif filename.startswith('_'.join(('PSSM',) + struc)):
                resultFile = strucDir / filename
        
        if resultFile:
            with io.open(resultFile, "r") as f:
                mutResidues = list(map(str.strip, f.readline().strip().split('\t')))
                for line in f:
                    ddgList = list(map(str.strip, line.split('\t')))
                    if len(ddgList) > 1:
                        wt = ddgList[0]
                        for mt, ddg in zip(mutResidues, ddgList[1:]):
                            processed[struc + (wt + mt,)] = float(ddg) if len(ddg) else 'X'
                            
            for mut in mutList:
                if struc + (mut,) not in processed:
                    processed[struc + (mut,)] = 'X'
        else:
            if len(mutList) == 1:
                processed[struc + (mutList.pop(),)] = 'X'
            elif len(mutList) > 1:
                for mut in mutList:
                    unprocessed[struc + (mut,)] = [mut]
    
    return processed, unprocessed

def read_protein_mutation_ddg (inPath, type = 'binding'):
    """Read protein mutation ∆∆G values from file.

    Args:
        inPath (Path): path to file containing mutations with ∆∆G values.
        type (str): type of ∆∆G values, 'binding' for interface ∆∆G, 'folding' for protein folding ∆∆G.

    Returns:
        dict

    """
    ddgDict = {}
    with io.open(inPath, "r", encoding="utf-8") as f:
        next(f)
        for line in f:
            linesplit = list( map ( str.strip, line.split('\t') ) )
            if len(linesplit) > 8:
                if linesplit[8] is not 'X':
                    protein, partner, pr_pos, pdbid, chainID, ch_pos, ch_mut, ch_partner, ddg = linesplit
                    if type is 'binding':
                        k = protein, partner, int(pr_pos), ch_mut[-1]
                        val = pdbid, chainID, ch_partner, ch_mut, float(ddg)
                    elif type is 'folding':
                        k = protein, int(pr_pos), ch_mut[-1]
                        val = pdbid, chainID, ch_mut, float(ddg) 
                    if k not in ddgDict:
                        ddgDict[k] = val
    return ddgDict

def read_chain_mutation_ddg (inPath, type = 'binding'):
    """Read PDB chain mutation ∆∆G values from file.

    Args:
        inPath (Path): path to file containing mutations with ∆∆G values.
        type (str): type of ∆∆G values, 'binding' for interface ∆∆G, 'folding' for protein folding ∆∆G.

    Returns:
        dict

    """
    ddgDict = {}
    with io.open(inPath, "r", encoding="utf-8") as f:
        next(f)
        for line in f:
            linesplit = list( map ( str.strip, line.split('\t') ) )
            if len(linesplit) > 8:
                protein, partner, pr_pos, pdbid, chainID, ch_pos, ch_mut, ch_partner, ddg = linesplit[:9]
                if type is 'binding':
                    k = (pdbid,) + tuple(sorted([chainID, ch_partner])) + (ch_mut,)
                elif type is 'folding':
                    k = pdbid, chainID, ch_mut
                ddgDict[k] = ddg
    return ddgDict

def copy_mutation_ddg (inPath1, inPath2, outPath, type = 'binding'):
    """Transfer mutation ∆∆G values from one file to another file based on similar mutation 
        structure mappings.

    Args:
        inPath1 (Path): path to file containing mutation ∆∆G values.
        inPath2 (Path): path to file containing mutations with missing ∆∆G values.
        outPath (Path): file path to save a copy of file in inPath2 with updated ∆∆G values.

    """
    ddg = read_chain_mutation_ddg (inPath1, type)
    write_mutation_ddg_tofile (ddg, inPath2, outPath, type)

def write_mutation_ddg_tofile (ddg, inPath, outPath, type = 'binding'):
    """Update file with mutation ∆∆G values.

    Args:
        ddg (dict): mutation ∆∆G values.
        inPath (Path): path to file whose mutations will be updated with ∆∆G values.
        outPath (Path): file path to save a copy of file in inPath with updated ∆∆G values.
        type (str): type of ∆∆G values, 'binding' for interface ∆∆G, 'folding' for protein folding ∆∆G.
    """
    
    with io.open(inPath, "r", encoding="utf-8") as f, io.open(outPath, "w") as fout:
        fout.write(f.readline().strip() + '\n')
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
            fout.write('\t'.join(map(str, strsplit)) + '\n')

def append_mutation_ddg_files (inPath1, inPath2, outPath):
    """Append two ∆∆G files together and save to anothor file.

    Args:
        inPath1 (Path): path to first input file.
        inPath2 (Path): path to second input file.
        outPath (Path): file path to save appended output to.
    """
    
    with io.open(outPath, "w") as fout:
        with io.open(inPath1, "r", encoding="utf-8") as f1:
            fout.write(f1.readline().strip() + '\n')
            for line in f1:
                fout.write(line)
        with io.open(inPath2, "r", encoding="utf-8") as f2:
            next(f2)
            for line in f2:
                fout.write(line)
