import io
import csv
import pandas as pd
from simple_tools import listStr_to_list

def parse_domine_DDI_file(inPath, outPath):
    """Read domain-domain interactions (DDIs) from DOMINE file.

    Args:
        inPath (str): DOMINE DDI file directory.
        outPath (str): file directory to save DDIs to.

    """
    interactions = pd.read_table(inPath, sep="|")
    DDIs = interactions.loc[(interactions["iPfam"]==1)|(interactions["3did"]==1), ["P1","P2"]]
    DDIs.rename(columns={"P1":"dom1","P2":"dom2"}, inplace=True)
    print('%d domain-domain interactions extracted from file ' % len(DDIs) +
          str(inPath))
    DDIs.to_csv(outPath, index=False, sep='\t')

def parse_DDI_files(_3didFile, domineFile, _3didOutPath, domineOutPath, outPath):
    """Read domain-domain interactions (DDIs) from DOMINE and 3did files.

    Args:
        _3didFile (str): 3did DDI file directory.
        domineFile (str): DOMINE DDI file directory.
        _3didOutPath (str): file directory to save 3did DDIs to.
        domineOutPath (str): file directory to save DOMINE DDIs to.
        outPath (str): file directory to save all DDIs to, with no repetition.

    """
    parse_3did_DDI_file(_3didFile, _3didOutPath)
    parse_domine_DDI_file(domineFile, domineOutPath)
    _3didDDIs = pd.read_table(_3didOutPath, sep='\t')
    _3didDDIs = _3didDDIs[["dom1","dom2"]]
    domineDDIs = pd.read_table(domineOutPath, sep='\t')
    DDIs = _3didDDIs.append(domineDDIs, ignore_index=True)
    DDIs = DDIs.drop_duplicates(subset=["dom1","dom2"]).reset_index(drop=True)
    # drop duplicates that are in reverse order
    DDIs["duplicate"] = False
    for i, row in DDIs.iterrows():
        isduplicate = list((DDIs["dom1"]==row.dom2) & (DDIs["dom2"]==row.dom1))
        if sum(isduplicate) == 1:
            DDIs.loc[isduplicate,"duplicate"] = (isduplicate.index(True)>i)
    DDIs = DDIs.loc[DDIs['duplicate']==False,["dom1","dom2"]]
    DDIs.to_csv(outPath, index=False, sep='\t')
    print('%d domain-domain interactions written to file ' % len(DDIs) +
          str(outPath))
    
def parse_3did_DDI_file(inPath, outPath):
    """Read domain-domain interactions (DDIs) from 3did file.

    Args:
        inPath (str): 3did DDI file directory.
        outPath (str): file directory to save DDIs to.

    """
    numLines = 10000
    c = -1
    DDIs = pd.DataFrame('NA', index=range(numLines),
                        columns=["dom1", "dom2", "dom1_name", "dom2_name"])
    with io.open(inPath, "r", encoding="utf-8") as f:
        reader = csv.reader(f)
        for row in reader:
            if len(row[0]) > 4:
                if row[0][:4] == '#=ID':
                    ddi = row[0].split('\t')
                    if len(ddi) == 5:
                        if ((ddi[3].find('@Pfam') > -1) and
                                (ddi[4].find('@Pfam') > -1)):
                            s1 = ddi[3].find('PF')
                            e1 = ddi[3].find('.')
                            s2 = ddi[4].find('PF')
                            e2 = ddi[4].find('.')
                            if (-1 < s1 < e1) and (-1 < s2 < e2):
                                c += 1
                                DDIs.loc[c, ["dom1", "dom2", "dom1_name", "dom2_name"]] = ddi[3][s1:e1], ddi[4][s2:e2], ddi[1].strip(), ddi[2].strip()
    print('%d domain-domain interactions extracted from file ' % (c + 1) +
          str(inPath))
    DDIs.to_csv(outPath, index=False, sep='\t')

def parse_3did_DDIinterface_file(inPath, outPath):
    """Read domain-domain interactions (DDIs) with interface residue interactions

    Args:
        inPath (str): 3did interface DDI file directory.
        outPath (str): file directory to save DDIs to.

    """
    numDDIs = 0
    with io.open(inPath, 'r', encoding='utf-8') as f:
        reader = csv.reader(f)
        for row in reader:
            if len(row[0]) > 4:
                if row[0][:4] == '#=ID':
                    numDDIs += 1
    c = -1
    DDIs = pd.DataFrame('NA', index=range(numDDIs),
                        columns=["dom1_name", "dom2_name","interface1_pos","interface2_pos"])
    with io.open(inPath, "r", encoding="utf-8") as f:
        reader = csv.reader(f)
        for row in reader:
            if len(row[0]) > 4:
                if row[0][:4] == '#=ID':
                    ddi = list(map(str.strip, row[0].split('\t')))
                    if len(ddi) == 3:
                        c += 1
                        print(c)
                        DDIs.loc[c, ["dom1_name", "dom2_name"]] = ddi[1], ddi[2]
                        interface1 = []
                        interface2 = []
                elif row[0][:4] == '#=IF':
                    resIx = list(map(str.strip, row[0].split('\t')))
                    if len(resIx) == 2:
                        resIx = list(map(str.strip, resIx[1].split(' ')))
                        resIx = [x.split('-') for x in resIx[1:]]
                        interface1.extend(map(int, [x for [x,_] in resIx]))
                        interface2.extend(map(int, [x for [_,x] in resIx]))
            elif row[0][:2] == '//':
                DDIs.loc[c, ["interface1_pos", "interface2_pos"]] = [sorted(list(set(interface1))),
                                                                           sorted(list(set(interface2)))]
                
    print('%d domain-domain interactions extracted from file ' % (c + 1) +
          str(inPath))
    DDIs.to_csv(outPath, index=False, sep='\t')

def read_DDIinterface_file(inPath):
    """Read DDIs with their interaction interfaces.

    Args:
        inPath (str): file directory containing DDIs with interface residue positions.
    
    Returns:
        DataFrame: DDIs with interface residue positions in list form.

    """
    DDIs = pd.read_table(inPath, sep='\t')
    DDIs["interface1_pos"] = DDIs["interface1_pos"].apply(lambda x: x[1:-1])
    DDIs["interface1_pos"] = DDIs["interface1_pos"].apply(listStr_to_list)
    DDIs["interface1_pos"] = DDIs["interface1_pos"].apply(lambda x: list(map(int,x)))
    DDIs["interface2_pos"] = DDIs["interface2_pos"].apply(lambda x: x[1:-1])
    DDIs["interface2_pos"] = DDIs["interface2_pos"].apply(listStr_to_list)
    DDIs["interface2_pos"] = DDIs["interface2_pos"].apply(lambda x: list(map(int,x)))
    return DDIs

def merge_duplicate_DDI_interfaces(inPath, outPath):
    """Merge interfaces of duplicate DDIs.

    Args:
        inPath (str): file directory containing DDIs with interface residue positions.
        outPath (str): file directory to save unique DDIs with interface residue positions.

    """
    DDIs = read_DDIinterface_file(inPath)
    mergedDDIs = pd.DataFrame('NA', index=list(range(len(DDIs))),
                                columns=("dom1",
                                         "dom2",
                                         "dom1_name",
                                         "dom2_name",
                                         "interface1_pos",
                                         "interface2_pos"))
    c = -1
    checked = pd.Series([False]*len(DDIs))
    for i, ddi in DDIs.iterrows():
        if not checked[i]:
            duplicates = (DDIs["dom1"]==ddi.dom1) & (DDIs["dom2"]==ddi.dom2)
            interface1_array = DDIs.loc[duplicates,"interface1_pos"].values
            interface1 = [pos for interface in interface1_array for pos in interface]
            interface2_array = DDIs.loc[duplicates,"interface2_pos"].values
            interface2 = [pos for interface in interface2_array for pos in interface]
            checked[duplicates] = True
            
            duplicates = (DDIs["dom2"]==ddi.dom1) & (DDIs["dom1"]==ddi.dom2)
            if duplicates.sum() > 0:
                interface1_array = DDIs.loc[duplicates,"interface2_pos"].values
                interface1b = [pos for interface in interface1_array for pos in interface]
                interface2_array = DDIs.loc[duplicates,"interface1_pos"].values
                interface2b = [pos for interface in interface2_array for pos in interface]
                checked[duplicates] = True
                interface1.extend(interface1b)
                interface2.extend(interface2b)
            c += 1
            mergedDDIs.loc[c, ["dom1",
                               "dom2",
                               "dom1_name",
                               "dom2_name"]] = (ddi.dom1,
                                                ddi.dom2,
                                                ddi.dom1_name,
                                                ddi.dom2_name)
            mergedDDIs.loc[c, ["interface1_pos",
                               "interface2_pos"]] = (sorted(list(set(interface1))),
                                                     sorted(list(set(interface2))))
    mergedDDIs = mergedDDIs[:c + 1]
    mergedDDIs.to_csv(outPath, index=False, sep='\t')
