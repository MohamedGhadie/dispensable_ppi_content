#----------------------------------------------------------------------------------------
# Modules for simple computations.
#----------------------------------------------------------------------------------------

import os
import re
from pathlib import Path
from random import choice

def valid_uniprot_id(s):
    """Check if a string is a valid UniProt accession ID.

    See: http://www.uniprot.org/help/accession_numbers

    Args:
        s (str): string to check.

    Returns:
        bool

    """
    e = '[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2}'
    return bool(re.match(e, s))

def create_dir (dir):
    """Create a file directory.

    Args:
        dir (Path): directory to create.

    """
    if len(dir) > 0:
        dir = Path(dir)
        if not dir.exists():
            os.makedirs(dir, exist_ok=True)

def is_numeric (s):
    """Check if string contains numeric characters only.

    Args:
        s (str): string to check.

    Returns:
        bool

    """
    try:
        float(s)
        return True
    except ValueError:
        return False

def hamming_dist (a, b):
    """Calculate hamming distance between two vectors.

    Args:
        a (list): first vector.
        b (list): second vector.

    Returns:
        int

    """
    assert len(a) == len(b), 'Vectors of unequal lengths passed to hamming distance'
    return sum(c1 != c2 for c1, c2 in zip(a, b))

def position_overlap (s1, e1, s2, e2):
    """Determine if two sequence positions overlap.
	
    Args:
        s1 (int): starting position for first sequence.
        e1 (int): ending position for first sequence.
        s2 (int): starting position for second sequence.
        e2 (int): ending position for second sequence.
    
    Returns:
        bool
    
    """ 
    return (e1 >= s2) and (e2 >= s1)

def segment_overlaps (seq, otherSeqs):
    """Calculate fraction of continuous positions in a sequence segment that overlap with 
        a list of other sequence segments.
	
    Args:
        seq1 (tuple): start and end positions in first sequence segment.
        otherSeqs (list): list of start and end position tuples for other sequence segments.
    
    Returns:
        list
    """
    return [segment_overlap(seq, s) for s in otherSeqs]

def segment_overlap (seq1, seq2):
    """Calculate fraction of continuous positions in a sequence segment that overlap with 
        another sequence segment.
	
    Args:
        seq1 (tuple): start and end positions in first sequence.
        seq2 (tuple): start and end positions in second sequence.
    
    Returns:
        numeric
    
    """
    start1, end1 = seq1
    start2, end2 = seq2
    if (end1 >= start2) and (end2 >= start1):
        return (min(end1, end2) - max(start1, start2) + 1) / (end1 - start1 + 1)
    else:
        return 0

def find_char (ch, s):
    """Find all occurence indices of a character in a string.

    Args:
        ch (str): character to find.
        s (str): string to search.

    Returns:
        list

    """
    return [i for i, ltr in enumerate(s) if ltr == ch]

def first_mismatch (s1, s2):
    """Find index of first mismatch between two strings.

    Args:
        s1 (str): first string to be compared.
        s2 (str): second string to be compared.

    Returns:
        int: index of first mismatch, -1 if no mismatch found.

    """
    maxIndex = min(len(s1),len(s2))
    mismatch = [i for i in range(maxIndex) if s1[i] != s2[i]]
    return mismatch[0] if len(mismatch) > 0 else -1

def find_substring (s1, s2):
    """Find all occurence indices of one string in another string.

    Args:
        s1 (str): string to find.
        s2 (str): string to search.

    Returns:
        list

    """
    s1 = s1.replace('*','\\*')
    return [m.start() for m in re.finditer('(?=%s)' % s1, s2)]

def find_masked_substring (s1, s2, maskPos):
    """Find all occurence indices of one string in another string, ignoring a single 
        position in the first string.

    Args:
        s1 (str): string to find.
        s2 (str): string to search.
        maskPos (int): index of position to ignore.

    Returns:
        list

    """
    if 0 <= maskPos < len(s1):
        s1 = s1[:maskPos] + '.' + s1[maskPos + 1:]
    s1 = s1.replace('*', '\\*')
    return [m.start() for m in re.finditer('(?=%s)' % s1, s2)]

def reverseTuples(s):
    """Reverse the order of each pair in a list of pairs.

    Args:
        s (list or str):  list of pairs to reverse. If comma-separated string, each pair
                            must be separated by '-'.

    Returns:
        list

    """
    if isinstance(s, str):
        s = map(str.strip, s.split(','))
    return list(map(reverseTuple, s))

def reverseTuple(s):
    """Reverse the order of a pair.

    Args:
        s (list, str):  pair to reverse. If string, pair must be separated by '-'.

    Returns:
        tuple

    """
    if isinstance(s, str):
        s = tuple(map(str.strip, s.split('-')))
    return tuple(reversed(s))

def str_to_tuples(s):
    """Convert a string or list of comma-separated pairs into a list of tuples

    Args:
        s (str, list): pairs to convert to tuples.

    Returns:
        list

    """
    if isinstance(s, str):
        if s == '-':
            return []
        else:
            s = map(str.strip, s.split(','))
    return list(map(str_to_tuple, s))

def str_to_tuple(s):
    """Convert a comma-separated pair into a tuple.

    Args:
        s (str): string to convert to tuple.

    Returns:
        tuple

    """
    if isinstance(s, str):
        s = tuple(map(str.strip, s.split('-')))
    return s

def listStr_to_list(s):
    """Split a comma-separated string into a list.

    Args:
        s (str): comma-separated string to split.

    Returns:
        list

    """
    if isinstance(s, str):
        s = list(map(str.strip, s.split(',')))
    return s

def list_to_str (ls, delm):
    """Join a multi-layer list of lists into a string using a series of delimiters.

    Args:
        ls (list): list with nested list elements.
        delm (list): list of delimiters to split multiple layers of the list.

    Returns:
        str

    """
    if len(ls) > 0:
        if len(delm) > 1:
            return delm[0].join([list_to_str(elm, delm[1:]) for elm in ls])
        else:
            return delm[0].join(map(str, ls))
    else:
        return '-'

def str_to_list (s, delm, dtype):
    """Split a string into a multi-layer list using a series of specified delimiters.

    Args:
        s (str): string to split.
        delm (list): list of delimiters to split multiple layers of the list.
        dtype (type): data type of output list elements.

    Returns:
        list

    """
    if s != '-':
        if len(delm) > 1:
            return [str_to_list(elm, delm[1:], dtype) for elm in s.split(delm[0])]
        else:
            return list(map(dtype, s.split(delm[0])))
    else:
        return []

def isolate_pairs (ls):
    """Extract pairs of elements from a list of tuples where each tuple contains two lists
        of the elements to be combined in pairs.

    Args:
        ls (list): list of tuples whose two lists of elements are to be paired up.

    Returns:
        list

    """
    isolatedPairs = []
    for i in ls:
        isolatedPairs.extend(extract_pairs(i))
    return isolatedPairs

def extract_pairs (ls):
    """Extract pairs of elements from a list of two lists.

    Args:
        ls (list): list whose nested list elements are to be paired up.

    Returns:
        list

    """
    lefts, rights = ls
    if (len(lefts) > 0) and (len(rights) > 0):
        return [(l,r) for l in lefts for r in rights]
    elif (len(lefts) > 0) and (len(rights) == 0):
        return [(l,[]) for l in lefts]
    elif (len(lefts) == 0) and (len(rights) > 0):
        return [([],r) for r in rights]
    else:
        return ([], [])

def merge_list_pairs (ls):
    """Merge a list of tuples of lists into one tuple. Returned tuple elements are lists 
        with no duplicates.

    Args:
        ls (list): list of tuples to be merged.

    Returns:
        tuple

    """
    elm1, elm2 = [], []
    for e1, e2 in ls:
        elm1.extend(e1)
        elm2.extend(e2)
    return sorted(set(elm1)), sorted(set(elm2))

def sample_random_pairs (ls, N):
    """Generate a list of random pairs sampled with replacement from a list of elements.
        Elements of each pair are non-identical.

    Args:
        ls (list): elements to sample pairs from.
        N (int): number of pairs to sample.

    Returns:
        list

    """
    pairs, n = [], 0
    while n < N:
        p1, p2 = choice(ls), choice(ls)
        if p1 != p2:
            pairs.append((p1, p2))
            n += 1
    return pairs
