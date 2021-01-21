"""
In Sequence pattern matching there are three methods of achieving this:
1. Naive Algorithm or method which scans from left to right finding all possible occurrences of the pattern
   at every sub-sequence of length equals to that of the pattern, at successive positions till the pattern is found
   or if all other occurrences are needed it scans through the sequence matching the pattern

2. Heuristic Algorithm one of which is the Boyer-moore method: A major disadvantage of the naive algorithm is 
   computational efficiency when the sequence is large. This algorithm tries to improve on the average 
   computational efficiency of pattern searching by trying to use the structure of the pattern to speed up 
   the search. 
   
   In the Boyer-Moore algorithm, the worst case scenario is the runtime similar to that of the naive algorithm
   The algorithm is based on two rules that allow to move forward more than one position in the target sequence 
   in some situations. Here the pattern searching is done from right to left. When there is a mismatch between 
   the target sequence and the pattern, two rules can be applied to check if the process can be more efficient, 
   by moving forward more than one position in the sequence
   
   One rule is the bad-character rule (bcr): It advances the pattern to the next occurrence of the symbol in the 
   sequence at the position of the mismatch. If no occurrences of that symbol exist in the pattern, we can move 
   forward the maximum number of symbols until the end of the pattern.
   The other rule is the good-suffix rule (gsr): It states that when there's a mismatch, we can move forward to 
   the next instance in the pattern of the part (suffix) that matched before (in the right) the mismatch.

"""


def naiveAlgo(seq, pat):
    pat_list = []
    for i in range(len(seq)):
        match = seq[i:i+len(pat)]
        if pat == match:
            pat_list.append(i)
    return pat_list


def Boyer_Moore(seq, alphabet, pat):
    res = []
    occ, s = preprocess(alphabet, pat)
    i = 0
    while i <= len(seq) - len(pat):
        j = len(pat)- 1
        while j >= 0 and pat[j] == seq[j+i]:
            j -= 1
        if j < 0:
            res.append(i)
            i += s[0]
        else:
            c = seq[j+1]
            i += max(s[j+1], j-occ[c])
    return res

def preprocess(alphabet, pat):
    """
    Implementation of the Boyer-Moore algorithm
    """
    occ = preprocess_bcr(alphabet, pat)
    s = preprocess_gsr(pat)
    return occ, s


def preprocess_bcr(alphabet, pat):
    occ = {}
    for symb in alphabet:
        occ[symb] = -1
    for j in range(len(pat)):
        c = pat[j]
        occ[c] = j
    return occ


def preprocess_gsr(pat):
    f = [0] * (len(pat)+1)
    s = [0] * (len(pat)+1)
    i = len(pat)
    j = len(pat)+1
    f[i] = j

    while i > 0:
        while j <= len(pat) and pat[i-1] != pat[j-1]:
            if s[j] == 0:
                s[j] = j-i
            j =  f[j]
        i -= 1
        j -= 1
        f[i] = j
    j = f[0]
    for i in range(len(pat)):
        if s[i] == 0:
            s[i] = j
        if i == j:
            j = f[j]
    return s


def RegEx_pattern(seq, pat):
    """
    This gets all pattern occurrences of the a given base pair pattern
    and returns the frequency
    """
    import re
    pat_list = []
    pattern = re.compile(pat)
    all_occurrences = re.finditer(pattern, seq)
    for match in all_occurrences:
        pat_list.append(match.span()[0])
    print(pat_list)