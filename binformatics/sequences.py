
def generate_genes(genes, num=100, sample=250):
    """
    Generates gene sequences of length
    """
    import  random
    gene_seq = genes.strip() * num
    gene_seq = random.sample(gene_seq, sample,) # samples 250 sequences from the list 
    return ''.join(gene_seq).upper()


def validate_seq(gene_seq):
    """
    Checks and validates if the gene sequence is a gene composed of AGCT
    Returns True if gene seqeunce is valid or otherwise
    """
    base_counts = gene_seq.count('A') + gene_seq.count('G') + \
                gene_seq.count('T') + gene_seq.count('C')
    
    valid = len(gene_seq)
    return valid == base_counts 


def frequency(gene_seq):
    """This counts the number of AGTC in the gene_seq
    """
    from collections import Counter
    
    assert validate_seq(gene_seq), 'Invalid gene sequence'
    gene_counts = Counter(gene_seq)

    gene_counts = {i:j for i,j in sorted(gene_counts.items(), key=lambda x: x[0])}
    return gene_counts


def complementaryDNA(gene_seq):
    """
    Returns the complementary DNA of a gene sequence
    """
    assert validate_seq(gene_seq), 'Invalid gene sequence'
    dna_map = {'A':'T', 'G':'C', 'C':'G', 'T':'A'}
    cDNA = ''.join(dna_map[base] for base in gene_seq)
    return cDNA


def replication(gene_seq):
    """This function performs replication of a DNA Strand. 
    In replication, the double helix strands of the DNA are unwound and in each strand, 
    a complementary strand (daughter strand) is replicated, forming a parent-daughter strand
    annealing. 
    Although the other strand is replicated, we are assuming replication of one strand here.
    """
    assert validate_seq(gene_seq), 'Invalid DNA sequence'
    first_strand = gene_seq
    second_strand = complementaryDNA(gene_seq)
    return first_strand, second_strand[::-1] #reversing into 5! ->3!


def reverse_transcription(rna_seq):
    count = rna_seq.count('A') + rna_seq.count('G') + rna_seq.count('U') + rna_seq.count('C')
    assert count == len(rna_seq), 'Invalid RNA sequence'
    dna_seq = rna_seq.replace('U', 'T')
    return dna_seq


def nucleotide_base_count(seq):
    """
    This counts the number of purine and pyrimidine bases in the RNA or DNA sequence
    """
    if 'T' in seq:
        count = seq.count('A') + seq.count('G') + seq.count('T') + seq.count('C')
        if len(seq) != count:
            print('Invalid sequence')
        else:
            purines = seq.count('A') + seq.count('G')
            pyrimidines = seq.count('T') + seq.count('C')
            print(f'Purines: {purines}, Pyrimidines: {pyrimidines}')
    elif 'U' in seq:
        count = seq.count('A') + seq.count('G') + seq.count('U') + seq.count('C')
        if len(seq) != count:
            print('Invalid sequence')
        else:
            purines = seq.count('A') + seq.count('G')
            pyrimidines = seq.count('U') + seq.count('C')
            print(f'Purines: {purines}, Pyrimidines: {pyrimidines}')


def transcription(gene_seq):
    assert validate_seq(gene_seq), 'Invalid gene sequence'
    return gene_seq.replace('T', 'U')


def gc_content(gene_seq):
    """Returns the proportion of Guanine-Cytosine content in the gene sequence
    """
    assert validate_seq(gene_seq), 'Invalid gene sequence'
    freq = frequency(gene_seq)

    gc_content = (freq['G'] + freq['C']) / sum(freq.values())
    return gc_content


def three_letter_names():
    """Returns the three letter word for the one letter encoded protein
    """
    maps = {
        'M' : 'MET', 'P' : 'PHE', 'L' : 'LEU', 'G':'GLY',
        'A' : 'ALA', 'V' : 'VAL', 'I' : 'ILE', 'R':'ARG',
        'C' : 'CYS', 'K' : 'LYS', 'G' : 'GLU', 'D':'ASP',
        'Q' : 'GLN', 'H' : 'HIS', 'T' : 'THR', 'P':'PRO',
        'S' : 'SER', 'Y' : 'TYR', 'N' : 'ASN', 'W':'TRP'
        }
    return maps


def orfs(rna_seq, n_rf, protein_map):
    """In Translation, there are three possible ways a ribosome parses through an mRNA sequence
    These parses are called reading frames.
    On parsing the reading frames, when an AUG codon is met, translation starts, till a stop codon 
    the ribomsome sees a stop codon, leading to its termination.
    A reading frame with sufficient length and without a stop codon is called an open reading frame (ORF)

    Args
    ====
    rna_seq = transcribed DNA sequence (ie the mRNA sequence)
    n_rf = the three possible reading frames for one strand of transcribed DNA

    Returns
    =======
    The proteins formed from each reading frame
    """
    protein_3_map = three_letter_names()
    proteins = ""
    protein3 = ""
    
    start = False  # initiation not started
    for i in range(n_rf, len(rna_seq)-2, 3):
        codon = rna_seq[i:i+3]
        if codon == 'AUG':
            start = True # switching initiation on
        if start and len(codon) == 3:
            protein = protein_map[codon]
            if protein == '_': # stop translation if a stop codon is reached
                break
            proteins += '-' + protein
    start = False # turning initiation off as soon as translation is terminated or end of gene reached

    # converting to three letter protein name
    for one_letter in proteins.split('-'):
        if one_letter in protein_3_map:
            protein3 += '-' + protein_3_map[one_letter].capitalize()
    return protein3


def translation(gene_seq):
    """
    Translates the RNA to proteins
    Args:
    =====
    gene_seq = DNA sequence
    """
    from ast import literal_eval

    # checking if DNA sequence is valid
    assert validate_seq(gene_seq), 'Invalid DNA sequence'
    
    # DNA double strands
    first_strand, second_strand = replication(gene_seq)
    
    # transcribing to RNA sequences
    rna_seq1 = transcription(first_strand)
    rna_seq2 = transcription(second_strand)

    # importing protein map
    with open('binformatics/proteins.txt', 'r+') as file:
        protein_map = file.readlines()
        protein_map = ''.join(protein_map).strip()
        protein_map = literal_eval(protein_map)
        possible_proteins = [] # gets all possible proteins at all ORFs (open reading frames)
        final_prots = [] # get all possible proteins that aren't empty

        for n_rf in range(3):
            proteins1 = orfs(rna_seq1, n_rf, protein_map)
            proteins2 = orfs(rna_seq2, n_rf, protein_map)
            possible_proteins.append(proteins1)
            possible_proteins.append(proteins2)
        
        # removing empty proteins
        for prot in possible_proteins:
            if len(prot) > 0:
                final_prots.append(prot)
        # final_prots.sort()
        return final_prots

