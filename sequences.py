
def generate_genes(genes, num=100):
    """
    Generates gene sequences of length
    """
    from random import sample
    gene_seq = genes.strip() * num
    gene_seq = sample(gene_seq, 250) # samples 500 sequences from the list 
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


def transcription(gene_seq):
    assert validate_seq(gene_seq), 'Invalid gene sequence'
    return gene_seq.replace('T', 'U')


def gc_content(gen_seq):
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
    
    # transcribing to RNA sequences
    rna_seq = transcription(gene_seq)

    def orfs(rna_seq, n_rf):
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
        for i in range(n_rf, len(rna_seq)-3, 3):
            codon = rna_seq[i:i+3]
            if codon == 'AUG':
                start = True
            if start:
                protein = protein_map[codon]
                if protein == '_': # stop translation if a stop codon is reached
                    break
                proteins += '-' + protein
        start = False # turning initiation off as soon as translation is terminated

        # converting to three letter protein name
        for one_letter in proteins.split('-'):
            if one_letter in protein_3_map:
                protein3 += '-' + protein_3_map[one_letter]
        return protein3

    # importing protein map
    with open('proteins.txt', 'r+') as file:
        protein_map = file.readlines()
        protein_map = ''.join(protein_map).strip()
        protein_map = literal_eval(protein_map)
        possible_proteins = [] # gets all possible proteins at all ORFs (open reading frames)

        # ORF 1
        for n_rf in range(3):
            proteins = orfs(rna_seq, n_rf)
            possible_proteins.append(proteins)

        return possible_proteins







if __name__ == '__main__':
    genes = 'AGTC'
    gene_seq = generate_genes(genes)
    
    freq = frequency(gene_seq)
    print(freq)
    
    cdna = complementaryDNA(gene_seq)
    rna = transcription(gene_seq)
    proteins = translation(gene_seq)
    # print(proteins)
    # print(rna[:50])
    # print(gene_seq[:50])
    # print(cdna[:50])
    gc_prop = gc_content(gene_seq)
    print(gc_prop)