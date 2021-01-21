from binformatics import patterns
from binformatics import sequences



if __name__ == '__main__':
    genes = 'CTGA'
    gene_seq = sequences.generate_genes(genes, num=250, sample=800)
    rna = sequences.transcription(gene_seq)


    pat='GCUA'
    # b = 'ATAGAACCAATGAACCATGATGAACCATGGATACCCAACCACC'
    res = patterns.naiveAlgo(rna, pat)
    print(res)
    res = patterns.Boyer_Moore(seq=rna, alphabet='ACUG', pat=pat)
    print(res)
    res = patterns.RegEx_pattern(seq=rna, pat=pat)