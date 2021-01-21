
from binformatics import sequences


if __name__ == '__main__':
    genes = 'AGTC'
    gene_seq = sequences.generate_genes(genes, num=250, sample=700) 
    
    freq = sequences.frequency(gene_seq)
    print(freq)

    cdna = sequences.complementaryDNA(gene_seq)
    rna = sequences.transcription(gene_seq)
    dna = sequences.reverse_transcription(rna)
    gc_prop = sequences.gc_content(gene_seq)
    proteins = sequences.translation(gene_seq)
    sequences.nucleotide_base_count(gene_seq)
    sequences.nucleotide_base_count(rna)


    # prints
    print(proteins)
    print(rna[:50])
    print(dna[:50])
    print(gene_seq[:50])
    print(cdna[:50])
    print(gc_prop)