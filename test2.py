from Bio.Seq import Seq
from binformatics import sequences


gene_seq = sequences.generate_genes('ATCG', sample=150)
my_seq= Seq(gene_seq)
print(my_seq)