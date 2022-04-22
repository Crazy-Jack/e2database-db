mysql \
  --user=genome \
  -N \
  --host=genome-mysql.cse.ucsc.edu \
  -A \
  -D hg19 \
  -e "select ensGene.name, name2, chrom, strand, txStart, txEnd, value from ensGene, ensemblToGeneName where ensGene.name = ensemblToGeneName.name" > \
  output.txt