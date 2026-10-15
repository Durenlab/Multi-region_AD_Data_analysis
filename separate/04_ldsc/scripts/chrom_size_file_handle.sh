wget http://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/hg19.chrom.sizes

awk '$1 ~ /^chr([1-9]|1[0-9]|2[0-2])$/' \
  hg19.chrom.sizes > hg19.autosomes.chrom.sizes
