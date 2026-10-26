#!/usr/bin/env bash
set -euo pipefail
: "${REGA_SEPARATE_ROOT:?Set REGA_SEPARATE_ROOT}"
: "${MAGMA_OUTPUT_ROOT:?Set MAGMA_OUTPUT_ROOT to a new output directory}"
: "${AD_GWAS_GRCH38_FILE:?Set AD_GWAS_GRCH38_FILE}"
: "${MAGMA_GENE_LOC:?Set MAGMA_GENE_LOC}"
: "${LDSC_CHAIN_FILE:?Set LDSC_CHAIN_FILE}"
: "${LDSC_BIM_DIR:?Set LDSC_BIM_DIR}"
[[ ! -e "$MAGMA_OUTPUT_ROOT" ]] || { echo "MAGMA_OUTPUT_ROOT already exists; choose a fresh directory" >&2; exit 1; }
command -v CrossMap >/dev/null || { echo "CrossMap must be available on PATH" >&2; exit 1; }
magma="${MAGMA_BIN:-${REGA_SEPARATE_ROOT}/external_inputs/magma}"
[[ -x "$magma" ]] || { echo "MAGMA executable missing: $magma" >&2; exit 1; }
# step1 prepare MAGMA input GWAS file (hg19, autochromasome)
#############################################################################################
################################ 00 Lift the GWAS data to hg19 ##############################
#############################################################################################
outdir="${MAGMA_OUTPUT_ROOT}/data/gwas/liftover_hg19"
mkdir -p "$outdir"
gwas38="$AD_GWAS_GRCH38_FILE"

chain_file="${LDSC_CHAIN_FILE:?Set LDSC_CHAIN_FILE}"

ls -lh "$chain_file"

zcat "$gwas38" | awk -v OFS='\t' '
NR==1{
 for(i=1;i<=NF;i++){
   if($i=="variant_id") snp=i
   if($i=="p_value") p=i
   if($i=="chromosome") chr=i
   if($i=="base_pair_location") bp=i
   if($i=="n_cases") cases=i
   if($i=="n_controls") controls=i
 }
 next
}
$snp!="NA" && $snp!="." && $p!="NA" && $p>0 && $p<=1 && $chr!="X" && $chr!="Y" && $chr!="MT" {
 n = $cases + $controls
 if(n > 0){
   print "chr"$chr, $bp-1, $bp, $snp, $p, n
 }
}' > "$outdir/AD_Bellenguez_GRCh38.for_liftover.bed"

#check
head "$outdir/AD_Bellenguez_GRCh38.for_liftover.bed"
wc -l "$outdir/AD_Bellenguez_GRCh38.for_liftover.bed"

# Lift GWAS positions to hg19.
CrossMap bed "$chain_file" \
 "$outdir/AD_Bellenguez_GRCh38.for_liftover.bed" \
 "$outdir/AD_Bellenguez_GRCh38_to_hg19.crossmap.bed"
  
# Inspect the lifted output.
head "$outdir/AD_Bellenguez_GRCh38_to_hg19.crossmap.bed"
wc -l "$outdir/AD_Bellenguez_GRCh38_to_hg19.crossmap.bed"

ls -lh "$outdir"

# now prepare MAGMA input
# step 1 make input version of the file
awk -v OFS='\t' '
BEGIN{
 print "SNP","CHR","BP","P","N"
}
{
 chr=$1
 gsub(/^chr/,"",chr)
 bp=$3
 print $4, chr, bp, $5, $6
}' "$outdir/AD_Bellenguez_GRCh38_to_hg19.crossmap.bed" \
> "$outdir/AD_Bellenguez_hg19_MAGMA_input.tsv"

# step 2 debuplicate
awk -v OFS='\t' '
NR==1{
 header=$0
 next
}
{
 snp=$1
 p=$4+0
 if(!(snp in best) || p < best[snp]){
   best[snp]=p
   line[snp]=$0
 }
}
END{
 print header
 for(snp in line) print line[snp]
}' "$outdir/AD_Bellenguez_hg19_MAGMA_input.tsv" \
> "$outdir/AD_Bellenguez_hg19_MAGMA_input.dedup.tsv"

# step 3 check chromosome
head "$outdir/AD_Bellenguez_hg19_MAGMA_input.dedup.tsv"
wc -l "$outdir/AD_Bellenguez_hg19_MAGMA_input.tsv"
wc -l "$outdir/AD_Bellenguez_hg19_MAGMA_input.dedup.tsv"
cut -f2 "$outdir/AD_Bellenguez_hg19_MAGMA_input.dedup.tsv" | tail -n +2 | sort -n | uniq -c

#step 4 only autosome 
awk -v OFS='\t' '
NR==1{
 print
 next
}
$2 ~ /^[0-9]+$/ && $2 >= 1 && $2 <= 22 {
 print
}' "$outdir/AD_Bellenguez_hg19_MAGMA_input.dedup.tsv" \
> "$outdir/AD_Bellenguez_hg19_MAGMA_input.dedup.autosome.tsv"

#############################################################################################
############################# 01 chromosome-wise annotation input ###########################
#############################################################################################
# Use the lifted GWAS file and supplied protein-coding gene coordinates.
gwas_hg19="${outdir}/AD_Bellenguez_hg19_MAGMA_input.dedup.autosome.tsv"
gene_loc="$MAGMA_GENE_LOC"

ref_dir="${LDSC_BIM_DIR:?Set LDSC_BIM_DIR}"

magma_input_dir="${MAGMA_OUTPUT_ROOT}/data/annotation_input"
magma_output_dir="${MAGMA_OUTPUT_ROOT}/results/chromosome"

mkdir -p "$magma_input_dir"
mkdir -p "$magma_output_dir"

for i in {1..22}; do
  ref="${ref_dir}/1000G.EUR.QC.${i}"
  bim="${ref}.bim"

  snploc="${magma_input_dir}/chr${i}.magma.input.snp.chr.pos.txt"
  pval3="${magma_input_dir}/chr${i}.magma.input.pval3.tsv"

  echo "Preparing chromosome ${i}"

  awk -v FS='\t' -v OFS='\t' -v chr="$i" '
    NR==FNR{
      keep[$2]=1
      next
    }
    FNR==1{
      for(j=1;j<=NF;j++){
        if($j=="SNP") s=j
        if($j=="CHR") c=j
        if($j=="BP") b=j
        if($j=="P") p=j
        if($j=="N") n=j
      }
      next
    }
    $c==chr && keep[$s] && $p>0 && $p<=1 && $n>0 {
      print $s,$c,$b
    }
  ' "$bim" "$gwas_hg19" > "$snploc"

  awk -v FS='\t' -v OFS='\t' -v chr="$i" '
    NR==FNR{
      keep[$2]=1
      next
    }
    FNR==1{
      for(j=1;j<=NF;j++){
        if($j=="SNP") s=j
        if($j=="CHR") c=j
        if($j=="P") p=j
        if($j=="N") n=j
      }
      next
    }
    $c==chr && keep[$s] && $p>0 && $p<=1 && $n>0 {
      print $s,$p,$n
    }
  ' "$bim" "$gwas_hg19" > "$pval3"

  echo -n "  snploc rows: "
  wc -l < "$snploc"

  echo -n "  pval rows:   "
  wc -l < "$pval3"
done


#############################################################################################
############################# 02 chromosome-wise gene-based analysis ########################
#############################################################################################

for i in {1..22}; do
  ref="${ref_dir}/1000G.EUR.QC.${i}"

  snploc="${magma_input_dir}/chr${i}.magma.input.snp.chr.pos.txt"
  pval3="${magma_input_dir}/chr${i}.magma.input.pval3.tsv"

  out_prefix="${magma_output_dir}/chr${i}"

  echo "Running MAGMA chromosome ${i}"

  "$magma" \
    --annotate window=10 \
    --snp-loc "$snploc" \
    --gene-loc "$gene_loc" \
    --out "$out_prefix"

  "$magma" \
    --bfile "$ref" \
    --pval "$pval3" use=1,2 ncol=3 \
    --gene-annot "${out_prefix}.genes.annot" \
    --out "$out_prefix"
done

#############################################################################################
###################### 03 merge chromosome-wise gene-based analysis results #################
#############################################################################################

magma_merged_output_dir="${MAGMA_OUTPUT_ROOT}/results/merged"
mkdir -p "$magma_merged_output_dir"

#Create comma-separated input lists:
chr_result_out_list=$(ls "$magma_output_dir"/chr{1..22}.genes.out | tr '\n' ',' | sed 's/,$//')
chr_result_raw_list=$(ls "$magma_output_dir"/chr{1..22}.genes.raw | tr '\n' ',' | sed 's/,$//')

#Merge
"$magma" \
  --meta genes="$chr_result_out_list" \
  --out "$magma_merged_output_dir/AD_Bellenguez_hg19_merged"

"$magma" \
  --meta raw="$chr_result_raw_list" \
  --out "$magma_merged_output_dir/AD_Bellenguez_hg19_merged"
  
# check
ls -lh "$magma_merged_output_dir"

head "$magma_merged_output_dir/AD_Bellenguez_hg19_merged.genes.out"
head "$magma_merged_output_dir/AD_Bellenguez_hg19_merged.genes.raw"

wc -l "$magma_merged_output_dir/AD_Bellenguez_hg19_merged.genes.out"
wc -l "$magma_merged_output_dir/AD_Bellenguez_hg19_merged.genes.raw"
