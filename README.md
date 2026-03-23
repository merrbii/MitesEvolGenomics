# Population Genomics Pipeline for Seven Mite Species

This document describes the bioinformatics pipelines used to perform population genomic analyses on seven mite species using PacBio HiFi long reads. The species analyzed are:

| Code | Species |
|------|---------|
| Ass  | *Atropacarus striculus* |
| Hga  | *Hermannia gibba* |
| Hls  | *Hydrozetes lacustris* |
| Hrs  | *Hypochthonius rufulus* |
| Hti  | *Hydrozetes thienemanni* |
| Ppr  | *Platynothrus peltifer* |
| Sms  | *Steganacarus magnus* |

---

## Table of Contents

1. [Read Preprocessing](#1-read-preprocessing)
2. [Read Mapping](#2-read-mapping)
3. [Quality Control and Filtering](#3-quality-control-and-filtering)
4. [Read Group Assignment](#4-read-group-assignment)
5. [Variant Calling (GATK)](#5-variant-calling-gatk)
6. [Variant Extraction and Filtration](#6-variant-extraction-and-filtration)
7. [Phasing](#7-phasing)
8. [Phased Reference Genomes and Haplotype-Resolved Reads](#8-phased-reference-genomes-and-haplotype-resolved-reads)
9. [Haplotype Trees](#9-haplotype-trees)
10. [Heterozygosity and Diversity Statistics](#10-heterozygosity-and-diversity-statistics)
11. [Linkage Disequilibrium Analysis](#11-linkage-disequilibrium-analysis)
12. [Structural Variant Analysis (GraffiTE)](#12-structural-variant-analysis-graffite)
13. [Functional Annotation of Structural Variants](#13-functional-annotation-of-structural-variants)

---

## 1. Read Preprocessing

### 1.1 Convert BAM to FASTQ

Convert PacBio HiFi BAM files to FASTQ format.

```bash
ls *.bam | parallel --jobs 7 --eta \
  'samtools fastq -@ 10 {} -0 /home/merrbii/Scratch/sexAsex/popGen/fq/{.}.fq.gz'
```

### 1.2 Remove SMRTbell Barcodes

Strip barcodes using `lima`.

```bash
for i in *.gz; do
  b=$(echo ${i} | sed 's/.fq.gz//g')
  echo "lima -j 1 ${i} \
    /home/merrbii/Scratch/sexAsex/popGen/fq/SMRTbell_Barcoded_Adapter_Plate_3.0_bc2001-bc2096.fasta \
    ${b}.bcFree.fq.gz"
done | parallel --jobs 40 --eta
```

### 1.3 Remove PCR/Optical Duplicates

Use `pbmarkdup` to remove duplicate reads.

```bash
for i in *.gz; do
  base=$(echo ${i} | sed 's/.fq.gz//g')
  echo "pbmarkdup -j 5 ${i} ${base}.uniq.fq \
    --log-file ${base}.uniq.txt \
    --dup-file ${base}.dups.fa \
    --ignore-read-names; \
    bgzip ${base}.uniq.fq; \
    bgzip ${base}.dups.fa; \
    rm ${base}.uniq.fq ${base}.dups.fa"
done | parallel --jobs 16 --eta
```

---

## 2. Read Mapping

### 2.1 Build Winnowmap Repetitive k-mer Index

Generate repetitive k-mer lists for each reference genome using `meryl`.

```bash
for i in Ass Hga Hls Hrs Hti Ppr Sms; do
  meryl count k=15 threads=50 output merylDB \
    ~/Scratch/sexAsex/hyphy/genomes/${i}.fasta
  meryl print greater-than distinct=0.9998 merylDB > ${i}.repetitive_k15.txt
done
```

### 2.2 Map Reads with Winnowmap

Align HiFi reads to species-specific reference genomes using `winnowmap`, then sort with `samtools`.

```bash
for i in *gz; do
  b=$(echo ${i} | sed 's/.fq.gz//g')
  sp=$(echo ${i} | cut -f1 -d".")
  echo "winnowmap -W ${sp}.repetitive_k15.txt --MD -Y -ax map-pb -t 15 \
    ~/Scratch/sexAsex/hyphy/genomes/${sp}.fasta ${i} \
    | samtools sort -O bam -o ../bams/${b}.bam -"
done | parallel --jobs 2 --eta
```

---

## 3. Quality Control and Filtering

### 3.1 Index BAM Files

```bash
ls *.bam | parallel --jobs 30 --eta 'samtools index {}'
```

### 3.2 Generate Read Quality Plots

Use `NanoPlot` to visualize alignment quality metrics (downsampled to 10,000 reads).

```bash
ls *.bam | parallel --jobs 10 --eta \
  "NanoPlot -t 5 --color yellow --bam {} --downsample 10000 \
   -o nanoplots/{.}.bamplots_downsampled"
```

### 3.3 Filter Alignments

Retain only primary alignments with mapping quality ≥ 30 and read length > 5 kb.

```bash
ls *.bam | parallel --jobs 20 --eta \
  "samtools view -e 'length(seq)>4999' -F4 -q 30 -b {} -o {.}.filtered.bam"
```

### 3.4 Index Filtered BAMs

```bash
ls *.filtered.bam | parallel --jobs 60 --eta 'samtools index {}'
```

### 3.5 Calculate Coverage Statistics

Compute per-base mean coverage and standard deviation for each BAM file.

```bash
ls *.bam | parallel --jobs 20 --eta \
  "samtools depth -a {} | awk '{sum+=\$3; sumsq+=\$3*\$3} END { \
    print \"{}\", sum/NR, sqrt(sumsq/NR - (sum/NR)**2)}'" >> coverage.txt
```

### 3.6 Count Total Mapped Reads

```bash
for i in *.bam; do
  reads=$(samtools view -c ${i})
  echo ${i} ${reads} >> numberOfReads.txt
done
```

### 3.7 Merge Coverage and Read Count Statistics

```bash
paste \
  <(cat coverage.txt | sort | tr -s " " | tr " " "\t") \
  <(cat numberOfReads.txt | sort | tr -s " " | tr " " "\t") \
  > coverage.stats.txt
```

### 3.8 Subsample High-Coverage Samples to ~40×

For samples with mean coverage > 41×, downsample to 40× (two independent subsamples per sample for reproducibility testing).

```bash
cat coverage.stats.txt | sort -V -k2,2 | \
  awk '{ if($2>41) print $0"\t"(($4*40)/$2)/$4 }' | \
  while read sample a b c cov; do
    for i in {1..2}; do
      out=$(echo ${sample} | sed 's/.bam//g')
      echo "samtools view -@ 5 -b --subsample ${cov} --subsample-seed $RANDOM ${sample} \
        > subsampled/${out}.sub.${i}.bam; \
        samtools index -@ 2 subsampled/${out}.sub.${i}.bam"
    done
  done > subsampling.cmd.txt

parallel --jobs 10 --eta < subsampling.cmd.txt
```

> **Note:** Subsampling was performed in duplicate (sub.1 and sub.2) to verify that downsampling does not introduce bias. PCA confirmed that subsampled replicates superimpose, so one replicate per sample was retained for downstream analyses.

---

## 4. Read Group Assignment

### 4.1 Prepare Read Group Information

Create a tab-separated file mapping BAM filenames to sample names.

```bash
paste <(ls *.bam) <(ls *.bam | sed 's/.bam//g') > rgid.txt
```

> Manually edit `rgid.txt` to add correct sample names in the third column if needed.

### 4.2 Add Read Groups with Picard

```bash
parallel -a rgid.txt --colsep '\t' --jobs 25 --eta \
  java -jar ~/Software/picard.jar AddOrReplaceReadGroups \
    I={1} O={2}.R.bam \
    RGID={3} RGLB={3} RGPL=pb RGPU={3} RGSM={3}
```

### 4.3 Index Read-Group-Annotated BAMs

```bash
ls *R.bam | parallel --progress --eta --jobs 30 'samtools index {}'
```

---

## 5. Variant Calling (GATK)

### 5.1 Per-Sample Variant Calling (HaplotypeCaller in GVCF Mode)

Call variants per sample using GATK `HaplotypeCaller` in GVCF mode (diploid).

```bash
for i in *R.bam; do
  sp=$(echo ${i} | cut -f1 -d".")
  b=$(echo ${i} | sed 's/.bam//g')
  echo "~/programs/jdk-17.0.12/bin/java -jar \
    ~/programs/gatk-4.6.0.0/gatk-package-4.6.0.0-local.jar HaplotypeCaller \
    -R /home/merrbii/Scratch/sexAsex/hyphy/genomes/${sp}.fasta \
    -I ${i} -O ${b}.g.vcf.gz -ERC GVCF -ploidy 2"
done | parallel --jobs 25 --eta
```

### 5.2 Create Per-Species GVCF Lists

```bash
for sp in Ass Hga Hls Hrs Hti Ppr Sms; do
  ls ${sp}*.gz.tbi | sed 's/.tbi//g' > ${sp}.gVCF.list
done
```

### 5.3 Joint Genotyping (GenotypeGVCFs)

Jointly genotype all samples per species, emitting all sites (variant and invariant).

```bash
for sp in Ass Hga Hls Hrs Hti Ppr Sms; do
  echo "~/programs/jdk-17.0.12/bin/java \
    -DGATK_STACKTRACE_ON_USER_EXCEPTION=true \
    -Djava.io.tmpdir=/home/merrbii/Scratch/temp/ \
    -jar ~/programs/gatk-4.6.0.0/gatk-package-4.6.0.0-local.jar GenotypeGVCFs \
    -R /home/merrbii/Scratch/sexAsex/hyphy/genomes/${sp}.fasta \
    -V ${sp}.all.g.vcf.gz \
    -O VCFs/${sp}.all.vcf.gz --all-sites"
done | parallel --jobs 6 --eta --tmpdir ~/Scratch/temp/
```

---

## 6. Variant Extraction and Filtration

### 6.1 Extract Biallelic SNPs

```bash
# Extract SNPs from the all-sites VCF
for sp in Ass Hga Hls Hrs Hti Ppr Sms; do
  echo "~/programs/jdk-17.0.12/bin/java \
    -DGATK_STACKTRACE_ON_USER_EXCEPTION=true \
    -Djava.io.tmpdir=/home/merrbii/Scratch/temp/ \
    -jar ~/programs/gatk-4.6.0.0/gatk-package-4.6.0.0-local.jar SelectVariants \
    -V ${sp}.all.vcf.gz --select-type-to-include SNP \
    -O variants/${sp}.varSNPs.vcf.gz"
done | parallel --jobs 9 --eta --tmpdir ~/Scratch/temp/

# Retain only biallelic SNPs
for sp in Ass Hga Hls Hrs Hti Ppr Sms; do
  echo "bcftools view -m2 -M2 ${sp}.varSNPs.vcf.gz -Oz -o ${sp}.varSNPs.biall.vcf.gz"
done | parallel --jobs 10 --eta
```

### 6.2 Extract QC Metrics for Threshold Selection

```bash
for sp in Ass Hga Hls Hrs Hti Ppr Sms; do
  echo "bcftools query \
    -f '%INFO/DP\t%INFO/QD\t%INFO/FS\t%INFO/SOR\t%INFO/MQ\t%INFO/MQRankSum\t%INFO/ReadPosRankSum\n' \
    ${sp}.varSNPs.biall.vcf.gz > ${sp}.varSNPs.biall.info_QCs.txt"
done | parallel --jobs 9 --eta
```

> **Note:** Visualize these distributions in R to determine appropriate filtering thresholds. In this analysis, GATK hard-filtering recommendations were adopted as they were consistent with the observed distributions.

### 6.3 Apply GATK Hard Filters

```bash
for sp in Ass Hga Hls Hrs Hti Ppr Sms; do
  echo "~/programs/jdk-17.0.12/bin/java \
    -DGATK_STACKTRACE_ON_USER_EXCEPTION=true \
    -Djava.io.tmpdir=/home/merrbii/Scratch/temp/ \
    -jar ~/programs/gatk-4.6.0.0/gatk-package-4.6.0.0-local.jar VariantFiltration \
    -V ${sp}.varSNPs.biall.vcf.gz \
    -filter \"QD < 2.0\" --filter-name \"QD2\" \
    -filter \"SOR > 3.0\" --filter-name \"SOR3\" \
    -filter \"FS > 60.0\" --filter-name \"FS60\" \
    -filter \"MQ < 40.0\" --filter-name \"MQ40\" \
    -filter \"MQRankSum < -5.0\" --filter-name \"MQRankSum-5\" \
    -filter \"ReadPosRankSum < -5.0\" --filter-name \"ReadPosRankSum-5\" \
    -O ${sp}.varSNPs.biall_F.vcf.gz"
done | parallel --jobs 10 --eta --tmpdir ~/Scratch/temp/
```

### 6.4 Select PASS Variants Only

```bash
for sp in Ass Hga Hls Hrs Hti Ppr Sms; do
  echo "~/programs/jdk-17.0.12/bin/java \
    -DGATK_STACKTRACE_ON_USER_EXCEPTION=true \
    -Djava.io.tmpdir=/home/merrbii/Scratch/temp/ \
    -jar ~/programs/gatk-4.6.0.0/gatk-package-4.6.0.0-local.jar SelectVariants \
    -V ${sp}.varSNPs.biall_F.vcf.gz \
    -select 'vc.isNotFiltered()' \
    -O ${sp}.varSNPs.biallpass.vcf.gz"
done | parallel --jobs 6 --eta --tmpdir ~/Scratch/temp/
```

### 6.5 Extract VCFtools QC Metrics

```bash
for sp in Ass Hga Hls Hrs Hti Ppr Sms; do
  echo "/home/merrbii/Scratch/sexAsex/popGen/bams/gVCFs/VCFs/src/get.VCF.QCs.sh \
    ${sp}.varSNPs.biallpass.vcf.gz"
done | parallel --jobs 10 --eta
```

Fix allele frequency header:

```bash
for sp in Ass Hga Hls Hrs Hti Ppr Sms; do
  sed "s/{FREQ}/$(printf 'Ref_Allele\tAlt_Allele')/g" \
    ${sp}.varSNPs.biallpass_QC/${sp}.varSNPs.biallpass.freq.frq > tmp.txt \
    && mv tmp.txt ${sp}.varSNPs.biallpass_QC/${sp}.varSNPs.biallpass.freq.frq
done
```

Check per-individual missingness:

```bash
for i in *_QC/*varSNPs.biallpass.missing.indv.imiss; do
  echo $i
  cat $i
done
```

> **Note:** Visualize all QC distributions in R to decide on per-species depth and missingness thresholds.

### 6.6 Apply Depth, Missingness, and MAC Filters (VCFtools)

Filtering thresholds are stored in a tab-separated file `vcftools.filters.txt` (columns: species, min-meanDP, max-meanDP). Individuals to be excluded are listed in `remove.indv.txt`.

```bash
parallel -a vcftools.filters.txt --colsep '\t' --jobs 10 --eta \
  vcftools --gzvcf {1}.varSNPs.biallpass.vcf.gz \
    --min-alleles 2 --max-alleles 2 \
    --recode --recode-INFO-all \
    --out {1}.varSNPs.biallpass.filtered \
    --max-meanDP {3} --min-meanDP {2} \
    --minDP 5 --max-missing 0.8 --mac 1 \
    --remove remove.indv.txt \
    "2>" {1}.varSNPs.biallpass.filtered.log
```

Rename and compress:

```bash
for i in *.recode.vcf; do
  b=$(echo $i | sed 's/.recode//g')
  mv $i $b
done

ls *vcf | parallel --jobs 9 --eta 'bgzip {}'
```

### 6.7 Extract and Filter Invariant Sites

Extract invariant (monomorphic) sites separately for downstream diversity calculations.

```bash
# Extract invariant sites
for sp in Ass Hga Hls Hrs Hti Ppr Sms; do
  echo "~/programs/jdk-17.0.12/bin/java \
    -DGATK_STACKTRACE_ON_USER_EXCEPTION=true \
    -Djava.io.tmpdir=/home/merrbii/Scratch/temp/ \
    -jar ~/programs/gatk-4.6.0.0/gatk-package-4.6.0.0-local.jar SelectVariants \
    -V ${sp}.all.vcf.gz --select-type-to-include NO_VARIATION \
    -O invariants/${sp}.invarSNPs.vcf.gz"
done | parallel --jobs 9 --eta --tmpdir ~/Scratch/temp/

# Filter invariant sites (depth and missingness)
parallel -a ../variants/vcftools.filters.txt --colsep '\t' --jobs 9 --eta \
  vcftools --gzvcf {1}.invarSNPs.vcf.gz \
    --recode --recode-INFO-all \
    --out {1}.invarSNPs.filtered \
    --minDP 5 --max-missing 0.8 \
    "2>" {1}.invarSNPs.biallpass.filtered.log

# Rename and compress
for i in *.recode.vcf; do
  b=$(echo $i | sed 's/.recode//g')
  mv $i $b
done

ls *vcf | parallel --jobs 10 --eta 'bgzip {}'
```

---

## 7. Phasing

### 7.1 Phase Filtered Variant Sites with WhatsHap

Use long-read BAMs to perform read-backed phasing of biallelic SNPs.

```bash
for sp in Ass Hga Hls Hrs Hti Ppr Sms; do
  echo "whatshap phase \
    -o ${sp}.varSNPs.biallpass.filtered.phased.vcf.gz \
    --reference=/home/merrbii/Scratch/sexAsex/hyphy/genomes/${sp}.fasta \
    ${sp}.varSNPs.biallpass.filtered.vcf.gz \
    $(ls ${sp}*.bam | tr '\n' ' ')"
done | parallel --jobs 3 --eta
```

### 7.2 Generate Homozygous-Only VCFs for Phylogenetic Analysis

Mask heterozygous genotypes as missing for tree-building purposes.

```bash
for i in *biallpass.vcf.gz; do
  b=$(echo ${i} | sed 's/.vcf.gz//g')
  echo "bcftools filter -e 'FMT/GT=\"het\"' --set-GTs . ${i} \
    -Oz -o ${b}.hom.vcf.gz"
done | parallel --jobs 9 --eta
```

### 7.3 QC and Filter Homozygous VCFs

```bash
# Extract QC metrics
for sp in Ass Hga Hls Hrs Hti Ppr Sms; do
  echo "/home/merrbii/Scratch/sexAsex/popGen/bams/gVCFs/VCFs/src/get.VCF.QCs.sh \
    ${sp}.varSNPs.biallpass.hom.vcf.gz"
done | parallel --jobs 10 --eta

# Fix frequency header
for sp in Ass Hga Hls Hrs Hti Ppr Sms; do
  sed "s/{FREQ}/$(printf 'Ref_Allele\tAlt_Allele')/g" \
    ${sp}.varSNPs.biallpass.hom_QC/${sp}.varSNPs.biallpass.hom.freq.frq > tmp.txt \
    && mv tmp.txt ${sp}.varSNPs.biallpass.hom_QC/${sp}.varSNPs.biallpass.hom.freq.frq
done

# Apply filters (thresholds from R visualization; stored in vcftools.filters.txt)
parallel -a vcftools.filters.txt --colsep '\t' --jobs 9 --eta \
  vcftools --gzvcf {1}.varSNPs.biallpass.hom.vcf.gz \
    --min-alleles 2 --max-alleles 2 \
    --recode --recode-INFO-all \
    --out {1}.varSNPs.biallpass.hom.filtered \
    --max-meanDP {3} --min-meanDP {2} \
    --minDP 5 --max-missing 0.75 --maf 0.05 \
    "2>" {1}.varSNPs.biallpass.filtered.log
```

---

## 8. Phased Reference Genomes and Haplotype-Resolved Reads

### 8.1 Evaluate Phasing Statistics

The phasing strategy first evaluated phasing statistics using the reference samples to assess the size and continuity of haplotype blocks. Variants were then incorporated into the reference genome to generate two phased haplotype assemblies (haplotype 1 and haplotype 2). Sequencing reads were mapped competitively against both haplotypes, and reads were subsequently assigned to the haplotype to which they mapped best, thereby separating the data into haplotype-specific read sets. Check haplotype block sizes (phase sets) for each reference sample.

```bash
for i in *.gz; do
  sm=$(bcftools query -l ${i} | grep "ref")
  echo "whatshap stats ${i} \
    --gtf ${sm}.blocks.gtf \
    --block-list ${sm}.blocks.gtf \
    --tsv ${sm}.stats.tsv \
    --sample ${sm}"
done | parallel --jobs 9 --eta
```

### 8.2 Generate Haplotype-Resolved Reference Genomes

Introduce phased variants into the reference to produce two haplotype assemblies per species.

```bash
for i in *.gz; do
  for hap in 1 2; do
    sm=$(bcftools query -l ${i} | grep "ref")
    sp=$(echo ${i} | cut -f1 -d".")
    echo "bcftools consensus -s ${sm} \
      -f /home/merrbii/Scratch/sexAsex/hyphy/genomes/${sp}.fasta \
      -H ${hap} \
      -o /home/merrbii/Scratch/sexAsex/popGen/2ndRound/phasedRefs/${sp}.hap${hap}.fasta ${i}"
  done
done | parallel --jobs 9 --eta
```

### 8.3 Build k-mer Indices for Phased References

```bash
for i in *.fasta; do
  b=$(echo ${i} | cut -f1,2 -d".")
  echo $b
  meryl count k=15 threads=50 output merylDB ${i}
  meryl print greater-than distinct=0.9998 merylDB > ${b}.repetitive_k15.txt
done
```

### 8.4 Competitive Mapping to Both Haplotypes

Map each sample's reads to both haplotype references.

```bash
for i in *gz; do
  for hap in 1 2; do
    b=$(echo ${i} | sed 's/.fq.gz//g')
    sp=$(echo ${i} | cut -f1 -d".")
    echo "winnowmap \
      -W /home/merrbii/Scratch/sexAsex/popGen/2ndRound/phasedRefs/${sp}.hap${hap}.repetitive_k15.txt \
      -ax map-pb -t 15 \
      /home/merrbii/Scratch/sexAsex/popGen/2ndRound/phasedRefs/${sp}.hap${hap}.fasta ${i} \
      | samtools sort -O bam -o ../phasedBams/${b}.hap${hap}.bam -"
  done
done > competitiveMapping.cmd.txt

cat competitiveMapping.cmd.txt | parallel --jobs 2 --eta
```

### 8.5 Compute Weighted Alignment Score (WAS) and Assign Reads to Haplotypes

```bash
# Compute WAS for each BAM
for i in *bam; do
  b=$(echo $i | sed 's/.bam/.tsv/g')
  echo ./get_WAS.sh ${i} ${b} 10
done | parallel --jobs 10 --eta

# Assign reads to haplotypes using R script
Rscript process_reads.R
```

> get_WAS.sh can be found [here](https://github.com/hewm2008/PopLDdecay) and process_reads.R can be found [here](https://github.com/hewm2008/PopLDdecay).

### 8.6 Extract Haplotype-Assigned Reads

```bash
for i in *.fq.gz; do
  for h in hap1 hap2; do
    b=$(echo $i | sed 's/.fq.gz//g')
    echo "seqtk subseq ${i} ${b}.${h}.txt | gzip > ${b}.${h}.fq.gz"
  done
done | parallel --jobs 25 --eta
```

> After extraction, repeat the mapping, variant calling, and filtering pipeline (Sections 2–6) on the haplotype-resolved reads.

---

## 9. Haplotype Trees

### 9.1 Prepare Phased Consensus Sequences

Link the final filtered VCFs (variant + invariant sites) from the phased read analysis:

```bash
ln -s /home/merrbii/Scratch/sexAsex/popGen/phasedCalls/phasedBams/gVCFs/VCFs/finalVCFs/*.filtered.vcf.gz* .
```

### 9.2 Subset to Large Haplotype Blocks (≥ 100 kb)

Extract sites falling within haplotype blocks of at least 100 kb.

```bash
for i in Ass Hga Hls Hrs Hti Ppr Sms; do
  echo "bcftools view \
    -R ../haploblocks/${i}.blocks.atleast100kb.bed \
    -Oz -o ${i}.100kb.haploblocks.var.invar.filtered.vcf.gz \
    ${i}.var.invar.filtered.vcf.gz"
done | parallel --jobs 9 --eta

ls *100kb.haploblocks.var.invar.filtered.vcf.gz | parallel --jobs 10 tabix {}
```

### 9.3 Generate Per-Sample FASTA Consensus Sequences

```bash
for i in *.100kb.*.vcf.gz; do
  for sm in $(bcftools query -l ${i} | tr "\n" " "); do
    sp=$(echo ${i} | cut -f1 -d".")
    echo "bcftools consensus -s ${sm} -p ${sm} -M N -a N -I \
      -f /home/merrbii/Scratch/sexAsex/hyphy/genomes/${sp}.fasta \
      -o /home/merrbii/Scratch/sexAsex/popGen/MEanalysis/consensusSeqs/${sm}.fasta ${i}"
  done
done | parallel --jobs 50 --eta
```

### 9.4 Organize and Split by Chromosome

Concatenate all sample FASTAs per species, linearize, and split by chromosome.

```bash
# Concatenate and linearize
for i in Ass Hga Hls Hrs Hti Ppr Sms; do
  mkdir ${i}
  cat ${i}*.fasta | awk '{
    if(NR==1) {print $0}
    else {if($0 ~ /^>/) {print "\n"$0} else {printf $0}}
  }' > ${i}/${i}.ni.fasta
done

# Get chromosome names
for i in Ass Hga Hls Hrs Hti Ppr Sms; do
  grep ">" /home/merrbii/Scratch/sexAsex/hyphy/genomes/${i}.fasta \
    | sed 's/>//g' > ${i}/${i}.seqName.txt
done

# Split multi-FASTA by chromosome
for i in Ass Hga Hls Hrs Hti Ppr Sms; do
  (cd "$i" && parallel --will-cite --no-notice -j 50 \
    'grep -A 1 --no-group-separator "{}\b" '"${i}.ni.fasta"' > '"${i}"'_{}.fasta' \
    :::: "${i}.seqName.txt")
done
```

### 9.5 Window-Based Alignment Splitting (100 kb Windows)

```bash
for i in Ass Hga Hls Hrs Hti Ppr Sms; do
  mv ${i}/${i}.ni.fasta ${i}/${i}.ni.fa
done

for i in Ass Hga Hls Hrs Hti Ppr Sms; do
  (cd "$i" && parallel --will-cite --no-notice -j 40 \
    'msa_split "{}" --windows 100000,0 --out-root {.}' ::: *.fasta)
done
```

### 9.6 Build Maximum Likelihood Trees (IQ-TREE2)

Infer trees per 100 kb window using ModelFinder and 1000 ultrafast bootstrap replicates.

```bash
for dir in Ass Hga Hls Hrs Hti Ppr Sms; do
  (cd "$dir" && parallel --will-cite --no-notice -j 10 \
    'iqtree2 -nt 8 -s {} -m MFP -B 1000 --redo' ::: *.fa)
done
```

### 9.7 Post-Processing Trees

Move failed runs aside and clean tip labels.

```bash
# Identify and move failed runs
for i in Ass Hga Hls Hrs Hti Ppr Sms; do
  (cd "$i" && mkdir -p failed && \
   for f in *.log; do
     grep -q "ERROR" "$f" && mv "${f%.log}"* failed/
   done)
done

# Strip chromosome prefix from tip labels
for i in *.fa; do
  sp=$(echo $i | cut -f1 -d"_")
  b=$(echo $i | cut -f1 -d"." | cut -f2- -d"_")
  echo "cat ${i}.contree | sed 's/${b}//g' > ${i}.contree2"
done | parallel --jobs 50

# Collect cleaned trees
for i in Ass Hga Hls Hrs Hti Ppr Sms; do
  mkdir -p trees/${i}
  cp ${i}/*.contree2 trees/${i}/
done
```

---

## 10. Heterozygosity and Diversity Statistics

### 10.1 Concatenate Variant and Invariant Sites

Merge the filtered variant and invariant VCFs to produce an all-sites VCF required for unbiased diversity estimation.

```bash
for i in *gz; do
  sp=$(echo $i | cut -f1 -d".")
  echo "bcftools concat --allow-overlaps \
    invariants/${sp}.invarSNPs.filtered.recode.vcf.gz \
    variants/${sp}.varSNPs.biallpass.filtered.vcf.gz \
    | bcftools view -e 'GT=\".\"' -Oz \
    -o finalVCFs/${sp}.var.invar.filtered.vcf.gz"
done | parallel --jobs 9 --eta
```

### 10.2 Generate Genotype Files for `genomics_general`
All the scripts can be found [here](https://github.com/simonhmartin/genomics_general)
```bash
for i in *gz; do
  sp=$(echo $i | cut -f1 -d".")
  echo "python3 ~/programs/genomics_general/VCF_processing/parseVCF.py \
    -i ${i} -o ../geno/${sp}.var.invar.filtered.geno.gz"
done | parallel --jobs 9 --eta
```

### 10.3 Create Population Files

```bash
for i in *gz; do
  sp=$(echo $i | cut -f1 -d".")
  bcftools query -l $i | awk -v var="$sp" '{print $0"\t"var}' > ../geno/${sp}.pop.txt
done
```

### 10.4 Compute Windowed Diversity (π) — Unphased

Site-based windows of 25,000 sites with a minimum of 10,000 sites.

```bash
for i in *gz; do
  sp=$(echo $i | cut -f1 -d".")
  echo "python3 ~/programs/genomics_general/popgenWindows.py \
    --windType sites -w 25000 -m 10000 --roundTo 6 \
    -g ${sp}.var.invar.filtered.geno.gz \
    -o ${sp}.var.invar.filtered.csv.gz \
    -f phased -T 4 \
    -p ${sp} --popsFile ${sp}.pop.txt --writeFailedWindows"
done | parallel --jobs 3 --eta
```

### 10.5 Compute Windowed Diversity (π) — Phased (Between Haplotypes)

```bash
for i in *gz; do
  sp=$(echo $i | cut -f1 -d".")
  echo "python3 ~/programs/genomics_general/popgenWindows.py \
    --windType sites -w 25000 -m 10000 --roundTo 6 \
    -g ${sp}.var.invar.filtered.geno.gz \
    -o ${sp}.var.invar.filtered.csv.gz \
    -f phased -T 5 \
    -p ${sp}1 -p ${sp}2 \
    --popsFile ${sp}.pop.txt --writeFailedWindows"
done | parallel --jobs 9 --eta
```

### 10.6 Per-Individual Heterozygosity

```bash
for i in *.gz; do
  sp=$(echo $i | cut -f1 -d".")
  echo "./get_het_counts.sh ${i} > ../het/${sp}.het.txt"
done | parallel --jobs 3 --eta
```

### 10.7 Tajima's D

Calculate genome-wide Tajima's D in 25 kb windows.

```bash
parallel -a <(ls *gz | cut -f1 -d".") --colsep '\t' --jobs 10 --eta \
  vcftools --gzvcf {1}.var.invar.filtered.vcf.gz \
    --TajimaD 25000 \
    --out ../het/{1}.var.invar.filtered.tajD \
    "2>" ../het/{1}.var.invar.filtered.tajD.log
```

### 10.8 Site Frequency Spectrum (SFS)

Compute SFS from unphased variant-only VCFs.

```bash
for sp in Ass Hga Hls Hrs Hti Ppr Sms; do
  ../src/get_sfs_profiles.sh \
    ../unphasedVCFs/filteredVCFs/${sp}.varSNPs.biallpass.filtered.vcf.gz ${sp}
done
```

> Continue downstream SFS visualization and analysis in R using [get_sfs_profiles.sh](https://github.com/hewm2008/PopLDdecay).

---

## 11. Linkage Disequilibrium Analysis

LD decay is estimated using [PopLDdecay](https://github.com/hewm2008/PopLDdecay).

### 11.1 Calculate LD Decay at MAF ≥ 0.2

```bash
# MAF ≥ 0.2
for i in $(ls *vcf.gz | cut -f1 -d"."); do
  echo "~/programs/PopLDdecay-3.43/bin/PopLDdecay \
    -InVCF ${i}.varSNPs.biallpass.filtered.vcf.gz \
    -OutStat ${i}.LDdecay.maf0.2 -MAF 0.2"
done | parallel --jobs 10 --eta
```

### 11.2 Generate Multi-Species LD Decay Plots

```bash
# MAF 0.2
paste <(ls *maf0.2.stat.gz) <(ls *maf0.2.stat.gz | cut -f1 -d".") > multi.list.maf0.2
perl ~/programs/PopLDdecay-3.43/bin/Plot_MultiPop.pl \
  -InList multi.list.maf0.2 -output Multi.LDdecay.maf0.2
```

---

## 12. Structural Variant Analysis (GraffiTE)

Structural variants (SVs) and transposable element (TE) insertions are detected and genotyped using [GraffiTE](https://github.com/cgroza/GraffiTE), which internally uses Sniffles2 and a TE consensus library.

### 12.1 Environment Setup

```bash
# Install Nextflow
curl -s https://get.nextflow.io | bash
chmod +x nextflow
mkdir -p $HOME/.local/bin/
mv nextflow $HOME/.local/bin/

# Ensure Apptainer/Singularity is in PATH
PATH=$PATH:/home/merrbii/programs/apptainer/bin
unset JAVA_HOME
```

### 12.2 Prepare Input File Lists

```bash
# Long read input lists
for i in Ass Hga Hls Hrs Hti Ppr Sms; do
  ls cleanReads/${i}*fq.gz > ${i}.longreads.txt
done

for i in Ass Hga Hls Hrs Hti Ppr Sms; do
  paste ${i}.longreads.txt \
    <(cat ${i}.longreads.txt | cut -f2 -d"/" | sed 's/.uniq.bcFree.fq.gz//g') \
    | tr "\t" "," > tmp.txt && mv tmp.txt ${i}.longreads.txt
done

for i in Ass Hga Hls Hrs Hti Ppr Sms; do
  awk '{print $0",hifi"}' ${i}.longreads.txt > tmp.txt \
    && mv tmp.txt ${i}.longreads.txt
done

# BAM input lists
for i in Ass Hga Hls Hrs Hti Ppr Sms; do
  ls /home/merrbii/Scratch/sexAsex/popGen/bams/${i}*bam > ${i}.bamlist.txt
done

for i in Ass Hga Hls Hrs Hti Ppr Sms; do
  paste ${i}.bamlist.txt \
    <(cat ${i}.bamlist.txt | cut -f8 -d"/" | sed 's/.uniq.bcFree.bam//g') \
    | tr "\t" "," > tmp.txt && mv tmp.txt ${i}.bamlist.txt
done
```

### 12.3 Run GraffiTE

```bash
for i in Ass Hga Hls Hrs Hti Ppr Sms; do
  nextflow run /home/merrbii/programs/GraffiTE/main.nf \
    --reference ~/Scratch/sexAsex/hyphy/genomes/${i}.fasta \
    --bams ${i}.bamlist.txt \
    --TE_library Acari.deNovo-repeats.fa \
    -with-singularity ~/programs/graffite_latest.sif \
    --cores 40 \
    --longreads false \
    --graph_method graphaligner \
    --reads ${i}.longreads.txt \
    --out ${i}
  rm -rf work/
done
```

### 12.4 Annotate and Tabulate Sniffles2 SV Calls

```bash
# Add allele count/frequency annotations
for i in Ass Hga Hls Hrs Hti Ppr Sms; do
  bcftools +fill-tags ${i}.sniffles2_variants.vcf \
    -Oz -o ${i}.sniffles2_variants.annotated.vcf.gz \
    -- -t AN,AC,AF,F_MISSING,MAF,NS
done

# Extract tabular summary
for i in Ass Hga Hls Hrs Hti Ppr Sms; do
  bcftools query \
    -f'%CHROM\t%POS\t%ID\t%QUAL\t%INFO/SVTYPE\t%INFO/SVLEN\t%INFO/AC\t%INFO/F_MISSING\t%INFO/NS\t%INFO/AN\t%INFO/AF\t%INFO/MAF[\t%DR\t%DV]\n' \
    ${i}.sniffles2_variants.annotated.vcf.gz | \
    awk 'BEGIN{OFS="\t"}{
      n=(NF-12)/2;
      for(j=1;j<=12;j++) printf "%s\t",$j;
      for(k=1;k<=n;k++){
        dr=$(12+2*k-1); dv=$(12+2*k);
        dr=(dr==".")?0:dr; dv=(dv==".")?0:dv;
        printf "%s",dr+dv;
        if(k<n) printf "\t"
      }
      print ""
    }' > ${i}.sniffles2_variants.table
done
```

### 12.5 Process GraffiTE TE-Associated SVs

```bash
# Standardize non-biallelic genotypes to missing
for i in Ass Hga Hls Hrs Hti Ppr Sms; do
  bcftools +setGT ${i}.GraffiTE.merged.genotypes.vcf.gz \
    -- -t q -n './.' \
    -e 'GT="1/1" | GT="0/0" | GT="0/1" | GT="1/0"' \
    | bcftools view -Oz -o ${i}.GraffiTE.merged.genotypes.biall.vcf.gz
done

# Add annotations
for i in Ass Hga Hls Hrs Hti Ppr Sms; do
  bcftools +fill-tags ${i}.GraffiTE.merged.genotypes.biall.vcf.gz \
    -Oz -o ${i}.GraffiTE.merged.genotypes.annotated.biall.vcf.gz \
    -- -t AN,AC,AF,F_MISSING,MAF,NS
done

# Extract tabular summary
for i in Ass Hga Hls Hrs Hti Ppr Sms; do
  bcftools query \
    -f'%CHROM\t%POS\t%ID\t%QUAL\t%INFO/SVTYPE\t%INFO/SVLEN\t%INFO/AC\t%INFO/F_MISSING\t%INFO/NS\t%INFO/AN\t%INFO/AF\t%INFO/MAF\t%INFO/DP\n' \
    ${i}.GraffiTE.merged.genotypes.annotated.biall.vcf.gz \
    > ${i}.GraffiTE.merged.genotypes.annotated.biall.table
done
```

> **Note:** Visualize depth, MAF, and missingness distributions in R to determine appropriate filtering thresholds.

### 12.6 Filter SVs

Filter Sniffles2 calls and GraffiTE TE calls using curated ID lists generated from R-based QC visualization.

```bash
# Filter Sniffles2 SVs
for i in Ass Hga Hls Hrs Hti Ppr Sms; do
  bcftools view --include ID==@unique_ids_${i}.txt \
    ${i}.sniffles2_variants.annotated.vcf.gz \
    -Oz -o ${i}.sniffles2_variants.annotated.filtered.vcf.gz
done

# Filter GraffiTE TE SVs
for i in Ass Hga Hls Hrs Hti Ppr Sms; do
  bcftools view --include ID==@unique_ids_${i}.txt \
    ${i}.GraffiTE.merged.genotypes.annotated.biall.vcf.gz \
    -Oz -o ${i}.GraffiTE.merged.genotypes.annotated.biall.filtered.vcf.gz
done
```

### 12.7 Calculate Heterozygosity for TE and Non-TE SVs

```bash
# TE-associated SVs
for i in *gz; do
  b=$(echo ${i} | cut -f1 -d".")
  echo ${b}
  ../calculateHet.sh ${i} > ../het/${b}.TEs.het.txt
done

# Non-TE SVs
for i in *gz; do
  b=$(echo ${i} | cut -f1 -d".")
  echo ${b}
  ../calculateHet.sh ${i} > ../het/${b}.nonTEs.het.txt
done
```

---

## 13. Functional Annotation of Structural Variants

### 13.1 Run InterProScan on Predicted Proteins

```bash
# Link and clean protein sequences (remove stop codon asterisks)
ln -s /home/merrbii/Scratch/sexAsex/hyphy/pep/tmp/*pep.fa .

for i in *.fa; do
  b=$(echo $i | cut -f 1,2 -d".")
  cat ${i} | sed 's/*//g' > ${b}.asterixFree.fa
done

conda deactivate

for i in *asterixFree.fa; do
  b=$(echo $i | cut -f1 -d".")
  mkdir $b
  ~/Scratch/interpro/interproscan-5.69-101.0/interproscan.sh \
    -i ${i} -d ${b} -cpu 25 -T temp --goterms -dp
done
```

### 13.2 Extract SV Coordinates

```bash
for i in *ariants.annotated.filtered.vcf.gz; do
  b=$(echo ${i} | sed 's/_variants.annotated.filtered.vcf.gz//g')
  bcftools query -f '%CHROM\t%POS\t%END\t%ID\n' $i \
    > SVannotation/tmp/$b.CHR_POS_END_ID.table
done

for i in *table; do
  sortBed -i ${i} > ${i}.sorted
done
```

### 13.3 Link Gene Annotations

```bash
ln -s ~/sciebo/Postdoc/sexAsex/genomes/geneAnno/*.gff geneAnnotations
```

### 13.4 Intersect SVs with Gene Features

Classify each SV as intergenic, genic, exonic, or intronic based on overlap with gene annotations.

```bash
# Intergenic SVs (no overlap with any gene feature)
for i in *.sorted; do
  b=$(echo $i | cut -f1 -d ".")
  bedtools intersect -a ${i} \
    -b <(cat geneAnnotations/${b}.gff3) -wao \
    | sort -k1,1V -k2,2n | sort \
    | awk '{if ($14 == 0) print $0"\tintergenic"}' \
    | cut -f4,15 | sort | uniq \
    > ${i}.intergenic
done

# Genic SVs (overlap with gene features)
for i in *.sorted; do
  b=$(echo $i | cut -f1 -d ".")
  bedtools intersect -a ${i} \
    -b <(cat geneAnnotations/${b}.gff3) -wao \
    | sort -k1,1V -k2,2n | sort \
    | awk '{if ($14 > 0 && $7 == "gene") print $0}' \
    | cut -f4,13 \
    | awk '{ split($2, a, ";"); sub(/^ID=/, "", a[1]); sub(/\.t.*/, "", a[1]);
             print $1"\t"a[1]"\tgenic" }' \
    | sort | uniq \
    > ${i}.genes
done

# Exonic SVs
for i in *.sorted; do
  b=$(echo $i | cut -f1 -d ".")
  bedtools intersect -a ${i} \
    -b <(cat geneAnnotations/${b}.gff3) -wao \
    | sort -k1,1V -k2,2n | sort \
    | awk '{if ($14 > 0 && $7 == "exon") print $0}' \
    | cut -f4,13 \
    | awk '{ split($2, a, ";"); sub(/^ID=/, "", a[1]); sub(/\.t.*/, "", a[1]);
             print $1"\t"a[1]"\texonic" }' \
    | sort | uniq \
    > ${i}.exons
done

# Intronic SVs
for i in *.sorted; do
  b=$(echo $i | cut -f1 -d ".")
  bedtools intersect -a ${i} \
    -b <(cat geneAnnotations/${b}.gff3) -wao \
    | sort -k1,1V -k2,2n | sort \
    | awk '{if ($14 > 0 && $7 == "intron") print $0}' \
    | cut -f4,13 \
    | awk '{ split($2, a, ";"); sub(/^ID=/, "", a[1]); sub(/\.t.*/, "", a[1]);
             print $1"\t"a[1]"\tintronic" }' \
    | sort | uniq \
    > ${i}.introns
done
```

> Continue with GO term enrichment analysis using **topGO** in R, leveraging the InterProScan GO annotations and the gene lists generated above.

---

## Software Versions

| Tool | Version / Source |
|------|-----------------|
| samtools | — |
| lima | — |
| pbmarkdup | — |
| meryl | — |
| winnowmap | — |
| NanoPlot | — |
| Picard | — |
| GATK | 4.6.0.0 |
| JDK | 17.0.12 |
| WhatsHap | — |
| bcftools | — |
| VCFtools | — |
| PLINK / PLINK2 | — |
| IQ-TREE2 | — |
| msa_split (PHAST) | — |
| PopLDdecay | 3.43 |
| GraffiTE | — |
| Nextflow | — |
| Sniffles2 | (via GraffiTE) |
| seqtk | — |
| BEDTools | — |
| InterProScan | 5.69-101.0 |
| genomics_general | — |
| GNU Parallel | — |