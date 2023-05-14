# Short-read SV discovery
The scripts below outline how we're going to call structural variants for this course. The 'realigned' bam files were used for SV discovery.

For these exercises, we are going to run the SV discovery tools Delly and Smoove. Do not worry about Manta as we'll go over those results together. Another note is that we're going to be running our analysis for chromosome 1 only. This saves time and computational requirements. All scripts have been simplified from Merot et al. 2023. For original source code, please check out [here](https://github.com/clairemerot/SR_SV/).  

## Delly
Here, we modified used the script for calling SVs using Delly as outlined [here](https://github.com/clairemerot/SR_SV/blob/main/01_scripts/06_delly_CD.sh). The version of Delly used here is v1.1.6.  

To call SVs with Delly, Merot et al. chose to group populations. What are some potential reasons why this would be done?

As a side note, sample names may not always be straightforward. As an example, sample ID9 belongs in population CN while IN9 and CN7 belong in population ID.  

```
dir=/home/jwold/
ref=/home/jwold/GCF_018398675.1_ASM1839867v1_genomic.fna

for pop in "CD" "CN" "ID" "IN"
    do
    echo "Calling SVs for population ${pop}..."
    if [ "$pop" = "CN" ]
        then
        delly call -t ALL -o ${dir}delly/calls/${pop}.bcf \
            -g ${ref} \
            -q 20 -s 15 \
            ${dir}alignments/CN10_chr1.bam \
            ${dir}alignments/CN11_chr1.bam \
            ${dir}alignments/CN12_chr1.bam \
            ${dir}alignments/CN14_chr1.bam \
            ${dir}alignments/CN15_chr1.bam \
            ${dir}alignments/CN5_chr1.bam \
            ${dir}alignments/CN6_chr1.bam \
            ${dir}alignments/ID9_chr1.bam
    elif [ "$pop" = "ID" ]
        then
        delly call -t ALL -o ${dir}delly/calls/${pop}.bcf \
            -g ${ref} \
            -q 20 -s 15 \
            ${dir}alignments/CN7_chr1.bam \
            ${dir}alignments/ID13_chr1.bam \
            ${dir}alignments/ID14_chr1.bam \
            ${dir}alignments/ID1_chr1.bam \
            ${dir}alignments/ID2_chr1.bam \
            ${dir}alignments/ID3_chr1.bam \
            ${dir}alignments/ID4_chr1.bam \
            ${dir}alignments/ID7_chr1.bam \
            ${dir}alignments/IN9_chr1.bam
    elif [ "$pop" = "IN" ]
        then
        delly call -t ALL -o ${dir}delly/calls/${pop}.bcf \
            -g ${ref} \
            -q 20 -s 15 \
            ${dir}alignments/IN10_chr1.bam \
            ${dir}alignments/IN12_chr1.bam \
            ${dir}alignments/IN14_chr1.bam \
            ${dir}alignments/IN5_chr1.bam \
            ${dir}alignments/IN6_chr1.bam \
            ${dir}alignments/IN7_chr1.bam \
            ${dir}alignments/IN8_chr1.bam
    else
        delly call -t ALL -o ${dir}delly/calls/${pop}.bcf \
            -g ${ref} \
            -q 20 -s 15 \
            ${dir}alignments/${pop}*_chr1.bam
    fi
done
```
We now have four population-level calls, but we don't know how SVs called within each population may relate to other populations. To ensure we have a call set that is comparable across each of our populations, we must merge them.  
```
delly merge -o ${dir}delly_merged.bcf ${dir}calls/*.bcf
```
### Filtering of Delly SVs
We are now ready to filter our SV calls. Thus far, filtering for SVs called by short-read data can be regarded as 'context dependent'. Decisions are made depending on the tools used, characteristics of the sequence data (e.g., read legth, depth and insert length) and the contiguity of your reference assembly. Below we outline how Merot et al. performed filtering on the raw population call sets with Delly.  

To begin with, all breakend and translocation calls were removed with BCFtools.
```
bcftools view -i 'SVTYPE!="BND"' \
    -O v -o ${dir}delly/delly_noBND.vcf \
    ${dir}delly/delly_merged.bcf
```
But we're not done on our filtering journey just yet! Delly does not natively calculate SV length, so we have to infer it and add it to our new VCF manually. A slight variation on the method implemented by Merot et al. is presented below.  
```
# Estimating length and preparing to add it to the VCF
bcftools query -f '%CHROM\t%POS\t%INFO/END\n' ${dir}delly/delly_noBND.bcf | awk '{print $1"\t"$2"\t"$3-$2}' > ${dir}delly/size.annot
bgzip ${dir}delly/size.annot
tabix -s1 -b2 -e2 ${dir}delly/size.annot.gz

# Creating a new line for VCF header
echo -e '##INFO=<ID=SVLEN,Number=.,Type=Integer,Description="Difference in length between REF and ALT alleles">' > ${dir}delly/annot.hdr

# Annotating the VCF with the new field for SV length
bcftools annotate -a ${dir}delly/size.annot.gz -h ${dir}delly/annot.hdr -c CHROM,POS,INFO/SVLEN ${dir}delly/delly_noBND.vcf > ${dir}delly/delly_annot.vcf
```
With this step complete, we can now implement length thresholds for our SVs. Merot et al. chose to retain all SVs > 50bp and < 100kb in length.  
```
bcftools view -i 'SVLEN>=50 & SVLEN<=100000' -O v -o ${dir}delly/delly_size_filtered.vcf ${dir}delly/delly_annot.vcf 
```
Finally, we want to identify any SVs that were called in poorly assembled regions of the genome, or those with poorly resolved sequences. To do so, we'll leverage [this](https://github.com/clairemerot/SR_SV/blob/main/01_scripts/Rscripts/fix_sniffles_delly.R) custom Rscript from Merot et al. to infer SV sequences using the VCF and reference genome. This script works best on R v4.3 and requires a few dependencies to be installed before running.  

First start your R terminal, and install the required packages.  

```
R

if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("Rsamtools")
BiocManager::install("GenomicRanges")
```
Now we're ready to set file locations, source script and perform conversion.  
```
INPUT <- ("/home/jwold/delly/delly_size_filtered.vcf")
OUTPUT <- ("/home/jwold/delly/delly_size_filtered_withSEQ.vcf")
GENOME <- ("/home/jwold/reference/GCF_018398675.1_ASM1839867v1_genomic.fna")

source("/home/jwold/scripts/fix_sniffles_delly.R")
```
Now we're ready to run the function.  
```
fix_sniffles(input_vcf=INPUT, output_vcf=OUTPUT, refgenome=GENOME)
```
Finally, we can remove any SV calls that occurred in poorly assembled regions of the genome, or those with poorly resolved sequences.  

In this example, we are going to create 'blacklists'. This is essentially a list of sites that we will want to remove. We first start with summarising SV location, start and end points, type, length, the reference allele and alternate allele with bcftools.  
```
bcftools query -f '%CHROM %POS %INFO/END %INFO/SVTYPE %INFO/SVLEN %REF %ALT\n' ${dir}delly/delly_size_filtered_withSEQ.vcf > ${dir}delly/SV_data_with_seq.txt

```
In our example, we're going to remove any SV where there were > 10 uncalled basepairs.
```
grep -P "N{10,}" ${dir}delly/SV_data_with_seq.txt | awk '{print $1 "\t" $2 "\t" $6 "\t" $7}' > ${dir}delly/N10_blacklist.bed
echo "SVs excluded because of >10N" 
wc -l ${dir}delly/N10_blacklist.bed
```
As well as any SV where the reference sequence could not be resolved.
```
cat  ${dir}delly/SV_data_with_seq.txt | awk '{if ($6 == "N") print $1 "\t" $2 "\t" $6 "\t" $7;}' > ${dir}delly/N_blacklist.bed
echo "SVs excluded because absence of sequence ref" 
wc -l ${dir}delly/N_blacklist.bed
```
And any where the SV sequence could not be resolved.
```
cat  ${dir}delly/SV_data_with_seq.txt | awk '{print $1 "\t" $2 "\t" $6 "\t" $7}' | grep -P "<" > ${dir}delly/N_blacklist_bis.bed
echo "SVs excluded because absence of sequence alt" 
wc -l ${dir}delly/N_blacklist_bis.bed
```
We can now join all blacklists into one file for removal... But did you notice anything interesting about our blacklists?
```
cat ${dir}delly/N_blacklist.bed ${dir}delly/N_blacklist_bis.bed ${dir}delly/N10_blacklist.bed | sort -k1,1 -k2,2n > ${dir}delly/blacklist.bed
bgzip -c ${dir}delly/blacklist.bed > ${dir}delly/blacklist.bed.gz
tabix -s1 -b2 -e2 ${dir}delly/blacklist.bed.gz

# Finally removal of bad sites
bcftools view -T ^${dir}delly/blacklist.bed.gz ${dir}delly/delly_size_filtered_withSEQ.vcf > ${dir}delly/delly_Nfiltered.vcf
echo "SVs after filtration for N seq" 
grep -v ^\#\# ${dir}delly/delly_Nfiltered.vcf | wc -l
```
 * What is 'filtering' and why might we want to perform this step?  
 * How many SVs did this leave? Was there an SV type more heavily impacted by our quality filters?

## Smoove
While Delly was run for grouped populations, Smoove v0.2.8 was run by batching chromosomes (i.e., chr 1-10, 11-20, 21-30, & 31-40). Multi-threading with the `-p` flag can cause Smoove to crash, so I don't generally suggest running it.

Below is an example of how I would run SMOOVE in a similar way to Merot et al. for the whole genome. Since we have bam files representing chromosome 1 and we're only interested in chromosome 1, how might we augment this script?  
```
mkdir ${dir}smoove/calls

for chr in chr1-10 chr11-20 chr21-30 chr31-40
    do
    mkdir ${dir}smoove/calls/${chr}
    echo "Running smoove for ${chr}..."
    smoove call -x --name all_${chr} \
        --outdir ${dir}smoove/calls/${chr} \
        --exclude ${dir}retain_${chr}.bed \
        --fasta $ref --duphold --genotype \
        ${dir}alignments/*_chr1.bam
done
```
Unlike Delly, SMOOVE annotates the file with SV lengths. So we can immediately remove any breakend calls and filter for SVs that are >=50bp and <=100kb in length.
```
bcftools view \
    -i 'SVTYPE!="BND" & ((SVLEN>=50 && SVLEN<=100000) | (SVLEN<=-50 && SVLEN>=-100000))' \
    -O v -o ${dir}smoove/smoove_size_filtered.vcf \
    ${dir}smoove/smoove_merged.vcf
```
And, just like with the Delly SV call set, we need to remove sites called in poorly resolved regions of the genome. Again, we add sequence info inferred from the VCF and the reference genome.
```
R

INPUT <- ("/home/jwold/delly/delly_size_filtered.vcf")
OUTPUT <- ("/home/jwold/delly/delly_size_filtered_withSEQ.vcf")
GENOME <- ("/home/jwold/reference/GCF_018398675.1_ASM1839867v1_genomic.fna")

source("/home/jwold/scripts/fix_sniffles_delly.R")

fix_sniffles(input_vcf=INPUT, output_vcf=OUTPUT, refgenome=GENOME)
```
Summarise SVs for filtering.
```
bcftools query -f '%CHROM %POS %INFO/END %INFO/SVTYPE %INFO/SVLEN %REF %ALT\n' ${dir}smoove/smoove_withSEQ.vcf > ${dir}smoove/SV_data_with_seq.txt
```

And as before, we filtered out SVs in poorly assembled regions of the genome or those that could not be resolved using blacklists.
```
# SVs with >10 unknown basepairs
grep -P "N{10,}" ${dir}smoove/SV_data_with_seq.txt | awk '{print $1 "\t" $2 "\t" $6 "\t" $7}' > ${dir}smoove/N10_blacklist.bed
echo "SVs excluded because of >10N" 
wc -l ${dir}smoove/N10_blacklist.bed

# SVs where the reference sequences is unresolved
cat  ${dir}smoove/SV_data_with_seq.txt | awk '{if ($6 == "N") print $1 "\t" $2 "\t" $6 "\t" $7;}' > ${dir}smoove/N_blacklist.bed
echo "SVs excluded because absence of sequence ref" 
wc -l ${dir}smoove/N_blacklist.bed

# SVs where the variant sequence is unresolved
cat  ${dir}smoove/SV_data_with_seq.txt | awk '{if ($7 == "N") print $1 "\t" $2 "\t" $6 "\t" $7;}' > ${dir}smoove/N_blacklist_bis.bed
echo "SVs excluded because absence of sequence alt" 
wc -l ${dir}smoove/N_blacklist_bis.bed

# Now merging all blacklists into one and preparing for filtering
cat ${dir}smoove/N_blacklist.bed ${dir}smoove/N_blacklist_bis.bed ${dir}smoove/N10_blacklist.bed | sort -k1,1 -k2,2n > ${dir}smoove/blacklist.bed
bgzip -c ${dir}smoove/blacklist.bed > ${dir}smoove/blacklist.bed.gz
tabix -s1 -b2 -e2 ${dir}smoove/blacklist.bed.gz
```
Finally, filtering poor quality SVs with our final blacklist and BCFtools.
```
bcftools view -T ^${dir}smoove/blacklist.bed.gz ${dir}smoove/smoove_withSEQ.vcf > ${dir}smoove/smoove_withSEQ_Nfiltered.vcf
echo "SVs after filtration for N seq" 
grep -v ^\#\# ${dir}smoove/smoove_withSEQ_Nfiltered.vcf | wc -l
```
And with that, we have completed filtering for our SMOOVE data set!
## Manta
Due to time constraints, we have performed all SV discovery and filtering for Manta prior to this lecture. However, we have provided the steps used to generate the data we explore below. To begin, SVs called using Manta were jointly called across all chromosomes and chromosome 1 `NC_055067.1` was extracted to ensure they were comparable with the Delly and Smoove datasets above.  
```
bgzip chromosome_scaffolds.bed
tabix chromosome_scaffolds.bed.gz

configManta.py --referenceFasta ${ref} \
    --callRegions chromosome_scaffolds.bed.gz \
    --runDir manta/calls/
    --bam alignments/CD17_realigned.bam \
    --bam alignments/CD18_realigned.bam \
    --bam alignments/CD19_realigned.bam \
    --bam alignments/CD20_realigned.bam \
    --bam alignments/CD21_realigned.bam \
    --bam alignments/CD22_realigned.bam \
    --bam alignments/CD28_realigned.bam \
    --bam alignments/CD32_realigned.bam \
    --bam alignments/CN10_realigned.bam \
    --bam alignments/CN11_realigned.bam \
    --bam alignments/CN12_realigned.bam \
    --bam alignments/CN14_realigned.bam \
    --bam alignments/CN15_realigned.bam \
    --bam alignments/CN5_realigned.bam \
    --bam alignments/CN6_realigned.bam \
    --bam alignments/CN7_realigned.bam \
    --bam alignments/ID13_realigned.bam \
    --bam alignments/ID14_realigned.bam \
    --bam alignments/ID1_realigned.bam \
    --bam alignments/ID2_realigned.bam \
    --bam alignments/ID3_realigned.bam \
    --bam alignments/ID4_realigned.bam \
    --bam alignments/ID7_realigned.bam \
    --bam alignments/ID9_realigned.bam \
    --bam alignments/IN10_realigned.bam \
    --bam alignments/IN12_realigned.bam \
    --bam alignments/IN14_realigned.bam \
    --bam alignments/IN5_realigned.bam \
    --bam alignments/IN6_realigned.bam \
    --bam alignments/IN7_realigned.bam \
    --bam alignments/IN8_realigned.bam \
    --bam alignments/IN9_realigned.bam

bcftools view -t NC_055067.1 -O v -o ${dir}manta/manta_raw_chr1.vcf ${dir}manta/results/variants/diploidSV.vcf.gz
```
This call set was then filtered to exlude all translocations and breakend calls, SVs smaller than 50bp and larger than 100kb in length.
```
bcftools view -i '(SVTYPE!="TRA" & SVTYPE!="BND") & (SVLEN>=50 | SVLEN<=-50) & (SVLEN<=100000 | SVLEN>=-100000)'\
    -O v -o ${dir}manta/manta_noBND.vcf\
    ${dir}manta/manta_raw_chr1.vcf
```

With the help of a custom [Rscript](https://github.com/clairemerot/SR_SV/blob/main/01_scripts/Rscripts/fix_sniffles_manta.R) we identified any SVs called in poorly assembled regions of the genome (i.e., SVs called in regions of the genome with NNNN sequences).  However, this script uses packages that are most compatible with the most recent version of R (v4.3) and are a bit hungry on memory. 
```
# Start the R terminal
R

# Now install dependent packages
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("Rsamtools")
BiocManager::install("GenomicRanges")

# Ready to set file locations, source script and perform conversion
INPUT <- ("/home/jwold/manta/manta_noBND.vcf")
OUTPUT <- ("/home/jwold/manta/manta_withSEQ.vcf")
GENOME <- ("/home/jwold/reference/GCF_018398675.1_ASM1839867v1_genomic.fna")

source("/home/jwold/scripts/fix_sniffles_manta.R")

# Now we're ready to run the function
fix_sniffles(input_vcf=INPUT, output_vcf=OUTPUT, refgenome=GENOME)
```
With the reference used to fill in the `ALT` field of our VCF, we are now able to identify SVs called in poorly assembled regions with some bash wizardry. First, we identify sequences with more than 10 unresolved basepairs.
```
bcftools query -f '%CHROM %POS %INFO/END %INFO/SVTYPE %INFO/SVLEN %REF %ALT\n' ${dir}manta/manta_withSEQ.vcf > ${dir}manta/SV_data_with_seq.txt
grep -P "N{10,}" ${dir}manta/SV_data_with_seq.txt | awk '{print $1 "\t" $2 "\t" $6 "\t" $7}' > ${dir}manta/N10_blacklist.bed
echo "SVs excluded because of >10N" 
wc -l ${dir}manta/N10_blacklist.bed
```
We also exclude any SVs where the reference sequence could not be resolved.  
```
cat  ${dir}manta/SV_data_with_seq.txt | awk '{if ($6 == "N") print $1 "\t" $2 "\t" $6 "\t" $7;}' > ${dir}manta/N_blacklist.bed
echo "SVs excluded because absence of sequence ref" 
wc -l ${dir}manta/N_blacklist.bed
```
And finally we identified how many SVs had unresolved sequences.
```
cat  ${dir}manta/SV_data_with_seq.txt | awk '{print $1 "\t" $2 "\t" $6 "\t" $7}' | grep -P "<" > ${dir}manta/N_blacklist_bis.bed
echo "SVs excluded because absence of sequence alt" 
wc -l ${dir}manta/N_blacklist_bis.bed
```
The final blacklist was merged and prepared for filtering. 
```
cat ${dir}manta/N_blacklist.bed ${dir}manta/N_blacklist_bis.bed ${dir}manta/N10_blacklist.bed | sort -k1,1 -k2,2n > ${dir}manta/blacklist.bed
head ${dir}manta/blacklist.bed
bgzip -c ${dir}manta/blacklist.bed > ${dir}manta/blacklist.bed.gz
tabix -s1 -b2 -e2 ${dir}manta/blacklist.bed.gz
```
And these variants were finally removed.
```
bcftools view -T ^${dir}manta/blacklist.bed.gz ${dir}manta/manta_withSEQ.vcf > ${dir}manta/manta_withSEQ_Nfiltered.vcf
echo "SVs after filtration for N seq" 
grep -v ^\#\# ${dir}manta/manta_withSEQ_Nfiltered.vcf | wc -l 
```
## Summaries
Finally, we're ready to create a quick summary file for visualising the characteristics of our SVs in R. For this, we're going to be using BCFtools. Remember that the unfiltered SVs for Delly don't have length profiles, so there's a few extra steps.  
```
mkdir ~/SV_summaries

# Summary for Delly. 
printf "CHROM\tPOS\tSVLEN\tSVTYPE\tData_set\n" > ${dir}delly_summary.tsv
bcftools query -f '%CHROM\t%POS\t%INFO/END\t%SVTYPE\tDelly_unfiltered\n' ${dir}delly/delly_merged.bcf | awk '{print $1"\t"$2"\t"$3-$2"\t"$4"\t"$5}' >> ${dir}delly_summary.tsv
bcftools query -f '%CHROM\t%POS\t%SVLEN\t%SVTYPE\tDelly_filtered\n' ${dir}delly/delly_Nfiltered.vcf >> ${dir}delly_summary.tsv
```
Now we can move onto the summary for Smoove and Manta respectively. 
```
printf "CHROM\tPOS\tSVLEN\tSVTYPE\tData_set\n" > ${dir}smoove_summary.tsv
bcftools query -f '%CHROM\t%POS\t%SVLEN\t%SVTYPE\tSmoove_unfiltered\n' ${dir}smoove/smoove_merged.vcf >> ${dir}smoove_summary.tsv
bcftools query -f '%CHROM\t%POS\t%SVLEN\t%SVTYPE\tSmoove_filtered\n' ${dir}smoove/smoove_withSEQ_Nfiltered.vcf >> ${dir}smoove_summary.tsv
sed -i 's/-//g' ${dir}smoove_summary.tsv

printf "CHROM\tPOS\tSVLEN\tSVTYPE\tData_set\n" > ${dir}manta_summary.tsv
bcftools query -f '%CHROM\t%POS\t%SVLEN\t%SVTYPE\tManta_unfiltered\n' ${dir}manta/manta_raw_chr1.vcf >> ${dir}manta_summary.tsv
bcftools query -f '%CHROM\t%POS\t%SVLEN\t%SVTYPE\tManta_filtered\n' ${dir}manta/manta_withSEQ_Nfiltered.vcf >> ${dir}manta_summary.tsv
sed -i 's/-//g' ${dir}manta_summary.tsv
```