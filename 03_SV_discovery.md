# Short-read SV discovery
The scripts below outline how we're going to call structural variants for this course. The 'realigned' bam files were used for SV discovery.

For these exercises, we are going to run the SV discovery tools Delly and Smoove. Do not worry about Manta as we'll go over those results together. Another note is that we're going to be running our analysis for chromosome 1 only. This saves time and computational requirements. All scripts have been simplified from Merot et al. 2023. For original source code, please check out [here](https://github.com/clairemerot/SR_SV/).  

## Delly
Here, we modified used the script for calling SVs using Delly as outlined [here](https://github.com/clairemerot/SR_SV/blob/main/01_scripts/06_delly_CD.sh). The version of Delly used here is v1.1.6.  

SVs were called by grouping populations. It is important to note that ID9 belongs in population CN and that IN9 and CN7 belong in population ID.  

### Merot et al. 2023 called SVs in population batches. As we go throught this process, think about what are some of the benefits and caveats of this approach.  

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
We are now ready to filter our SV calls.  

### What is 'filtering' and why might we want to perform this step?  

  
```
bcftools filter \
    -i 'INFO/SVTYPE!="BND" & INFO/SVTYPE!="INS"' \
    -O b -o ${dir}delly/delly_noTRA.bcf \
    ${dir}delly/delly_merged.bcf
```
### How many SVs did this leave?

This custom [Rscript](https://github.com/clairemerot/SR_SV/blob/main/01_scripts/Rscripts/add_explicit_seq.r) was used to convert sequence IDs. The code sourced below can be found [here](https://github.com/clairemerot/SR_SV/blob/main/01_scripts/Rscripts/fix_sniffles.R).  
```
script=/home/jwold/scripts/
```

## Smoove
While Delly was run for grouped populations, Smoove v0.2.8 was run by batching chromosomes (i.e., chr 1-10, 11-20, 21-30, & 31-40). Multi-threading with the `-p` flag caused Smoove to crash, therefore this option was not used.  
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
        ${dir}alignments/realigned/*_realigned.bam
done

```
## Manta
SVs called using Manta were jointly called across all chromosomes.  
```
bgzip chromosome_scaffolds.bed
tabix chromosome_scaffolds.bed.gz

configManta.py --referenceFasta ${ref} \
    --callRegions chromosome_scaffolds.bed.gz \
    --runDir manta/calls/
    --bam alignments/realigned/CD17_realigned.bam \
    --bam alignments/realigned/CD18_realigned.bam \
    --bam alignments/realigned/CD19_realigned.bam \
    --bam alignments/realigned/CD20_realigned.bam \
    --bam alignments/realigned/CD21_realigned.bam \
    --bam alignments/realigned/CD22_realigned.bam \
    --bam alignments/realigned/CD28_realigned.bam \
    --bam alignments/realigned/CD32_realigned.bam \
    --bam alignments/realigned/CN10_realigned.bam \
    --bam alignments/realigned/CN11_realigned.bam \
    --bam alignments/realigned/CN12_realigned.bam \
    --bam alignments/realigned/CN14_realigned.bam \
    --bam alignments/realigned/CN15_realigned.bam \
    --bam alignments/realigned/CN5_realigned.bam \
    --bam alignments/realigned/CN6_realigned.bam \
    --bam alignments/realigned/CN7_realigned.bam \
    --bam alignments/realigned/ID13_realigned.bam \
    --bam alignments/realigned/ID14_realigned.bam \
    --bam alignments/realigned/ID1_realigned.bam \
    --bam alignments/realigned/ID2_realigned.bam \
    --bam alignments/realigned/ID3_realigned.bam \
    --bam alignments/realigned/ID4_realigned.bam \
    --bam alignments/realigned/ID7_realigned.bam \
    --bam alignments/realigned/ID9_realigned.bam \
    --bam alignments/realigned/IN10_realigned.bam \
    --bam alignments/realigned/IN12_realigned.bam \
    --bam alignments/realigned/IN14_realigned.bam \
    --bam alignments/realigned/IN5_realigned.bam \
    --bam alignments/realigned/IN6_realigned.bam \
    --bam alignments/realigned/IN7_realigned.bam \
    --bam alignments/realigned/IN8_realigned.bam \
    --bam alignments/realigned/IN9_realigned.bam
```
