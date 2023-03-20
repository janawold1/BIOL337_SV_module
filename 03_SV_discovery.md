# Short-read SV discovery
The 'realigned' bam files were used for SV discovery. Merot et al. called SVs in population batches.  

## Delly
Here, we modified used the script for calling SVs using Delly outlined [here](https://github.com/clairemerot/SR_SV/blob/main/01_scripts/06_delly_CD.sh). The version of Delly used here is v1.1.6.  

SVs were called by grouping populations. It is important to note that ID9 belongs in population CN and that IN9 and CN7 belong in population ID.  
```
dir=/media/jana/BigData/BIOL337/C_clupeaformis_SR/

for pop in CD CN ID IN
    do
    echo "Calling SVs for population ${pop}..."
    if [ "$pop" = "CN" ]
        then
        delly call -t ALL -o ${dir}delly/calls/${pop}.bcf \
            -g ${ref} \
            -q 20 -s 15 \
            ${dir}alignments/realigned/CN10_realigned.bam \
            ${dir}alignments/realigned/CN11_realigned.bam \
            ${dir}alignments/realigned/CN12_realigned.bam \
            ${dir}alignments/realigned/CN14_realigned.bam \
            ${dir}alignments/realigned/CN15_realigned.bam \
            ${dir}alignments/realigned/CN5_realigned.bam \
            ${dir}alignments/realigned/CN6_realigned.bam \
            ${dir}alignments/realigned/ID9_realigned.bam
    elif [ "$pop" = "ID" ]
        then
        delly call -t ALL -o ${dir}delly/calls/${pop}.bcf \
            -g ${ref} \
            -q 20 -s 15 \
            ${dir}alignments/realigned/CN7_realigned.bam \
            ${dir}alignments/realigned/ID13_realigned.bam \
            ${dir}alignments/realigned/ID14_realigned.bam \
            ${dir}alignments/realigned/ID1_realigned.bam \
            ${dir}alignments/realigned/ID2_realigned.bam \
            ${dir}alignments/realigned/ID3_realigned.bam \
            ${dir}alignments/realigned/ID4_realigned.bam \
            ${dir}alignments/realigned/ID7_realigned.bam \
            ${dir}alignments/realigned/IN9_realigned.bam
    elif [ "$pop" = "IN" ]
        then
        delly call -t ALL -o ${dir}delly/calls/${pop}.bcf \
            -g ${ref} \
            -q 20 -s 15 \
            ${dir}alignments/realigned/IN10_realigned.bam \
            ${dir}alignments/realigned/IN12_realigned.bam \
            ${dir}alignments/realigned/IN14_realigned.bam \
            ${dir}alignments/realigned/IN5_realigned.bam \
            ${dir}alignments/realigned/IN6_realigned.bam \
            ${dir}alignments/realigned/IN7_realigned.bam \
            ${dir}alignments/realigned/IN8_realigned.bam
    else
        delly call -t ALL -o ${dir}delly/calls/${pop}.bcf \
            -g ${ref} \
            -q 20 -s 15 \
            ${dir}alignments/realigned/${pop}*_realigned.bam
    fi
done
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