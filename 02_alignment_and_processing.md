# Read Alignment
Once trimming was complete for all 32 samples, reads were aligened with BWA and filtered for a minimum mapping quality of 10 with SAMtools.  

This alignment script was modified from this [shell script](https://github.com/clairemerot/wgs_sample_preparation/blob/master/01_scripts/02_bwa_mem_align_reads_PE_with_sample_file.sh) cited in Merot et al.  
```
ref=/media/jana/BigData/BIOL337/reference/GCF_018398675.1_ASM1839867v1_genomic.fna
fqs=/media/jana/BigData/BIOL337/C_clupeaformis_SR/trim_reads/
bam=/media/jana/BigData/BIOL337/C_clupeaformis_SR/alignments/bam/

for samp in ${fqs}*_trimmed_1.fastq.gz
    do
    base=$(basename $samp _trimmed_1.fastq.gz)
    infoline=$(zcat ${samp} | head -n 1)
    instrument=`echo ${infoline} | cut -d ' ' -f1`
    instrumentrun=`echo $infoline | cut -d ' ' -f2`

    #Now to incorporate this information into the alignment
    rgid="ID:${instrument}_${instrumentrun}"
    rgpl="PL:Illumina"
    rgsm="SM:${base}"
    ID="@RG\t${rgid}\t${rgpl}\t${rgsm}"

    echo "Aligning reads for ${base}..."
    bwa mem -t 16 -R ${ID} ${ref} ${samp} ${fqs}${base}_trimmed_2.fastq.gz | samtools view -Sb -q10 - > ${bam}${base}.bam
    
    echo "Sorting ${bam}${base}.bam..."
    samtools sort --threads 16 ${bam}${base}.bam > ${bam}${base}.sorted.bam
    
    echo "Indexing ${bam}${base}.sorted.bam..."
    samtools index ${bam}${base}.sorted.bam
done
```
There were a small handfull of samples that had issues parsing the `$infoline` above (CD22, CD28, CN10, CN11, and CD32). These samples had the forward (.1) or reverse (.2) read notation annotated to the end of the infoline, which was problematic for relating forward and reverse reads to one another. This was fixed as per below.  
```
for redo in ${fqs}redo/*_trimmed_1.fastq.gz
    do
    base=$(basename $fqs _trimmed_1.fastq.gz)
    echo "Fixing reads for ${base}..."
    zcat ${redo} | sed -e 's/.1 / /g' | gzip > ${fqs}redo/${base}_fix_1.fastq.gz
    zcat ${fqs}redo/${base}_trimmed_2.fastq.gz | sed -e 's/.2 / /g' | gzip > ${fqs}redo/${base}_fix_2.fastq.gz
done
```
Read alignments were then performed as per above.  

Once alignments were complete, alignment metrics were estimated using picard v3.0.0. The below script was modified from this [shell script](https://github.com/clairemerot/wgs_sample_preparation/blob/master/01_scripts/03_collect_metrics.sh).  
```
stats=/media/jana/BigData/BIOL337/C_clupeaformis_SR/alignments/bam_stats/

for file in ${bam}*.sorted.bam
    do
    base=$(basename ${file} .sorted.bam)
    echo "Computing alignment metrics for ${base}..."
    picard CollectAlignmentSummaryMetrics \
        -R ${ref} \
        -I ${file} \
        -O ${stats}${base}_alignment_metrics.txt

    echo "Computing insert size metrics for ${base}..."
    picard CollectInsertSizeMetrics \
        -I ${file} \
        -O ${stats}${base}_insert_size_metrics.txt \
        -H ${stats}${base}_insert_size_histogram.pdf
done
```
Finally, duplicates were removed using picard. This [shell script](https://github.com/clairemerot/wgs_sample_preparation/blob/master/01_scripts/04_remove_duplicates.sh) was modified as per below.  
```
nodup=/media/jana/BigData/BIOL337/C_clupeaformis_SR/alignments/nodup_bam/
stats=/media/jana/BigData/BIOL337/C_clupeaformis_SR/alignments/nodup_stats/

for file in ${bam}*.sorted.bam
    do
    base=$(basename ${file} .sorted.bam)
    echo "Removing duplicates for ${base}..."
    picard MarkDuplicates \
        -I ${file} \
        -O ${nodup}${base}_nodup.bam \
        -M ${stats}${base}_metrics.txt \
        -VALIDATION_STRINGENCY SILENT \
        -REMOVE_DUPLICATES true
done
```
Once duplicates were removed, bams were filtered to exclude unplaced scaffolds to expedite later steps.
```
dir=/media/jana/BigData/BIOL337/C_clupeaformis_SR/

for bam in ${nodup}*_nodup.bam
    do
    base=$(basename ${bam} _nodup.bam)
    echo "Subsetting file for ${base}..."
    samtools view --threads 16 \
        -L ${dir}chromosome_scaffolds.bed \
        -b ${bam} > ${nodup}${base}_subset.bam

    echo "Computing alignment metrics for ${base}..."
    picard CollectAlignmentSummaryMetrics \
        -R ${ref} \
        -I ${nodup}${base}_subset.bam \
        -O ${stats}${base}_alignment_metrics.txt

    echo "Computing insert size metrics for ${base}..."
    picard CollectInsertSizeMetrics \
        -I ${nodup}${base}_nodup.bam \
        -O ${stats}${base}_insert_size_metrics.txt \
        -H ${stats}${base}_insert_size_histogram.pdf
    
    echo "Indexing ${base}..."
    picard BuildBamIndex \
        -I ${nodup}${base}_subset.bam \
        -O ${nodup}${base}_subset.bam.bai
done
```
Upon completion of the subsetting of bams with samtools, files were indexed (as above) and a dictionary reference was constructed using picard. Scripts were modified from [Merot et al](https://github.com/clairemerot/wgs_sample_preparation/tree/master/01_scripts).  
```
picard CreateSequenceDictionary -R ${ref} -O ${ref}.dict
samtools faidx ${ref}
```
Realignment around INDELs was then performed using GATK.  
```
for bam in ${nodup}*_subset.bam
    do
    base=$(basename ${bam} _subset.bam)
    echo "Realigning INDELs for ${base}..."
    GenomeAnalysisTK -T RealignerTargetCreator \
        -R ${ref} \
        -I ${bam} \
        -o ${bam}.intervals
    GenomeAnalysisTK -T IndelRealigner \
        -R ${ref} \
        -I ${bam} \
        -targetIntervals ${bam}.intervals \
        --consensusDeterminationModel USE_READS \
        -o ${dir}alignments/realigned/${base}_realigned.bam
done
```
Finally, the [09_clip_overlap.sh](https://github.com/clairemerot/wgs_sample_preparation/blob/master/01_scripts/09_clip_overlap.sh) script provided by Merot et al. was modified to soft clip overlapping read ends using clipOverlap in BamUtil.  
```
for bam in ${dir}realigned/*_realigned.bam
    do
    base=$(basename ${bam} _realigned.bam)
    echo "Clipping overlapping reads for ${base}..."
    bam clipOverlap \
        --in ${bam} \
        --out ${dir}realigned/${base}.temp.bam \
        --unmapped --storeOrig OC --stats
    samtools view --threads 16 -hb -F 4 \
        ${dir}realigned/${base}.temp.bam > ${dir}realigned/${base}_noOverlap.bam
    samtools index ${dir}realigned/${base}_noOverlap.bam
    rm ${dir}realigned/${base}.temp.bam
done
```