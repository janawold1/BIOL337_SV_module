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
    rgpl="PL:${platform}"
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
Once alignments were complete, alignment metrics were estimated using picard v3.0.0. The below script was modified from this [shell script](https://github.com/clairemerot/wgs_sample_preparation/blob/master/01_scripts/03_collect_metrics.sh).  
```
stats=/media/jana/BigData/BIOL337/C_clupeaformis_SR/alignments/bam_stats/

for file in ${bam}*.sorted.bam
    do
    base=$(basename ${file} .sorted.bam)
    echo "Computing alignment metrics for ${base}..."
    picard CollectAlignmentSummaryMetrics \
        R=${ref} \
        I=${file} \
        O=${stats}${base}_alignment_metrics.txt

    echo "Computing insert size metrics for ${base}..."
    picard CollectInsertSizeMetrics \
        I=${file} \
        O=${stats}${base}_insert_size_metrics.txt \
        H=${stats}${base}_insert_size_histogram.pdf
done
```
Finally, duplicates were removed using picard. This [shell script](https://github.com/clairemerot/wgs_sample_preparation/blob/master/01_scripts/04_remove_duplicates.sh) was modified as per below.  
```
nodup=/media/jana/BigData/BIOL337/C_clupeaformis_SR/alignments/nodup_bam/
stats=/media/jana/BigData/BIOL337/C_clupeaformis_SR/alignments/nodup_stats/

for file in ${bam}*.sorted.bam
    do
    base=$(basename ${file} .sorted.bam)
    echo "Computing alignment metrics for ${base}..."
    picard MarkDuplicates \
        I=${file} \
        O=${nodup}${base}_nodup.bam \
        M=${stats}${base}_metrics.txt \
        VALIDATION_STRINGENCY=SILENT \
        REMOVE_DUPLICATES=true 
done
```