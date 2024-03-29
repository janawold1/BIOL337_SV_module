# Data download and read preparation
## Downloading short-read data
Short-read data for [Mérot et al. 2022](https://onlinelibrary.wiley.com/doi/10.1111/mec.16468) is accessible on the NCBI SRA database under project number PRJNA820751. For consistency and reproducibility, code used here is modified directly from the author's [github repo](https://github.com/clairemerot).  

While genome assemblies for both the '[normal](https://www.ncbi.nlm.nih.gov/data-hub/genome/GCF_018398675.1/)' and '[dwarf](https://www.ncbi.nlm.nih.gov/data-hub/genome/GCF_020615455.1/)' phenotype are available, the genome assembly for the normal phenotype was used for short-read based alignment. Here, the `to_download` file is the SR names from NCBI project PRJNA820751.  
```
work=/media/jana/BigData/BIOL337/C_clupeaformis_SR/

while read -r line
    do
    echo "Downloading ${line}..."
    fasterq-dump --verbose --split-files --progress -e 16 -O ${work}raw_reads/ ${line}
done < to_download
```
Sample IDs were then renamed for simplicity and to be consistent with Merot et al. Fastqc was then run on raw reads for an initial quality check.  
The `rename_list` file was constructed using metadata in the NCBI project PRJNA820751.
```
while read -r line
    do
    sr=$(echo $line | awk '{print $1}')
    id=$(echo $line | awk '{print $2}')
    echo "Renaming ${sr} to ${id}..."
    rename "s/${sr}/${id}/g" ${work}raw_reads/${sr}
    echo "Running FastQC for ${id}..."
    fastqc --threads 16 -o ${work}raw_fastqc ${work}raw_reads/${id}*
done < rename_list

multiqc --title "Raw read QC" ${dir}raw_fastqc/
```

## Read trimming
The for loop below was modified from [01_fastp.sh](https://github.com/clairemerot/wgs_sample_preparation/blob/master/01_scripts/01_fastp.sh). Because it is unclear whether there were any additional options used, we will run under default settings (as it was indicated in this script) and added a hard trim for 20bp and 5bp at the start and end of reads to remove noise.

```
mkdir -p ${work}{trim_reads,trim_reports}

for file in ${work}raw_reads/*_1.fastq
do
    base=$(basename $file _1.fastq)
    echo "Treating: ${base}..."
    fastp -w 6 \
        -i ${work}raw_reads/${base}_1.fastq \
        -I ${work}raw_reads/${base}_2.fastq \
        -f 20 /
        -t 5 /
        -o ${work}trim_reads/${base}_trimmed_1.fastq.gz \
        -O ${work}trim_reads/${base}_trimmed_2.fastq.gz \
        -j ${work}trim_reports/${base}.json \
        -h ${work}trim_reports/${base}.html
done
```

Once read trimming was complete, outputs were collated and examined using MultiQC.  
```
multiqc --title "Trimmed read QC" ${work}trim_reports/
```