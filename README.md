# HackBio-Workshop-Week-3

## The week's task was to reproduce a project on Somatic and Germline variant Identification from Tumor and normal Sample Pairs

### First, the datasets to be used were downloaded from Zenodo

## Command

echo -e "\n Downloading data... \n"

mkdir -p raw_data 
cd raw_data

wget https://zenodo.org/record/2582555/files/SLGFSK-N_231335_r1_chr5_12_17.fastq.gz
wget https://zenodo.org/record/2582555/files/SLGFSK-N_231335_r2_chr5_12_17.fastq.gz
wget https://zenodo.org/record/2582555/files/SLGFSK-T_231336_r1_chr5_12_17.fastq.gz
wget https://zenodo.org/record/2582555/files/SLGFSK-T_231336_r2_chr5_12_17.fastq.gz

## Reference sequence
echo -e "\n Downloading reference sequence... \n"

wget https://zenodo.org/record/2582555/files/hg19.chr5_12_17.fa.gz

## Quality check and adapter trimming
### Fastqc was used to check the quality of the datasets and trimmomatic was used to trim the adapter sequence

mkdir -p Fastqc_Reports

fastqc *.fastq.gz -o Fastqc_Reports

## Trimmomatic was downloaded from here

http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/Trimmomatic-0.39.zip

### Command

mkdir -p trimmed_reads

for sample in `cat list.txt`
do
       trimmomatic PE -threads 8 ${sample}_r1_chr5_12_17.fastq.gz ${sample}_r2_chr5_12_17.fastq.gz \
               trimmed_reads/${sample}_r1_paired.fq.gz trimmed_reads/${sample}_r1_unpaired.fq.gz \
               trimmed_reads/${sample}_r2_paired.fq.gz trimmed_reads/${sample}_r2_unpaired.fq.gz \
               ILLUMINACLIP:TruSeq3-PE.fa:2:30:10:8:keepBothReads \
               LEADING:3 TRAILING:10 MINLEN:25

        fastqc  trimmed_reads/${sample}_r1_paired.fq.gz  trimmed_reads/${sample}_r2_paired.fq.gz \
                 -o trimmed_reads/Fastqc_results
                 
done 

