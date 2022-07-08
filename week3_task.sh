
#!bin/bash

echo -e "\n Downloading data... \n"

mkdir -p raw_data 
cd raw_data

wget https://zenodo.org/record/2582555/files/SLGFSK-N_231335_r1_chr5_12_17.fastq.gz
wget https://zenodo.org/record/2582555/files/SLGFSK-N_231335_r2_chr5_12_17.fastq.gz
wget https://zenodo.org/record/2582555/files/SLGFSK-T_231336_r1_chr5_12_17.fastq.gz
wget https://zenodo.org/record/2582555/files/SLGFSK-T_231336_r2_chr5_12_17.fastq.gz

echo -e "\n Downloading reference sequence... \n"

wget https://zenodo.org/record/2582555/files/hg19.chr5_12_17.fa.gz

#To unzip reference
gunzip hg19.chr5_12_17.fa.gz

echo -e "\n Data Preprocessing... \n"

#Make directory for fastqc reports named Fastqc_Reports
mkdir -p Fastqc_Reports

#Checking the quality of the reads
fastqc *.fastq.gz -o Fastqc_Reports

#Multiqc reports
multiqc Fastqc_Reports -o Fastqc_Reports

mkdir -p trimmed_reads

#Use of trimmomatic to trim adapters and poor reads
# download Trimmomatic frm here http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/Trimmomatic-0.39.zip, unzip and copy the specific adapter into your directory
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

multiqc  trimmed_reads/Fastqc_results  -o trimmed_reads/Fastqc_results

#Index reference file	
bwa index hg19.chr5_12_17.fa 

#Perform alignment
bwa mem -R '@RG\tID:231335\tSM:Normal' hg19.chr5_12_17.fa trimmed_reads/SLGFSK-N_231335_r1_paired.fq.gz \
      trimmed_reads/SLGFSK-N_231335_r2_paired.fq.gz > Mapping/SLGFSK-N_231335.sam

bwa mem -R '@RG\tID:231336\tSM:Tumor' hg19.chr5_12_17.fa trimmed_reads/SLGFSK-T_231336_r1_paired.fq.gz \
       trimmed_reads/SLGFSK-T_231336_r2_paired.fq.gz > Mapping/SLGFSK-T_231336.sam

#Convert SAM to BAM and sort it 
for sample in `cat list.txt`
do
        samtools view -@ 20 -S -b Mapping/${sample}.sam | samtools sort -@ 32 > Mapping/${sample}.sorted.bam
 
        Index BAM file
        samtools index Mapping/${sample}.sorted.bam
done

#Filter BAM files
for sample in `cat list.txt`
do
        samtools view -q 1 -f 0x2 -F 0x8 -b Mapping/${sample}.sorted.bam > Mapping/${sample}.filtered1.bam
done

#To remove duplicates, the markdup command was used
for sample in `cat list.txt`
do
	samtools sort -n -o Mapping/${sample}.namesort.bam Mapping/${sample}.filtered1.bam
        samtools fixmate -m Mapping/${sample}.namesort.bam Mapping/${sample}.fixmate.bam
        samtools sort -@ 4 -o Mapping/${sample}.positionsort.bam Mapping/${sample}.fixmate.bam
        samtools markdup -@4 -r Mapping/${sample}.positionsort.bam Mapping/${sample}.markdup.bam
done

#To left align bam
for sample in `cat list.txt`
do      
        cat Mapping/${sample}.markdup.bam  | bamleftalign -f hg19.chr5_12_17.fa -m 5 -c > Mapping/${sample}.leftAlign.bam
done

#Recalibrate read mapping qualities
for sample in `cat list.txt`
do
        samtools calmd -@ 4 -b Mapping/${sample}.leftAlign.bam hg19.chr5_12_17.fa > Mapping/${sample}.recalibrate.bam
done

#Refilter read mapping qualities
for sample in `cat list.txt`
do
        bamtools filter -in Mapping/${sample}.recalibrate.bam -mapQuality "<=254" > Mapping/${sample}.refilter.bam
done

#Variant calling and classification
#Installation
wget https://sourceforge.net/projects/varscan/files/VarScan.v2.3.9.jar

#Convert data to pileup
mkdir Variants

for sample in `cat list.txt`
do
        samtools mpileup -f hg19.chr5_12_17.fa Mapping/${sample}.refilter.bam --min-MQ 1 --min-BQ 28 \
                > Variants/${sample}.pileup
done

#Variant calling
varscan somatic Variants/SLGFSK-N_231335.pileup \
        Variants/SLGFSK-T_231336.pileup Variants/SLGFSK \
        --normal-purity 1 --tumor-purity 0.5 --output-vcf 1 

#zip, index and merge vcf

#Install bcftools; 
conda install -c bioconda bcftools

#zip vcf files, index and merge
bgzip Variants/SLGFSK.snp.vcf > Variants/SLGFSK.snp.vcf.gz
bgzip Variants/SLGFSK.indel.vcf > Variants/SLGFSK.indel.vcf.gz
tabix Variants/SLGFSK.snp.vcf.gz
tabix Variants/SLGFSK.indel.vcf.gz
bcftools merge SLGFSK.snp.vcf.gz SLGFSK.indel.vcf.gz > SLGFSK.vcf

#Functional annotation using snpEff
#download snpEff
wget https://snpeff.blob.core.windows.net/versions/snpEff_latest_core.zip

#unzip file
unzip snpEff_latest_core.zip

#download snpEff database
snpEff download hg19

#to annotate using snpEff
snpEff hg19 Variants/merged.vcf > Varriants/merged.ann.vcf

