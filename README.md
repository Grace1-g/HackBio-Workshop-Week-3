# Somatic and Germline variant Identification from Tumor and normal Sample Pairs (Reproduced)

### The datasets to be used were downloaded from Zenodo

#### Command

       echo -e "\n Downloading data... \n"

       mkdir -p raw_data 
       cd raw_data

       wget https://zenodo.org/record/2582555/files/SLGFSK-N_231335_r1_chr5_12_17.fastq.gz
       wget https://zenodo.org/record/2582555/files/SLGFSK-N_231335_r2_chr5_12_17.fastq.gz
       wget https://zenodo.org/record/2582555/files/SLGFSK-T_231336_r1_chr5_12_17.fastq.gz
       wget https://zenodo.org/record/2582555/files/SLGFSK-T_231336_r2_chr5_12_17.fastq.gz

### Reference sequence
       echo -e "\n Downloading reference sequence... \n"

       wget https://zenodo.org/record/2582555/files/hg19.chr5_12_17.fa.gz
       
       #unzip reference
       gunzip hg19.chr5_12_17.fa.gz

### Pre-processing and Trimming (Quality check and adapter trimming)
#### Fastqc was used to check the quality of the datasets and trimmomatic was used to trim the adapter sequence

       mkdir -p Fastqc_Reports

       fastqc *.fastq.gz -o Fastqc_Reports
       multiqc Fastqc_Reports -o Fastqc_Reports
       
#### Trimmomatic adapter sequence was downloaded from here

	http://www.usadellab.org/cms/?page=trimmomatic

	http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/Trimmomatic-0.39.zip

#### Installation

       conda install -c bioconda trimmomatic --yes

#### Command

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
       multiqc  trimmed_reads/Fastqc_results  -o trimmed_reads/Fastqc_results

### Mapping
The sample sequences were mapped against the reference genome using bwa-mem and the output are in the sam (Single Alignment Map) format.

#### Installation
       conda install -y -c bioconda bwa
       conda install -c bioconda samtools
       conda install -c bioconda bamtools

#### Command
#### First, the reference genome was indexed followed by the alignment of the sample sequences 

       mkdir Mapping
       
       bwa index hg19.chr5_12_17.fa
       
       #Perform alignment
       bwa mem -R '@RG\tID:231335\tSM:Normal' hg19.chr5_12_17.fa trimmed_reads/SLGFSK-N_231335_r1_paired.fq.gz \
             trimmed_reads/SLGFSK-N_231335_r2_paired.fq.gz > Mapping/SLGFSK-N_231335.sam

       bwa mem -R '@RG\tID:231336\tSM:Tumor' hg19.chr5_12_17.fa trimmed_reads/SLGFSK-T_231336_r1_paired.fq.gz \
              trimmed_reads/SLGFSK-T_231336_r2_paired.fq.gz > Mapping/SLGFSK-T_231336.sam


       
### Conversion of the SAM file to BAM file, sorting and indexing
#### Command
       for sample in `cat list.txt`
       do
               samtools view -@ 20 -S -b Mapping/${sample}.sam | samtools sort -@ 32 > Mapping/${sample}.sorted.bam

               Index BAM file
               samtools index Mapping/${sample}.sorted.bam
       done
### Filtering of mapped reads
       for sample in `cat list.txt`
       do
               samtools view -q 1 -f 0x2 -F 0x8 -b Mapping/${sample}.sorted.bam > Mapping/${sample}.filtered1.bam
       done
       
To view the output of the results use :

       samtools flagstat <bam file>
       
### To remove duplicates, the markdup command was used

#### Command

       for sample in `cat list.txt`
       do
               samtools sort -n -o Mapping/${sample}.namesort.bam Mapping/${sample}.filtered1.bam
               samtools fixmate -m Mapping/${sample}.namesort.bam Mapping/${sample}.fixmate.bam
               samtools sort -@ 4 -o Mapping/${sample}.positionsort.bam Mapping/${sample}.fixmate.bam
               samtools markdup -@4 -r Mapping/${sample}.positionsort.bam Mapping/${sample}.markdup.bam
       done

       #rmdup can also be used
       samtools rmdup SLGFSK35.sorted.bam  SLGFSK35.rdup and samtools rmdup SLGFSK36.sorted.bam  SLGFSK36.rdup
       
 ### Left align BAM
       for sample in `cat list.txt`
       do      
               cat Mapping/${sample}.markdup.bam  | bamleftalign -f hg19.chr5_12_17.fa -m 5 -c > Mapping/${sample}.leftAlign.bam
       done
       
        #-c - compressed, -m - max-iterations      
 
 ### Recalibrate read mapping qualities
       for sample in `cat list.txt`
       do
               samtools calmd -@ 4 -b Mapping/${sample}.leftAlign.bam hg19.chr5_12_17.fa > Mapping/${sample}.recalibrate.bam
       done
       
 ### Refilter read mapping qualities
       for sample in `cat list.txt`
       do
               bamtools filter -in Mapping/${sample}.recalibrate.bam -mapQuality "<=254" > Mapping/${sample}.refilter.bam
       done
       
 ## Variant calling and classification 
 http://varscan.sourceforge.net/somatic-calling.html
 
#### Installation 
	wget https://sourceforge.net/projects/varscan/files/VarScan.v2.3.9.jar

#### Command
#### Convert data to pileup
       mkdir Variants

       for sample in `cat list.txt`
       do
               samtools mpileup -f hg19.chr5_12_17.fa Mapping/${sample}.refilter.bam --min-MQ 1 --min-BQ 28 \
                       > Variants/${sample}.pileup
       done
       
### Variant calling
       varscan somatic Variants/SLGFSK-N_231335.pileup \
               Variants/SLGFSK-T_231336.pileup Variants/SLGFSK \
               --normal-purity 1 --tumor-purity 0.5 --output-vcf 1 
               
               
### Zip, index and merge vcf
varscan generates two files (snp.vcf and indel.vcf). Each file is zipped and the two zipped files were merged with bcftools

	bgzip Variants/SLGFSK.snp.vcf > Variants/SLGFSK.snp.vcf.gz
	bgzip Variants/SLGFSK.indel.vcf > Variants/SLGFSK.indel.vcf.gz
	tabix Variants/SLGFSK.snp.vcf.gz
	tabix Variants/SLGFSK.indel.vcf.gz
	bcftools merge SLGFSK.snp.vcf.gz SLGFSK.indel.vcf.gz > merged.vcf

## Variant Annotation
### Functional Annotation using SnpEff

https://pcingola.github.io/SnpEff/examples/

#### Installation
       #download jar file
       wget https://snpeff.blob.core.windows.net/versions/snpEff_latest_core.zip

       # Unzip file
       unzip snpEff_latest_core.zip
		
       #download snpEff database
       snpEff download hg19	
       
#### Command
snpEff hg19 Variants/merged.vcf > Varriants/merged.ann.vcf
       

