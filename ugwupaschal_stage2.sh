#TASK 1

#Download the reference genome
mkdir ~/dc_workshop
cd ~/dc_workshop
mkdir ref_genome
wget https://zenodo.org/record/2582555/files/SLGFSK-N_231335_r1_chr5_12_17.fastq.gz
gunzip SLGFSK-N_231335_r1_chr5_12_17.fastq.gz
wget https://zenodo.org/record/2582555/files/SLGFSK-N_231335_r2_chr5_12_17.fastq.gz
gunzip SLGFSK-N_231335_r2_chr5_12_17.fastq.gz

#Download a set of trimmed FastQ files to work with
curl -L -o ref_genome/hg19.chr5_12_17.fa.gz https://zenodo.org/record/2582555/files/hg19.chr5_12_17.fa.gz

#Create directories for the results that will be generated as part of this workflow
mkdir results

#Index the reference genome for use by BWA
sudo apt-get -y install bwa
bwa index ref_genome/hg19.chr5_12_17.fa.gz

#Align reads to reference genome
sudo apt-get -y install samtools
bwa mem ref_genome/hg19.chr5_12_17.fa.gz SLGFSK-N_231335_r1_chr5_12_17.fastq SLGFSK-N_231335_r2_chr5_12_17.fastq > results/SLGFSK-N_231335.sam

#Convert the SAM file to BAM format using samtools program
ls -lh
samtools view -S -b results/SLGFSK-N_231335.sam > results/SLGFSK-N_231335.bam

#Next, we sort the the BAM file using the sort command from samtools
samtools  view results/SLGFSK-N_231335.bam | head -n 5
samtools sort results/SLGFSK-N_231335.bam -o results/SLGFSK-N_231335.sorted.bam
samtools  view results/SLGFSK-N_231335.sorted.bam | head -n 5
samtools flagstat results/SLGFSK-N_231335.sorted.bam

#Calculate the read coverage of positions in the genome
sudo apt-get -y install bcftools
gunzip ref_genome/hg19.chr5_12_17.fa.gz
bcftools mpileup -O b -o results/SLGFSK-N_231335.bcf -f ref_genome/hg19.chr5_12_17.fa --threads 8 -q 20 -Q 30 results/SLGFSK-N_231335.sorted.bam

#Detect the single nucleotide variants (SNVs)
bcftools call --ploidy 1 -m -v -o results/SLGFSK-N_231335.vcf results/SLGFSK-N_231335.bcf 

#Filter and report the SNV variants in variant calling format (VCF)
vcfutils.pl varFilter results/SLGFSK-N_231335.vcf  > results/SLGFSK-N_231335.vcf

#Explore the VCF format
less -S results/SLGFSK-N_231335.vcf

#Use the grep and wc commands you have learned to assess how many variants are in the vcf file.
grep -v "#" results/SLGFSK-N_231335.vcf | wc -l

#Assess the alignment (visualization) - optional step
samtools index results/SLGFSK-N_231335.sorted.bam

#Viewing with tview
samtools tview results/SLGFSK-N_231335.sorted.bam ref_genome/hg19.chr5_12_17.fa


#TASK 2

#Download the reference genome
mkdir ~/dc_workshop
cd ~/dc_workshop
mkdir ref_genome
wget https://zenodo.org/record/2582555/files/SLGFSK-T_231336_r1_chr5_12_17.fastq.gz
gunzip SLGFSK-T_231336_r1_chr5_12_17.fastq.gz
wget https://zenodo.org/record/2582555/files/SLGFSK-T_231336_r2_chr5_12_17.fastq.gz
gunzip SLGFSK-T_231336_r2_chr5_12_17.fastq.gz

#Download a set of trimmed FastQ files to work with
curl -L -o ref_genome/hg19.chr5_12_17.fa.gz https://zenodo.org/record/2582555/files/hg19.chr5_12_17.fa.gz

#Create directories for the results that will be generated as part of this workflow
mkdir results

#Index the reference genome for use by BWA
sudo apt-get -y install bwa
bwa index ref_genome/hg19.chr5_12_17.fa.gz

#Align reads to reference genome
sudo apt-get -y install samtools
bwa mem ref_genome/hg19.chr5_12_17.fa.gz SLGFSK-T_231336_r1_chr5_12_17.fastq SLGFSK-T_231336_r2_chr5_12_17.fastq > results/SLGFSK-T_231336.sam

#Convert the SAM file to BAM format using samtools program
ls -lh
samtools view -S -b results/SLGFSK-T_231336.sam > results/SLGFSK-T_231336.bam

#Next, we sort the the BAM file using the sort command from samtools
samtools  view results/SLGFSK-T_231336.bam | head -n 5
samtools sort results/SLGFSK-T_231336.bam -o results/SLGFSK-T_231336.sorted.bam
samtools  view results/SLGFSK-T_231336.sorted.bam | head -n 5
samtools flagstat results/SLGFSK-T_231336.sorted.bam

#Calculate the read coverage of positions in the genome
sudo apt-get -y install bcftools
gunzip ref_genome/hg19.chr5_12_17.fa.gz
bcftools mpileup -O b -o results/SLGFSK-T_231336.bcf -f ref_genome/hg19.chr5_12_17.fa --threads 8 -q 20 -Q 30 results/SLGFSK-T_231336.sorted.bam

#Detect the single nucleotide variants (SNVs)
bcftools call --ploidy 1 -m -v -o results/SLGFSK-T_231336.vcf results/SLGFSK-T_231336.bcf 

#Filter and report the SNV variants in variant calling format (VCF)
vcfutils.pl varFilter results/SLGFSK-T_231336.vcf  > results/SLGFSK-T_231336.vcf

#Explore the VCF format
less -S results/SLGFSK-T_231336.vcf

#Use the grep and wc commands you have learned to assess how many variants are in the vcf file.
grep -v "#" results/SLGFSK-T_231336.vcf | wc -l

#Assess the alignment (visualization) - optional step
samtools index results/SLGFSK-T_231336.sorted.bam

#Viewing with tview
samtools tview results/SLGFSK-T_231336.sorted.bam ref_genome/hg19.chr5_12_17.fa
