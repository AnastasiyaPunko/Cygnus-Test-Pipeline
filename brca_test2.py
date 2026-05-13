
"""# QC"""

fastqc BRCA_L000_R1.fq

trimmomatic SE -phred33 BRCA_L000_R1.fq BRCA_L000_R1_trimm.fq.gz ILLUMINACLIP:/home/punko_a/genom/TruSeq3-SE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36

cutadapt -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC -a "A{10}" -e 0.15 -m 36 -o BRCA_L000_R1_trimmed.fastq.gz BRCA_L000_R1_trimm.fq.gz

fastqc BRCA_L000_R1_trimmed.fastq.gz

bwa mem -R '@RG\tID:GSS5-0798\tSM:Hemato_pro_540_chip\tPL:ILLUMINA' -t 3 Homo_sapiens_assembly37.fasta BRCA_L000_R1_trimmed.fastq.gz -o BRCA.sam

samtools view -Sb BRCA.sam > BRCA.bam

samtools sort -m 3G -o BRCA_sort.bam BRCA.bam

samtools index BRCA_sort.bam



"""# Somatic"""

#Mutect2
gatk Mutect2 -R Homo_sapiens_assembly37.fasta -I BRCA_sort.bam --disable-read-filter MateOnSameContigOrNoMappedMateReadFilter -O BRCA_sort_somatic.vcf.gz

gatk FilterMutectCalls -R Homo_sapiens_assembly37.fasta -V BRCA_sort_somatic.vcf.gz -O BRCA_sort_somatic_filtered.vcf.gz


#VarScan2
samtools mpileup -f Homo_sapiens_assembly37.fasta BRCA_sort.bam > BRCA.mpileup

varscan mpileup2cns BRCA.mpileup --min-coverage 20 --min-var-freq 0.01 --p-value 0.05 --output-vcf 1 > BRCA_somatic_varscan.vcf



"""# Germline"""

#HaplotypeCaller
gatk --java-options "-Xmx4g" HaplotypeCaller -R Homo_sapiens_assembly37.fasta -I BRCA_sort.bam -O BRCA_sort_germline.vcf.gz -bamout BRCA_sort_var.bam


#VarScan2
samtools mpileup -f Homo_sapiens_assembly37.fasta -q 20 -Q 20 -B BRCA_sort.bam > BRCA_germline.mpileup

#SNP calling
varscan mpileup2snp BRCA_germline.mpileup --min-coverage 20 --min-reads2 4 --min-var-freq 0.20 --p-value 0.05 --output-vcf 1 > BRCA_varscan_snps.vcf

#InDel calling
varscan mpileup2indel BRCA_germline.mpileup --min-coverage 20 --min-reads2 4 --min-var-freq 0.20 --p-value 0.05 --output-vcf 1 > BRCA_varscan_indels.vcf

#Merge SNP+InDel
gatk MergeVcfs -I BRCA_varscan_snps.vcf -I BRCA_varscan_indels.vcf -O BRCA_varscan_all.vcf -D Homo_sapiens_assembly37.dict

#LoFreq
lofreq indelqual --dindel -f Homo_sapiens_assembly37.fasta -o BRCA_indelq.bam BRCA_sort.bam

samtools index BRCA_indelq.bam

lofreq call-parallel --pp-threads 3 -f Homo_sapiens_assembly37.fasta -o BRCA_lofreq.vcf BRCA_indelq.bam

lofreq filter -i BRCA_lofreq.vcf -o BRCA_lofreq_filtered.vcf --cov-min 20
