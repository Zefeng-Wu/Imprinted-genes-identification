# Imprinted-genes-identification
## 1. Resquencing data analysis and call SNPs
 ### 1.1 reads mapping
    bwa index Brapa_sequence_v3.0.fasta 
    bwa mem -t 40 -M  -R "@RG\tID:01\tLB:01\tPL:Illumina\tPU:01\tSM:01" reference_genome/Brapa_sequence_v3.0.fasta Cole_R428-02-    R01_good_1.fq Cole_R428-02-R01_good_2.fq > 3sam/Cole_R428-02-R01_good.sam
    bwa mem -t 40 -M  -R "@RG\tID:02\tLB:02\tPL:Illumina\tPU:02\tSM:02" reference_genome/Brapa_sequence_v3.0.fasta Cole_R428-02-R02_good_1.fq Cole_R428-02-R02_good_2.fq > 3sam/Cole_R428-02-R02_good.sam
 ### 1.2 Sam to bam
    cd ../3sam
    for m in $(ls *.sam); do samtools view  -Sb $m > ../4bam/${m%.sam}.bam;done 
 ### 1.3 Keep aligned reads
    for m in $(ls *.bam); do samtools view -bF 4 $m >../5aligned_bam/$m;done 
 ### 1.4 Sort bam
    for m in $(ls *.bam); do java -jar ~/soft/picard.jar  SortSam I= $m O=../6sorted_bam/$m SORT_ORDER=coordinate ;done
 ### 1.5 Mark duplications 
    for m in $(ls *.bam); do java -jar ~/soft/picard.jar  MarkDuplicates REMOVE_DUPLICATES=false MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=8000 INPUT=$m OUTPUT= ../7dup_marked_bam/$m METRICS_FILE=../7dup_marked_bam/bam.metrics ; done
 ### 1.6 Make bam index 
    for m in $(ls *.bam); do samtools index $m;done
 ### 1.7 Varients calling 
    gatk  HaplotypeCaller -R Brapa_sequence_v3.0.fasta -I ../7dup_marked_bam/Cole_R428-02-R01_good.bam -I ../7dup_marked_bam/Cole_R428-02-R02_good.bam -O ../8vcf/B.rapa.vcf
 ### 1.8 Selecet snp type variations
    gatk  SelectVariants -R ../reference_genome/Brapa_sequence_v3.0.fasta -variant B.rapa.vcf -O ../9snp/1SNP.vcf -select-type SNP 
 ### 1.9. Filter SNPs called from GATK pipline to keep the homozygou SNPs in each parental genome but different between two parental genomes.(Assumed to be diploid).
    python 2Filter.snp_from_gatk.py
 ### 1.10. Modified the SNP position in the vcf (important!) 
    python Filter.snp_modified.py
## 2. RNA-Seq data analysis and call SNPs
 ### 2.1 Reads mapping to reference genome using Hisat2 or other mapping tools.
    for m in $(ls *_1.fq.gz); do if [ ! -f ../2Hisat_result/${m%_good_1.fq.gz}.sam ]; then echo $m; hisat2 -p 40 -x  ../../rapa_data/reference_genome/Brapa_sequence_v3.0 -1 $m -2 ${m%1.fq.gz}2.fq.gz --no-mixed  --no-discordant  -S ../2Hisat_result/${m%_good_1.fq.gz}.sam; fi; done
 ### 2.2 Extract uniqe mapping reads
    for m in $(ls *.sam); do echo $m; grep 'NH:i:1' $m > ../3uniq_mapping/${m::-4}.uniq.sam;done
 ### 2.3 Sam to bam and add header to bam
    for m in $(ls *.sam); do echo $m; samtools view -S -b $m -T ../../rapa_data/reference_genome/Brapa_sequence_v3.0.fasta > ../4bam/${m%.uniq.sam}.bam;done
 ### 2.4 Sort bam files
    for m in $(ls *.bam); do samtools sort $m -o ../5sorted_bam/${m%.bam}.sorted.bam -O bam -@ 40 ;done
 ### 2.5 Merge all the RNA-Seq bam files    
    samtools merge ../all.bam *.bam -@ 40
 ### 2.6 Call SNPs with mplieup and bcftools
    samtools mpileup -v -u -I 6merge_bam/all.bam -o 7mpileup_out/1RNA_merge.vcf --reference ../rapa_data/reference_genome/Brapa_sequence_v3.0.fasta
    bcftools call -mv 1RNA_merge.vcf  -o 2calls.bcf --threads 40
 ### 2.7 Filter SNPs by sequence reads depth (threshold: 10)
    awk 'm=split($8,u,";"){y=split(u[1],x,"=")}{if(x[2]>=10) print $0}' 2calls.bcf  >3RNA_merge.vcf  
## 3.Obtain the common SNPs between RNA-Seq SNPs and resequenced SNPs, and keep only those SNPs within reference gene sets (from GFF file).
    Rscript obtining_common_snp.R
## 4. Count the allelic reads from RNA-Seq alignment file (bam formatted) at the SNPs sites derived from  the previous step.
    for m in $(ls ../5sorted_bam/*.bam);do echo $m; python ASE_count.py -b $m -s 4target.snp.txt --mq 0 --bq 0;done
    
    usage: ASE_count.py [-h] -b BAM_FILE -s SNP_FILE [-m MAP_QUAL]
                    [--bq BASE_QUAL]

    Split bam file and count allelic reads number based on known SNP sites

    optional arguments:
      -h, --help            show this help message and exit
      -b BAM_FILE           Input bam file (required)
      -s SNP_FILE           Input SNP file (a vcf formatted file or at least the first 5 columns of vcf file is required)
      -m MAP_QUAL, --mq MAP_QUAL
                        Set a mapping quality threshold (default:20)
      --bq BASE_QUAL        Set a base quality threshold (default:20)
 ## 5. Merge mutiple files
    awk '{for(i=1;i<=6;i++)printf("%s\t",$i);for(i=7;i<=NF;i+=6)printf("%s\t",$(i+4)"\t"$(i+5));print ""}' all.allel.counts.txt > 2all_allele_counts
    ls -1  ../5sorted_bam/*.txt  | tr '\n' '\0' |xargs -0 -n 1 basename | grep -wo "T.." | xargs # get the sample name

