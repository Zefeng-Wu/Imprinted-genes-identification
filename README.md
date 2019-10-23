# Imprinted-genes-identification
1. Resquencing mapping
## reads mapping
    bwa index Brapa_sequence_v3.0.fasta 
    bwa mem -t 40 -M  -R "@RG\tID:01\tLB:01\tPL:Illumina\tPU:01\tSM:01" reference_genome/Brapa_sequence_v3.0.fasta Cole_R428-02-    R01_good_1.fq Cole_R428-02-R01_good_2.fq > 3sam/Cole_R428-02-R01_good.sam
    bwa mem -t 40 -M  -R "@RG\tID:02\tLB:02\tPL:Illumina\tPU:02\tSM:02" reference_genome/Brapa_sequence_v3.0.fasta Cole_R428-02-R02_good_1.fq Cole_R428-02-R02_good_2.fq > 3sam/Cole_R428-02-R02_good.sam

## sam to bam
    cd ../3sam
    for m in $(ls *.sam); do samtools view  -Sb $m > ../4bam/${m%.sam}.bam;done 

# keep aligned reads
    for m in $(ls *.bam); do samtools view -bF 4 $m >../5aligned_bam/$m;done 

## sort bam
    for m in $(ls *.bam); do java -jar ~/soft/picard.jar  SortSam I= $m O=../6sorted_bam/$m SORT_ORDER=coordinate ;done

### mark duplication 
    for m in $(ls *.bam); do java -jar ~/soft/picard.jar  MarkDuplicates REMOVE_DUPLICATES=false MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=8000 INPUT=$m OUTPUT= ../7dup_marked_bam/$m METRICS_FILE=../7dup_marked_bam/bam.metrics ; done

## make bam index 
    for m in $(ls *.bam); do samtools index $m;done

## vcf calling 
    gatk  HaplotypeCaller -R Brapa_sequence_v3.0.fasta -I ../7dup_marked_bam/Cole_R428-02-R01_good.bam -I ../7dup_marked_bam/Cole_R428-02-R02_good.bam -O ../8vcf/B.rapa.vcf

## selecet snp type variations
    gatk  SelectVariants -R ../reference_genome/Brapa_sequence_v3.0.fasta -variant B.rapa.vcf -O ../9snp/1SNP.vcf -select-type SNP 

2.Filter SNPs called from GATK pipline to keep the homozygou SNPs in each parental genome but different between two parental genomes. 
the same position (Assumed to be diploid).

    python 2Filter.snp_from_gatk.py
2. 

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


