# Imprinted-genes-identification

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


