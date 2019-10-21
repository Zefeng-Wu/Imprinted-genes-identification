#!/usr/bin/env python
#coding:utf-8

"""
  Author:   --<Zefeng Wu>
  Email:    wuzefeng2008@gmail.com
  Purpose:  Extract allele-specific reads and counts at known SNP sites
  Created:  2019年10月17日
"""

import re
import pysam
import argparse

parser = argparse.ArgumentParser(description="""Split bam file and count allelic reads number based on known SNP sites""") 
parser.add_argument("-b", dest="bam_file", help="Input bam file (required)", required=True)
parser.add_argument("-s", dest="snp_file", help="Input SNP file (a vcf format is required)", required=True)
parser.add_argument("-m", "--mq", dest="map_qual", help="Set a mapping quality threshold (default:20)", default=20)
parser.add_argument("--bq", dest="base_qual", help="Set a base quality threshold (default:20)", default=20)

args = parser.parse_args()
snp_file = vars(args)['snp_file']
bam_file = vars(args)['bam_file']
map_qual = int(vars(args)['map_qual'])
base_qual = int(vars(args)['base_qual'])

# define two function for store snp info and call snp reads
def Load_snps(snp_file):
    snp_list=[]
    with open(snp_file,"r") as fh:
        for row in fh:
            if re.match("\#",row):
                continue            
            snp_info = row.strip().split()
            snp_list.append([snp_info[0], int(snp_info[1]), snp_info[3],snp_info[4]])
    fh.close()
    return sorted(snp_list)


def Reads_parse(snp, bamfh,  map_qual=20, base_qual=20, snp1_count=0, snp2_count=0):
    read_set1 = set()
    read_set2 = set()
    (chrom,bps,allele1,allele2)=snp
    found = False
    for pile in bamfh.pileup(region=chrom,start=bps-1,end=bps):
        if pile.pos == bps-1:
            pile_col=pile.pileups # <pysam.calignmentfile.PileupRead object> 
            found = True
            break    
    
    if found:
        for pile_read in pile_col:
            # check whether read is in get_reads or col_reads
            qname=pile_read.alignment.query_name
            if  not pile_read.query_position:
                continue
            
            if (pile_read.alignment.query_sequence[pile_read.query_position] == allele1 
                     and pile_read.alignment.mapping_quality >= map_qual 
                     and pile_read.alignment.query_qualities[pile_read.query_position] >= base_qual ) :
                read_set1.add(qname)
                snp1_count+=1
                
            if (pile_read.alignment.query_sequence[pile_read.query_position] == allele2 
                    and pile_read.alignment.mapping_quality >= map_qual 
                    and pile_read.alignment.query_qualities[pile_read.query_position] >= base_qual ) :
                    read_set2.add(qname)  
                    snp2_count+=1                
    return [read_set1,read_set2,snp1_count,snp2_count]




# Import bam file and snp file
snp_list = Load_snps(snp_file)
bamfh=pysam.AlignmentFile(bam_file,'rb')
outfile = open(bam_file+"_allelic_counts.txt","w")

# Count the allelic read number and store the name of SNP reads.
read_names_1 = set()
read_names_2 = set()

for snp in snp_list:
    print snp
    read_set=Reads_parse(snp,bamfh,map_qual=map_qual,base_qual= base_qual,snp1_count=0,snp2_count=0)
    outfile.write("\t".join([str(m) for m in snp])+"\t"+str(read_set[2])+"\t"+str(read_set[3])+"\n")
    read_names_1 |= read_set[0]
    read_names_2 |= read_set[1]
bamfh.close()

# Output the splited bam files based on previous reads names at SNP sites
bamfh = pysam.AlignmentFile(bam_file,'rb')   
outfh1 = pysam.Samfile(bam_file+"_ref.bam","wb",template=bamfh)
outfh2 = pysam.Samfile(bam_file+"_alt.bam","wb",template=bamfh)

index = 0
for read in bamfh.fetch(until_eof=True):
    index += 1
    if index%1000000 == 0:
        print "screened "+ str(index) + " reads of " + bam_file
    if read.query_name in read_names_1:
        outfh1.write(read)
    if read.query_name in read_names_2:
        outfh2.write(read)    
bamfh.close()
