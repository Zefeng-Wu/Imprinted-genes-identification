output = open("2SNP_filtered.vcf","w")
with open("1SNP.vcf") as fh:
    for row in fh:
        if row.startswith("#"):
            output.write(row)
        else:
            data = row.strip().split("\t")
            if len(data[3])==len(data[4]): #(bi-allelic snp)
                R01_allele1 = data[9].split(":")[0].split("/")[0]
                R01_allele2 = data[9].split(":")[0].split("/")[1]
                R02_allele1 = data[10].split(":")[0].split("/")[0]
                R02_allele2 = data[10].split(":")[0].split("/")[1]
                if R01_allele1==R01_allele2 and R02_allele1==R02_allele2 and R01_allele1!=R02_allele1:
                    if R01_allele1!="." and R02_allele1!=".":
                        output.write("\t".join(data)+"\n") 
print 'OK' 
    