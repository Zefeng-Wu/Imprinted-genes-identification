output = open("2.5SNP_filtered_modified.vcf","w")
with open("2SNP_filtered.vcf") as fh:
    for row in fh:
        if row.startswith("#"):
            output.write(row)
        else:
            data = row.strip().split("\t")
            if data[9].split(":")[0]=="0/0":
                output.write(row)
            else:
                R01 = data[4]
                R02 = data[3]
                data[3] = R01
                data[4] = R02
                output.write("\t".join(data)+"\n")
            
print 'OK' 
    