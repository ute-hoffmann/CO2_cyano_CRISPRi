header = True
with open("Wide_DESeq_result_table.tsv") as f:
	for line in f:
		if header:
			with open("ncRNAs_results.tsv", "w") as nc_f:
				nc_f.write(line)
			with open("CDS_results.tsv", "w") as CDS_f:
				CDS_f.write(line)
			header = False
			continue
		if line[0:3] == "nc_":
			with open("ncRNAs_results.tsv", "a") as nc_f:
				split_line = line.split("\t")
				locus_tag = split_line[0].split("_")[1]
				assemb_line = ""
				for i in range(1,len(split_line)+1):
					if i == 3:
						assemb_line += locus_tag
						assemb_line += "\t"
						continue
					assemb_line += split_line[i-1]
					if i != (len(split_line)):
						assemb_line += "\t"
				nc_f.write(assemb_line)
		elif line[0:4] == "ctrl":
			print("Control sgRNAs")
			continue
		else:
			with open("CDS_results.tsv", "a") as CDS_f:
				CDS_f.write(line)
