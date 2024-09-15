sgRNA_dict = {}

with open("Synechocystis_v2_trimmed.fasta") as sgRNA_library:
	for line in sgRNA_library:
		if line[0] == ">":
			target_name = line.strip(">").split("|")[0]
		else:
			if target_name in sgRNA_dict.keys():
				sgRNA_dict[target_name] += 1
			else: 
				sgRNA_dict[target_name] = 1
				
with open("number_sgRNAs_per_target.tsv", "w") as new_f:
	for target in sgRNA_dict.keys():
		new_f.write(str(target) + "\t" + str(sgRNA_dict[target]) + "\n")	
