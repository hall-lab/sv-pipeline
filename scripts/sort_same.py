import re

def parse_name(line):
	list1 = re.split("/", line)
	list2 = re.split("\.", list1[len(list1)-1])
	return list2[0]

cram_list = "cram_list.txt"
cn_list = "cn_hist_roots.txt"
index_list = "indices.txt"
original_manta_list = "original_manta_vcfs.txt"

crams = open(cram_list, "r")
cram_list = []
for line in crams:
	cram_list.append(line.rstrip())

cns = open(cn_list, "r")
cns_dict = {}
for line in cns:
	sample_name = parse_name(line)
	cns_dict[sample_name] = line.rstrip()

indices = open(index_list, "r")
index_dict = {}
for line in indices:
	sample_name = parse_name(line)
	index_dict[sample_name] = line.rstrip()

manta = open(original_manta_list, "r")
manta_dict = {}
for line in manta:
	sample_name = parse_name(line)
	manta_dict[sample_name] = line.rstrip()

cns_out = open("cn_hist_roots_ordered.txt", "w")
index_out = open("indices_ordered.txt", "w")
manta_out = open("manta_ordered.txt", "w")
for cram in cram_list:
	sample_name = parse_name(cram)
	cns_out.write(cns_dict[sample_name])
	cns_out.write("\n")
	index_out.write(index_dict[sample_name])
	index_out.write("\n")
	manta_out.write(manta_dict[sample_name])
	manta_out.write("\n")

cns_out.close()
index_out.close()
