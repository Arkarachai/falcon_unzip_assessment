import sys
import collections

label_column = 12
chr_column = 11
coordinate_column = 1

#fd = open(sys.argv[1])
fd = open('reduce_noheader_h_pmatch.variants')
lines = fd.readlines()
contig_name_set = set()
for line in lines:
    temp = line.strip().split('\t')
    contig_name_set.add(temp[label_column-1])
#fdd = open(sys.argv[2])
fdd = open('reduce_noheader_h.variants')


contig_location = collections.defaultdict(list)
chr_location = collections.defaultdict(set)

lines = fdd.readlines()
for line in lines:
    temp = line.strip().split('\t')
    current_contig, chromosome, coordinate = temp[label_column-1], temp[chr_column-1], int(temp[coordinate_column-1])
    contig_location[current_contig].append(coordinate)
    chr_location[current_contig].add(chromosome)

for key in contig_location:
    coordinate_list = contig_location[key]
    print '\t'.join(map(str, (key, max(coordinate_list) - min(coordinate_list), key in contig_name_set, chr_location[key])))