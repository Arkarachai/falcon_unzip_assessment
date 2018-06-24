import sys
from collections import defaultdict
import pdb


def inside_interval(chr_i, coordinate_i, coordinate_tuple):
    range_chr, range_min, range_max = coordinate_tuple
    return range_chr == chr_i and (range_max > coordinate_i > range_min)


primary_inputfile = sys.argv[1]  # 'verify_primary'
haplotig_inputfile = sys.argv[2]  # 'verify_haplotype'

EXTENSION = 0  # 20000
FIELD_OF_INTEREST_COLUMN = 13
COORDINATE_OF_INTEREST_COLUMN = 2
CHR_OF_INTEREST_COLUMN = 12
EXPECTED_N_COLUMN = 13

chr_coordinate_storage = defaultdict(lambda: defaultdict(list))
coordinate_range = defaultdict(dict)
with open(haplotig_inputfile) as fh:
    lines = fh.readlines()
    for line in lines:
        split_field = line.strip().split('\t')
        if len(split_field) != EXPECTED_N_COLUMN:
            continue
        storage_key = split_field[FIELD_OF_INTEREST_COLUMN - 1].split('|')[0]
        current_chr = split_field[CHR_OF_INTEREST_COLUMN-1]
        chr_coordinate_storage[storage_key][current_chr].append(int(split_field[COORDINATE_OF_INTEREST_COLUMN-1]))

        print '\t'.join(split_field[:FIELD_OF_INTEREST_COLUMN-1])+'\t'+storage_key+'\t'+'h'

for contig_key in chr_coordinate_storage:
    for chr_key in chr_coordinate_storage[contig_key]:
        coordinate_list = chr_coordinate_storage[contig_key][chr_key]
        max_i, min_i = max(coordinate_list), min(coordinate_list)
        origin_contig = contig_key.split('_')[0]
        coordinate_range[origin_contig][(contig_key, chr_key)] = (chr_key, min_i-EXTENSION, max_i+EXTENSION)

with open(primary_inputfile) as fh:
    lines = fh.readlines()
    for line in lines:
        split_field = line.strip().split('\t')
        if len(split_field) != EXPECTED_N_COLUMN:
            continue
        chr_i = split_field[CHR_OF_INTEREST_COLUMN-1]
        coordinate_i = int(split_field[COORDINATE_OF_INTEREST_COLUMN-1])
        origin_contig = split_field[FIELD_OF_INTEREST_COLUMN - 1].split('|')[0]
        if origin_contig in coordinate_range:
            for tuple_key in coordinate_range[origin_contig]:
                if inside_interval(chr_i, coordinate_i, coordinate_range[origin_contig][tuple_key]):
                    print '\t'.join(split_field[:FIELD_OF_INTEREST_COLUMN-1])+'\t'+tuple_key[0]+'\t'+'p'
                    break
