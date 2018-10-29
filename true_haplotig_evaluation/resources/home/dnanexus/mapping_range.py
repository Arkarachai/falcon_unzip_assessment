import sys
from collections import defaultdict


def inside_interval(chr_i, coordinate_i, ref_coordinate_tuple):
    range_chr, range_min, range_max = ref_coordinate_tuple
    return range_chr == chr_i and (range_max > coordinate_i > range_min)


primary_inputfile = sys.argv[1]  # 'reduce_noheader_p.snps'
haplotig_inputfile = sys.argv[2]  # 'reduce_noheader_h.snps'

EXTENSION = 0  # 20000
FIELD_OF_INTEREST_COLUMN = 12
COORDINATE_OF_INTEREST_COLUMN = 1
CHR_OF_INTEREST_COLUMN = 11
EXPECTED_N_COLUMN = 12

chr_coordinate_storage = defaultdict(lambda: defaultdict(list))
coordinate_range = defaultdict(dict)
with open(primary_inputfile) as fh:
    lines = fh.readlines()
    for line in lines:
        split_field = line.strip().split('\t')
        #if len(split_field) < EXPECTED_N_COLUMN:
        #    pass
        storage_key = split_field[FIELD_OF_INTEREST_COLUMN - 1].split('|')[0]
        current_chr = split_field[CHR_OF_INTEREST_COLUMN - 1]
        chr_coordinate_storage[storage_key][current_chr].append(int(split_field[COORDINATE_OF_INTEREST_COLUMN-1]))

# for contig_key in chr_coordinate_storage:
#     if len(chr_coordinate_storage[contig_key]) >= 2:
#         print contig_key, chr_coordinate_storage[contig_key].keys()
for contig_key in chr_coordinate_storage:
    for chr_key in chr_coordinate_storage[contig_key]:
        coordinate_list = chr_coordinate_storage[contig_key][chr_key]
        max_i, min_i = max(coordinate_list), min(coordinate_list)
        coordinate_range[contig_key][chr_key] = (chr_key, min_i-EXTENSION, max_i+EXTENSION)


with open(haplotig_inputfile) as fh:
    lines = fh.readlines()
    for line in lines:
        split_field = line.strip().split('\t')
        #if len(split_field) < EXPECTED_N_COLUMN:
        #    pass
        chr_i = split_field[CHR_OF_INTEREST_COLUMN-1]
        coordinate_i = int(split_field[COORDINATE_OF_INTEREST_COLUMN-1])
        origin_contig = split_field[FIELD_OF_INTEREST_COLUMN - 1].split('|')[0].split('_')[0]
        if origin_contig in coordinate_range:
            if chr_i in coordinate_range[origin_contig]:
                if inside_interval(chr_i, coordinate_i, coordinate_range[origin_contig][chr_i]):
                    print line.strip()

            # else:
            #     #print 'not in range', (chr_i, coordinate_i, coordinate_range[origin_contig])
            #     if chr_i == coordinate_range[origin_contig][0]:
            #         #print min(abs(coordinate_i-coordinate_range[origin_contig][1]),abs(coordinate_i-coordinate_range[origin_contig][2]))
            #         print 'match_chr_not_range',line.strip()#, coordinate_range[origin_contig]
