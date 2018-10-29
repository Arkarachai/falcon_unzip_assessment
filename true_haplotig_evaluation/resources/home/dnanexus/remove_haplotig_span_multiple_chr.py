import sys
from collections import defaultdict
import pdb


def inside_interval(chr_i, coordinate_i, coordinate_tuple):
    range_chr, range_min, range_max = coordinate_tuple
    return range_chr == chr_i and (range_max > coordinate_i > range_min)


interlace_p_h_inputfile = sys.argv[1]  # interlace_p_h_sorted_correctseg_min10

EXTENSION = 0  # 20000
FIELD_OF_INTEREST_COLUMN = 13
COORDINATE_OF_INTEREST_COLUMN = 2
CHR_OF_INTEREST_COLUMN = 12
EXPECTED_N_COLUMN = 14

chr_coordinate_storage = defaultdict(lambda: defaultdict(list))
coordinate_range = defaultdict(dict)
with open(interlace_p_h_inputfile) as fh:
    lines = fh.readlines()
    for line in lines:
        split_field = line.strip().split('\t')
        storage_key = split_field[FIELD_OF_INTEREST_COLUMN - 1].split('|')[0]
        current_chr = split_field[CHR_OF_INTEREST_COLUMN-1]
        chr_coordinate_storage[storage_key][current_chr].append(int(split_field[COORDINATE_OF_INTEREST_COLUMN-1]))

with open(interlace_p_h_inputfile) as fh:
    lines = fh.readlines()
    for line in lines:
        split_field = line.strip().split('\t')
        storage_key = split_field[FIELD_OF_INTEREST_COLUMN - 1].split('|')[0]
        if len(chr_coordinate_storage[storage_key]) == 1:
            print line.strip()
