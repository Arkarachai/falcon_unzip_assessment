import sys
from collections import defaultdict

input_file = sys.argv[1]
col_num = int(sys.argv[2])
frequency_min = int(sys.argv[3])

frequency_dict = defaultdict(int)

fd = open(input_file)
lines = fd.xreadlines()
for line in lines:
    temp = line.strip().split()
    frequency_dict[temp[col_num-1]] += 1

fd = open(input_file)
lines = fd.xreadlines()
for line in lines:
    temp = line.strip().split()
    if frequency_dict[temp[col_num-1]] >= frequency_min:
        print line.strip()

