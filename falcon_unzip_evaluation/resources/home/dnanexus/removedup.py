import sys
import collections

input_file = sys.argv[1]
dup_column = map(int, sys.argv[2].split(','))
dup_column = map(lambda x: x-1, dup_column)

storage_list = []
with open(input_file) as fd:
    lines = fd.readlines()
    for line in lines:
        temp = line.strip().split()
        dup_value_tmp = [temp[i] for i in dup_column]
        storage_list.append(tuple(dup_value_tmp))

dup_collector = [item for item, count in collections.Counter(storage_list).items() if count > 1]
dup_collector = set(dup_collector)

with open(input_file) as fd:
    lines = fd.readlines()
    for line in lines:
        temp = line.strip().split()
        dup_value_tmp = [temp[i] for i in dup_column]
        if tuple(dup_value_tmp) not in dup_collector:
            print line.strip()
