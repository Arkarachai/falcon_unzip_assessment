import sys
from collections import Counter

input_name = sys.argv[1]
contig_name_column = int(sys.argv[2])-1
fd = open(input_name)
lines = fd.readlines()
duplist_collector = []
dupset = set()
overlapset = set()
previous_contig = 'PLACEHODLER_NAME'
for line in lines:
    temp = line.strip().split()
    current_contig = (temp[contig_name_column])
    if current_contig != previous_contig:
        duplist_collector.append(previous_contig)  # last contig name for each chunk
    previous_contig = current_contig
member_counter = Counter(duplist_collector)
for member in member_counter:
    if member_counter[member] > 1:
        dupset.add(member)
for dup_member in dupset:
    indices = [i for i, x in enumerate(duplist_collector) if x == dup_member]
    b_indices, e_indices = min(indices), max(indices)
    overlap_member = duplist_collector[b_indices:e_indices+1]
    overlapset.update(overlap_member)

fd = open(input_name)
lines = fd.readlines()
for line in lines:
    temp = line.strip().split()
    current_contig = (temp[contig_name_column])
    if current_contig not in overlapset:
        print line.strip()
