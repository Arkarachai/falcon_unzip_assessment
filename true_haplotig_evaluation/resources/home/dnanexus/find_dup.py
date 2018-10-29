import sys
fd = open(sys.argv[1])
lines = fd.readlines()
previous_chr = ''
for line in lines:
    temp = line.strip().split()
    current_chr = (temp[1], temp[3])
    current_line = line.strip()
    if current_chr != previous_chr:
        print line.strip()
    previous_chr = current_chr

print line.strip()