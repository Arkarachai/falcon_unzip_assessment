import sys
filename = sys.argv[1]

pseudo_autosomal = [(60001, 2699520), (154931044, 155260560)]
fd = open(filename)
lines = fd.readlines()
for line in lines:
    temp = line.strip().split()
    if temp[11] == 'X' or temp[11] == 'chrX':
        if (int(temp[1]) > pseudo_autosomal[0][0] and int(temp[1]) < pseudo_autosomal[0][1]) or \
                (int(temp[1]) > pseudo_autosomal[1][0] and int(temp[1]) < pseudo_autosomal[1][1]):
            print line.strip()
    else:
        print line.strip()