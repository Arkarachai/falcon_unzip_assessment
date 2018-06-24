import sys
filename = sys.argv[1]

pseudo_autosomal = [(10001, 2781479), (155701383, 156030895)]
fd = open(filename)
lines = fd.readlines()
for line in lines:
    temp = line.strip().split()
    if temp[11] != 'chrX':
        print line.strip()
    else:
        if (int(temp[1]) > pseudo_autosomal[0][0] and int(temp[1]) < pseudo_autosomal[0][1]) or \
            (int(temp[1]) > pseudo_autosomal[1][0] and int(temp[1]) < pseudo_autosomal[1][1]):
            print line.strip()