from __future__ import print_function

import sys


def genotype_column(record):

    tempgenotype = map(int, record.split(':')[0].split('/'))
    tempgenotype.sort()
    return tempgenotype

count_heterozygous_snp = 0
count_informative_snp = 0

for line in sys.stdin:
    if line.strip().startswith('#'):  # skip vcf header
        continue
    column_split = line.strip().split('\t')
    if column_split[9][0] == '.' or column_split[10][0] == '.' or column_split[11][0] == '.':  # skip empty genotype
        continue

    kid_info = column_split[9]
    dad_info = column_split[10]
    mom_info = column_split[11]

    genotype_kid = genotype_column(kid_info)
    genotype_dad = genotype_column(dad_info)
    genotype_mom = genotype_column(mom_info)

    if genotype_kid[0] == genotype_kid[1]:  # remove homozygote
        continue

    count_heterozygous_snp += 1

    if genotype_dad == genotype_mom:  # remove dad == mom
        continue

    if (genotype_kid[0] in genotype_dad and genotype_kid[1] in genotype_mom) or \
            (genotype_kid[1] in genotype_dad and genotype_kid[0] in genotype_mom):
        pass
    else:  # get rid of germline mutation
        continue

    # making dict of genotype
    alternative_form = column_split[4].split(',')
    num_alternative = len(alternative_form)
    if num_alternative == 1:
        decode_dict = {0: column_split[3], 1: column_split[4]}
    elif num_alternative == 2:
        decode_dict = {0: column_split[3], 1: alternative_form[0], 2: alternative_form[1]}
    elif num_alternative == 3:
        decode_dict = {0: column_split[3], 1: alternative_form[0], 2: alternative_form[1], 3: alternative_form[2]}
    elif num_alternative == 4:
        decode_dict = {0: column_split[3], 1: alternative_form[0], 2: alternative_form[1], 3: alternative_form[2], 4: alternative_form[3]}
    else:  # get rid of situation when there are more than 4 variant
        continue

    possible_P = set(genotype_kid).intersection(set(genotype_dad))
    possible_M = set(genotype_kid).intersection(set(genotype_mom))

    if len(possible_P) == 1 and len(possible_M) == 1:
        P_descence = decode_dict[int(list(possible_P)[0])]
        M_descence = decode_dict[int(list(possible_M)[0])]
    elif len(possible_P) == 2 and len(possible_M) == 1:
        M_descence = decode_dict[int(list(possible_M)[0])]
        P_descence = decode_dict[int(list(possible_P-possible_M)[0])]
    elif len(possible_P) == 1 and len(possible_M) == 2:
        P_descence = decode_dict[int(list(possible_P)[0])]
        M_descence = decode_dict[int(list(possible_M-possible_P)[0])]
    else:
        continue
    allgenotype = []
    allgenotype.extend(genotype_kid)
    allgenotype.extend(genotype_dad)
    allgenotype.extend(genotype_mom)

    count_informative_snp += 1

    print('{0}\t{1}\t{2}\t{3}'.format(P_descence, M_descence, '\t'.join(column_split), allgenotype))

print("count_heterozygous_snp ", count_heterozygous_snp, file=sys.stderr)
print("count_informative_snp ", count_informative_snp, file=sys.stderr)
