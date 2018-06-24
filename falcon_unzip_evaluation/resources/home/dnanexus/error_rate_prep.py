import sys
import collections
import random


def error_rate(alist):
    if len(set(alist)) == 1:
        return 0
    else:
        totallen = len(alist)
        acount = []
        for member in set(alist):
            acount.append(alist.count(member))
        return 1.0*min(acount)/totallen


def paternal_rate(alist):
    return 1.0*alist.count('d')/len(alist)


def error_pattern(alist):
    if len(set(alist)) == 1:
        return 0, set()
    mergeword = ''.join(alist)
    switch_counter = 0
    position = 1
    switch_position = set()
    currentword = mergeword[0]
    for character in mergeword[1::]:
        if character != currentword:
            switch_counter += 1
            switch_position.add(position)
        currentword = character
        position += 1
    return switch_counter, switch_position


def error_string(alist):

    acount = []
    for member in set(alist):
        acount.append((alist.count(member), member))
    correct_pattern = max(acount)[1]
    return_list = []
    for i in alist:
        if i == correct_pattern:
            return_list.append('c')
        else:
            return_list.append('i')
    return ''.join(return_list)


def error_string_breakdown(alist):
    splitlist = alist.split('c')
    itemizelist = filter(None, splitlist)
    return collections.Counter(itemizelist)


def parental_reconcile(original_parental, contig_type):
    if contig_type == 'h':
        return original_parental[0]
    else:
        if original_parental == 'mom':
            return 'd'
        elif original_parental == 'dad':
            return 'm'

if __name__ == "__main__":

    haplotype_dict = collections.defaultdict(list)
    haplotype_dict_location = collections.defaultdict(list)

    fd = open(sys.argv[1])
    lines = fd.readlines()
    for line in lines:
        temp = line.strip().split('\t')
        parental = parental_reconcile(temp[0], temp[-1])
        contig_origin = temp[-1]
        location = int(temp[1])
        haplotype_dict[(temp[-3].replace('chr', ''), temp[-2])].append((parental, contig_origin))
        haplotype_dict_location[(temp[-3].replace('chr', ''), temp[-2])].append(location)

    allkey = haplotype_dict.keys()
    allkey.sort()

    print '\t'.join(['chr', 'contig', 'N', 'error_rate', 'snp_span_length', 'parent_switch_n',
                     'contig_origin_switch_n', 'parent_pattern', 'contig_origin_pattern',
                     'pattern_switch_position', 'contig_origin_switch_position', 'parent_contig_subset',
                     'haplotig_daddy_rate', 'error_string', 'error_count', 'count_i', 'count_ii', 'count_iii',
                     'count_iiii', 'count_i>5', 'phasing_type'])

    for key in allkey:

        snp_span_list = haplotype_dict_location[key]
        snp_span_length = max(snp_span_list)-min(snp_span_list)
        parent_pattern = ''.join(map(lambda x: x[0], haplotype_dict[key]))
        contig_origin_pattern = ''.join(map(lambda x: x[1], haplotype_dict[key]))
        parent_pattern_switch_counter = error_pattern(parent_pattern)[0]
        parent_pattern_switch_position = error_pattern(parent_pattern)[1]
        contig_origin_switch_counter = error_pattern(contig_origin_pattern)[0]
        contig_origin_switch_position = error_pattern(contig_origin_pattern)[1]
        parent_contig_subset = parent_pattern_switch_position < contig_origin_switch_position
        error_string_breakdown_count = error_string_breakdown(error_string(parent_pattern))
        low_count = error_string_breakdown_count['i']+error_string_breakdown_count['ii']
        high_count = error_string(parent_pattern).count('iiiii')
        intermediate_count = error_string_breakdown_count['iii']+error_string_breakdown_count['iiii']
        if high_count > 0 and low_count > 0:
            mis_phase_type = 'mix'
        elif high_count > 0 and low_count == 0 and intermediate_count == 0:
            mis_phase_type = 'mis_join'
        elif high_count == 0 and low_count > 0 and intermediate_count == 0:
            mis_phase_type = 'random_error'
        elif intermediate_count > 0:
            mis_phase_type = 'inbetween'
        elif error_rate(parent_pattern) == 0:
            mis_phase_type = 'no_error'
        else:
            mis_phase_type = 'impossible_condition'

        print '\t'.join(map(str, [
            '\t'.join(map(str, key)), len(parent_pattern), error_rate(parent_pattern),
            snp_span_length, parent_pattern_switch_counter, contig_origin_switch_counter,
            parent_pattern, contig_origin_pattern, str(sorted(list(parent_pattern_switch_position))),
            str(sorted(list(contig_origin_switch_position))), parent_contig_subset, paternal_rate(parent_pattern),
            error_string(parent_pattern), error_string(parent_pattern).count('i'),
            error_string_breakdown_count['i'], error_string_breakdown_count['ii'],
            error_string_breakdown_count['iii'], error_string_breakdown_count['iiii'],
            error_string(parent_pattern).count('iiiii'), mis_phase_type
            ]))
