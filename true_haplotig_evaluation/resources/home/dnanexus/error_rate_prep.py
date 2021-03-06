from __future__ import print_function
import sys
import collections
import random


def normalized_coor_float(snp_span_list):
    length_snp = max(snp_span_list)-min(snp_span_list)
    return map(lambda x: 100.0*(x-min(snp_span_list))/length_snp, snp_span_list)


def error_pattern(alist):
    if len(set(alist)) == 1:
        return 0, set()
    switch_counter = 0
    position = 1
    switch_position = set()
    currentword = alist[0]
    for character in alist[1::]:
        if character != currentword:
            switch_counter += 1
            switch_position.add(position)
        currentword = character
        position += 1
    return switch_counter, switch_position


def error_string(alist, parent_label):

    correct_pattern = parent_label
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


def transfer_label_advance(alist, jitter_value):
    '''
    value=random.random()
    while value<0.7 and value>0.3:
        value = random.random()
    '''
    value = jitter_value
    if value < 0.5:
        value = 1-value

    blist = []
    for i in alist:
        if i == 'd':
            blist.append(value)
        elif i == 'm':
            blist.append(1-value)
    return blist


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
    parent_label = 'd' if sys.argv[2] == 'paternal' else 'm'
    print('correct label = {0}'.format(parent_label), file=sys.stderr)
    lines = fd.readlines()
    for line in lines:
        temp = line.strip().split('\t')
        parental = parental_reconcile(temp[0], 'h')
        location = int(temp[1])
        haplotype_dict[(temp[-2].replace('chr', ''), temp[-1])].append(parental)
        haplotype_dict_location[(temp[-2].replace('chr', ''), temp[-1])].append(location)

    allkey = haplotype_dict.keys()
    allkey.sort()

    jitter_value = 0

    fd_aggregate_parent_switch_mis_join = open('aggregate_parent_switch_mis_join', 'w')
    fd_aggregate_parent_switch_random_error = open('aggregate_parent_switch_random_error', 'w')
    fd_boundary_size_vs_location = open('boundary_size_vs_location', 'w')
    fd_boundary_size_vs_location.writelines('\t'.join(['class', 'n_start_position', 'n_end_position', 'distance', 'norm_distance', 'edge'])+'\n')


    print('\t'.join(['chr', 'contig', 'N', 'error_rate', 'snp_span_length', 'parent_switch_n',
                     'parent_pattern', 'pattern_switch_position',
                     'error_string', 'error_count', 'count_i', 'count_ii', 'count_iii',
                     'count_iiii', 'count_i>5', 'phasing_type']))

    for key in allkey:

        snp_span_list = haplotype_dict_location[key]
        snp_span_length = max(snp_span_list)-min(snp_span_list)
        parent_pattern = ''.join(haplotype_dict[key])
        parent_pattern_switch_counter, parent_pattern_switch_position = error_pattern(parent_pattern)
        error_string_sequence = error_string(parent_pattern, parent_label)
        error_string_breakdown_count = error_string_breakdown(error_string_sequence)
        low_count = error_string_breakdown_count['i']+error_string_breakdown_count['ii']
        high_count = error_string_sequence.count('iiiii')
        intermediate_count = error_string_breakdown_count['iii']+error_string_breakdown_count['iiii']
        if high_count > 0 and low_count > 0:
            mis_phase_type = 'mix'
        elif high_count > 0 and low_count == 0 and intermediate_count == 0:
            mis_phase_type = 'mis_join'
        elif high_count == 0 and low_count > 0 and intermediate_count == 0:
            mis_phase_type = 'random_error'
        elif intermediate_count > 0:
            mis_phase_type = 'inbetween'
        elif error_string_sequence.count('i') == 0:
            mis_phase_type = 'no_error'
        else:
            mis_phase_type = 'impossible_condition'

        print('\t'.join(map(str, [
            '\t'.join(map(str, key)), len(parent_pattern), (error_string_sequence.count('i')*1.0)/len(parent_pattern),
            snp_span_length, parent_pattern_switch_counter,
            parent_pattern,  str(sorted(list(parent_pattern_switch_position))),
            error_string_sequence, error_string_sequence.count('i'),
            error_string_breakdown_count['i'], error_string_breakdown_count['ii'],
            error_string_breakdown_count['iii'], error_string_breakdown_count['iiii'],
            error_string_breakdown_count['iiiii'], mis_phase_type
            ])))

        # make switching contig plot for i5 up
        if mis_phase_type == 'mis_join':
            for i in zip(transfer_label_advance(parent_pattern, jitter_value), normalized_coor_float(snp_span_list)):
                fd_aggregate_parent_switch_mis_join.writelines('\t'.join([key[1], str(i[0]), str(i[1])])+'\n')
            jitter_value += 0.005

        # make switching contig plot for i1 and i2
        if mis_phase_type == 'random_error':
            for i in zip(transfer_label_advance(parent_pattern, jitter_value), normalized_coor_float(snp_span_list)):
                fd_aggregate_parent_switch_random_error.writelines('\t'.join([key[1], str(i[0]), str(i[1])])+'\n')
            jitter_value += 0.005


        # get boundary position, size, normalized size
        if mis_phase_type == 'mis_join' or mis_phase_type == 'random_error':
            start_switch_position_normalize = [normalized_coor_float(snp_span_list)[index-1] for index in sorted(list(parent_pattern_switch_position))]
            end_switch_position_normalize = [normalized_coor_float(snp_span_list)[index] for index in sorted(list(parent_pattern_switch_position))]
            for start_norm_position, end_norm_position, index in zip(start_switch_position_normalize, end_switch_position_normalize, sorted(list(parent_pattern_switch_position))):
                if start_norm_position == 0 or end_norm_position == 100:
                    edge_label = "edge"
                else:
                    edge_label = "not_edge"
                norm_distance = normalized_coor_float(snp_span_list)[index] - normalized_coor_float(snp_span_list)[index-1]
                true_distance = snp_span_list[index] - snp_span_list[index-1]
                fd_boundary_size_vs_location.writelines('{0}\t{1}\t{2}\t{3}\t{4}\t{5}'.format(
                    mis_phase_type, start_norm_position, end_norm_position, true_distance, norm_distance, edge_label)+'\n')


    fd_aggregate_parent_switch_mis_join.close()
    fd_aggregate_parent_switch_random_error.close()
    fd_boundary_size_vs_location.close()
