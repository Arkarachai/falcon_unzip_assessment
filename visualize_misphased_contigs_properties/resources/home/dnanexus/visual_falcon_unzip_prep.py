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


def normalized_coor_float(snp_span_list):
    length_snp = max(snp_span_list)-min(snp_span_list)
    return map(lambda x: 100.0*(x-min(snp_span_list))/length_snp, snp_span_list)


def error_pattern(alist):
    if len(set(alist)) == 1:
        return 0, set()
    mergeword = ''.join(alist)
    switch_counter = 0
    position = 1
    end_switch_position = set()
    currentword = mergeword[0]
    for character in mergeword[1::]:
        if character != currentword:
            switch_counter += 1
            end_switch_position.add(position)
        currentword = character
        position += 1
    return switch_counter, end_switch_position


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

    jitter_value = 0

    fd_aggregate_parent_switch_mis_join = open('aggregate_parent_switch_mis_join', 'w')
    fd_aggregate_parent_switch_random_error = open('aggregate_parent_switch_random_error', 'w')
    fd_boundary_size_vs_location = open('boundary_size_vs_location', 'w')
    fd_boundary_size_vs_location.writelines('\t'.join(['class', 'n_start_position', 'n_end_position', 'distance', 'norm_distance', 'edge'])+'\n')


    for key in allkey:

        # setting and calculate error rates
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
