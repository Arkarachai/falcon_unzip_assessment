import sys
import collections
import random
import pandas
import pdb

REVERSE_PARENT_DICT = {'d': 'm', 'm': 'd'}


def _write_file(string, file):
    file.writelines(string+'\n')


def dataframe_filtering(df, parent, contig_origin):
    newdf = df[(df['parent'] == parent) & (df['contig_origin'] == contig_origin)]
    if newdf.shape[0] > 0:
        return newdf.to_string(index=False, header=False, columns=['chromosome', 'start', 'stop'])
    else:
        return '{0}\t{1}\t{2}'.format('chr1', 0, 1)


def making_visual_track(contig_id, contig_ucsc_dict):

    df = pandas.DataFrame(contig_ucsc_dict[contig_id])
    df.columns = ['chromosome', 'start', 'stop', 'parent', 'contig_origin', 'contig_id']
    min_coordinate = min(df['start'])
    max_coordinate = max(df['stop'])
    with open(contig_id+'.ucsc.track', 'w') as fdd:
        _write_file('browser position {0}:{1}-{2}'.format(df['chromosome'][0], min_coordinate, max_coordinate), fdd)
        _write_file('track name={0}_haplotig  description="haplotype contig" color=0,0,0,'.format(contig_id), fdd)
        _write_file('{0}\t{1}\t{2}'.format(df['chromosome'][0], min_coordinate, max_coordinate), fdd)
        _write_file('track name=Paternal_haplotig description="Paternal SNPs in haplotig" color=0,0,255,', fdd)
        _write_file(dataframe_filtering(df, 'd', 'h'), fdd)
        _write_file('track name=Maternal_primary description="Maternal SNPs in primary contig" color=0,0,255,', fdd)
        _write_file(dataframe_filtering(df, 'm', 'p'), fdd)
        _write_file('track name=Maternal_haplotig description="Maternal SNPs in haplotig" color=255,0,00,', fdd)
        _write_file(dataframe_filtering(df, 'm', 'h'), fdd)
        _write_file('track name=Paternal_primary description="Paternal SNPs in primary contig" color=255,0,00,', fdd)
        _write_file(dataframe_filtering(df, 'd', 'p'), fdd)


def parental_reconcile(original_parental, contig_type):
    if contig_type == 'h':
        return original_parental[0]
    else:
        if original_parental == 'mom':
            return 'd'
        elif original_parental == 'dad':
            return 'm'

if __name__ == "__main__":

    interested_contigs = set(sys.argv[2].replace(' ', '').split(','))
    print('interested_contigs ', interested_contigs)
    haplotype_dict = collections.defaultdict(list)
    haplotype_dict_location = collections.defaultdict(list)
    contig_ucsc_dict = collections.defaultdict(list)

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

    intersect_keys = []
    for key in allkey:
        if key[1] in interested_contigs:
            intersect_keys.append(key)

    print('All present keys:')
    print(map(lambda x: x[1], intersect_keys))

    for key in intersect_keys:
        snp_span_list = haplotype_dict_location[key]
        snp_span_length = max(snp_span_list)-min(snp_span_list)
        parent_pattern = ''.join(map(lambda x: x[0], haplotype_dict[key]))
        contig_origin_pattern = ''.join(map(lambda x: x[1], haplotype_dict[key]))

        for i in range(len(snp_span_list)):  # note that we can do this before for key in allkey loop
            simplify_parent = parent_pattern[i]
            contig_origin = contig_origin_pattern[i]
            if contig_origin == 'p':
                original_parent = REVERSE_PARENT_DICT[simplify_parent]
            else:
                original_parent = simplify_parent
            contig_ucsc_dict[key[1]].append(['chr'+str(key[0]), snp_span_list[i], snp_span_list[i]+1,
                                             original_parent, contig_origin, key[1]])

    for intersect_key in intersect_keys:
        making_visual_track(intersect_key[1], contig_ucsc_dict)
