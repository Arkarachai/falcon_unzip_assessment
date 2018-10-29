import sys
import collections

illumina_data = sys.argv[1]  # classifier_placer_holder
pacbio_data = sys.argv[2]  # reduce_noheader_h_pmatch.snps, reduce_noheader_p_pmatch.snps

Illumina_dict = collections.defaultdict(list)
fd = open(illumina_data)
lines = fd.xreadlines()
for line in lines:
    temp = line.strip().split('\t')
    Illumina_dict[(temp[2], temp[3])] = (temp)

fd.close()


fdd = open(pacbio_data)
lines = fdd.xreadlines()
for line in lines:
    temp = line.strip().split('\t')
    ref, alt, chrom, start = temp[1], temp[2], temp[10], temp[0]
    if ref == '.' or alt == '.':  # screen for SNPs only
        continue
    if (chrom, start) in Illumina_dict:
        Illumina_record = Illumina_dict[(chrom, start)]
        father_illumina, mother_illumina = Illumina_record[0], Illumina_record[1]
        if (alt == father_illumina) and (alt != mother_illumina):
            inherit = 'dad'
            print inherit + '\t' + line.strip()
        elif (alt == mother_illumina) and (alt != father_illumina):
            inherit = 'mom'
            print inherit + '\t' + line.strip()
        elif (alt == mother_illumina) and (alt == father_illumina):
            print 'problematic_match_both'+'\t'+line.strip()+'\t'+'#'+'\t'.join(Illumina_record)
        elif (alt != mother_illumina) and (alt != father_illumina):
            print 'problematic_no_match'+'\t'+line.strip()+'\t'+'#'+'\t'.join(Illumina_record)
