import os, csv
import zlib 
import argparse
from Bio.SeqIO import parse 


def allele_hash(seq):
    """CRC32 unsigned integer from allele nucleotide sequence.
    The "& 0xffffffff" is used to generate an unsigned integer and
    will generate the same number for both Python 2 and 3
    (https://docs.python.org/2/library/zlib.html#zlib.crc32).

    seq (str): nucleotide string
    returns: int: CRC32 checksum as unsigned 32bit integer
    """
    seq = str(seq).encode()
    return (zlib.crc32(seq) & 0xffffffff)




sample_file = snakemake.input['assembly']
sample = os.path.basename(sample_file).replace('.fa', '').replace('_contigs_filter', '')
print("sample name: %s"%(sample))


# contigsInfo = snakemake.input['contigsInfo']
contigsInfo_dir = os.path.dirname(snakemake.input['schema'])
contigsInfo = contigsInfo_dir+'/'+snakemake.params['genus']+'/results/results_contigsInfo.tsv'



alleles_position = {}

with open(contigsInfo) as f:
    header = f.readline().rstrip().split('\t')
    for line in f:
        line = line.rstrip().split('\t')
        print(line[:3])
        name = line[0].replace('_contigs_filter', '')
        found= False
        if name == sample:
            found= True
            print("sample name in matrice: %s"%(name))
            print("at this step sample name should be identical")
            name = line[0].replace('_contigs_filter', '') # catch the name value
            print(name)
            i=0
            for position in line[1:]:
                position = position.split('&')
                i+=1
                if len(position) == 3:
                    ref = position[0]
                    strand = position[2]
                    if ref not in alleles_position:
                        alleles_position[ref] = {'1':{}, '0':{}}
                    #n the reverse strand (represented by a 0 signal). 1 means that the CDS is encoded in the direct strand.
                    alleles_position[ref][strand][header[i]] = position[1]
                else:
                    print('no allele found for this gene')                   
            break

if not found:
    raise ValueError('A very bad thing happened!!! Sample is not found in chewbbaca output')


hash_dict = {}
for seq_record in parse(sample_file, "fasta"):
    ref = seq_record.id
    if ref in alleles_position:
        print(len(seq_record.seq))
        for gene, position in alleles_position[ref]['1'].items():
            position_start = int(position.split('-')[0])
            position_end = int(position.split('-')[1])
            seq = seq_record.seq[position_start:position_end]
            seq_hash = allele_hash(seq)
            # if str(seq_hash) == '0':
                # print("forward: "+position)
                # print(seq)
            hash_dict[gene] = seq_hash
        
        for gene, position in alleles_position[ref]['0'].items():
            position_start = int(position.split('-')[0])
            position_end = int(position.split('-')[1])
            seq = seq_record.seq[position_start:position_end].reverse_complement()
            seq_hash = allele_hash(seq)
            # if str(seq_hash) == '0':
                # print("reverse: "+position)
                # print(seq)
            hash_dict[gene] = seq_hash

with open(snakemake.output['profile'], 'w') as out:
    writer = csv.writer(out, delimiter='\t')
    n_header, n_line = ['#FILE'], [name]
    print(name)
    for h in header[1:]:
        n_header.append(h)
        if h in hash_dict:
            n_line.append(hash_dict[h])
        else:
            n_line.append('-')
    writer.writerow(n_header)
    writer.writerow(n_line)

