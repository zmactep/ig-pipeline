import os
import argparse
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


def main():
    parser = argparse.ArgumentParser(description='Cuts the region from FASTA based on KABAT and region number')
    parser.add_argument('--input', nargs=1, help='input fasta')
    parser.add_argument('--kabat', nargs=1, help='input kabat')
    parser.add_argument('--outdir', nargs=1, help='output dir')
    parser.add_argument('--reg_num', nargs=1, help='0-based region number')
    args, unknown = parser.parse_known_args()

    cut(args.input.pop(0), args.outdir.pop(0), args.kabat.pop(0), int(args.reg_num.pop(0)))


def cut(input, outdir, kabat, reg_num):
    with open(kabat, 'rU') as kabat_file:
        ref_dict = dict([parse_kabat_line(line) for line in kabat_file])

    out = open(os.path.join(outdir, 'output.fasta'), 'w')
    for record in SeqIO.parse(input, "fasta"):
        if record.id in ref_dict:
            regions_list = ref_dict[record.id]
            if len(regions_list) > 2 * reg_num:
                start_index = int(regions_list[2 * reg_num]) - 1
                stop_index = int(regions_list[2 * reg_num + 1])
                sub_seq = Seq(str(record.seq)[start_index: stop_index])
                cut_record = SeqRecord(sub_seq, id=record.id, name=record.name, description=record.description)
                SeqIO.write(cut_record, out, "fasta")
            else:
                print('Region %s is out of bound for %s' % (reg_num, record.id))
        else:
            print('%s not found in kabat' % record.id)


def parse_kabat_line(line):
    """ return: key = sequence_name, value = list of regions """
    tokens = line.split()
    key = tokens.pop(0).strip()
    return key, tokens


if __name__ == "__main__":
    main()

