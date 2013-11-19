from http.cookiejar import request_host
import os
import argparse
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from parsers import kabat


def main():
    parser = argparse.ArgumentParser(description='Cuts the region from FASTA based on KABAT and region number')
    parser.add_argument('--input', nargs=1, required=True, help='input fasta')
    parser.add_argument('--kabat', nargs=1, required=True, help='input kabat')
    parser.add_argument('--outdir', nargs=1, required=True, help='output dir')
    parser.add_argument('--reg_num', nargs=1, required=True, help='0-based region number')
    args, unknown = parser.parse_known_args()

    cut(args.input.pop(0), args.outdir.pop(0), args.kabat.pop(0), int(args.reg_num.pop(0)))


def cut(input, outdir, kabat_filename, reg_num):
    with open(kabat_filename, 'rU') as kabat_file:
        ref_dict = kabat.parse(kabat_file)

    out = open(os.path.join(outdir, 'output.fasta'), 'w')
    for record in SeqIO.parse(input, "fasta"):
        if record.id in ref_dict:
            regions_list = ref_dict[record.id]
            if len(regions_list) > reg_num:
                start_index, stop_index = regions_list[reg_num]  # indices still 1-based!
                sub_seq = record.seq[start_index - 1: stop_index]
                cut_record = SeqRecord(sub_seq, id=record.id, name=record.name, description=record.description)
                SeqIO.write(cut_record, out, "fasta")
            else:
                print('Region %s is out of bound for %s' % (reg_num, record.id))
        else:
            print('%s not found in kabat' % record.id)


if __name__ == "__main__":
    main()

