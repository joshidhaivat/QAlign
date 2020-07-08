from pathlib import Path
print('Running' if __name__ == '__main__' else 'Importing', Path(__file__).resolve())

from argparse import ArgumentParser
import sys
from qalign.convert import convert
import pdb

def main():
  parser = ArgumentParser()
  subs = parser.add_subparsers()
  #pdb.set_trace()

  if sys.argv[1]=='convert':

    s_parser = subs.add_parser('convert')

    s_parser.add_argument('--input_fasta')
    s_parser.add_argument('--outdir')
    s_parser.add_argument('--qlevel', type=int, default=2, help="2 for Q2, 3 for Q3")
    s_parser.add_argument('--rc', type=int, default=1, help="1 for rev-complement, 0 to disable")
    s_parser.add_argument('--kmerpath', default="None", help="absolute path for kmer model")

    args = parser.parse_args()

    #pdb.set_trace()
    convert(
      args.input_fasta,
      args.outdir,
      args.qlevel,
      True if args.rc==1 else False,
      None if args.kmerpath=="None" else args.kmerpath
      )

  else:
    #TODO
    pass

if __name__ == "__main__":
  main()
