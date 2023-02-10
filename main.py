import numpy as np
from Bio.SeqIO import parse
from src.motif_extraction import extract_motifs
import os
from argparse import ArgumentParser

parser = ArgumentParser()

parser.add_argument('-contexts', help='fna to enrich', type=str, required=True)
parser.add_argument('-ref', help='ref',type=str, required=True)
parser.add_argument('-max_motifs', help='max motifs', default=20, type=int)
parser.add_argument('-min_conf', help='min conf', default=1000)
parser.add_argument('-savepath', default=None, type=str)

args = parser.parse_args()


f = parse(args.contexts, format='fasta')
seqs = [str(rec.seq) for rec in f]

for rec in parse(args.ref, format='fasta'):
    reference = str(rec.seq)
    break



if args.savepath is None:
    import datetime
    
    sp = str(datetime.datetime.now()
        ).replace(' ', '_').replace(':', '').replace('-', '_').split('.')[0]
    outdir = 'Results_' + sp

else:
    outdir  = args.savepath

try:
    os.mkdir(outdir)
except:
    raise FileExistsError('The specified output dir already exists!')


DETAILED_MOTIF_SET = extract_motifs(
    seqs, 
    reference, 
    outdir, 
    args.max_motifs,
    args.min_conf, 
    threads=60
)

print('\n---------RESULTS---------\n')
for d in DETAILED_MOTIF_SET:
    print(d)


with open(f'{args.savepath}/results.txt', 'w') as f:
    for d in DETAILED_MOTIF_SET:
        f.write(f'{d}\n')

    