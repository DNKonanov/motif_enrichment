from random import sample, random
from Bio.SeqIO import parse
from argparse import ArgumentParser
import tqdm

parser = ArgumentParser()

parser.add_argument('-genome', required=True, help='reference genome', type=str)
parser.add_argument('-outdir', required=True, help='outdir', type=str)
parser.add_argument('-err_rate', default=0.01, help='False Positive Rate (default 0.01)', type=float)
parser.add_argument('-replicates', default=10, help='the number of replicates to generate (default 10)', type=int)

args = parser.parse_args()

f = parse(args.genome, format='fasta')

for rec in f:
    reference = str(rec.seq)
    break



from random import sample, random

test_size = 5
test_lens = (4,5,6)


od = args.outdir

import os
try:
    os.mkdir(args.outdir)
except FileExistsError:
    raise FileExistsError

for replicate in range(args.replicates):
    print(f'Replicate {replicate} processing...')
    target = []
    
    for l in test_lens:

        letters = set(['A','G','T','C'])

        motifs = list(set([reference[i:i+l] for i in range(len(reference)-l)]))    
   
        motif_freqs = {m:0 for m in motifs}

        for i in range(len(reference)-l):
            motif_freqs[reference[i:i+l]] += 1

        vector = sorted(
            [(motif_freqs[m], m) for m in motif_freqs if ('A' in m or 'C' in m)]
        )   


        target += sample(vector[:50], k=2)
        target += sample(vector[50:-50], k=2)

    
    with open(f'{od}/target_motifs_{replicate}.txt', 'w') as f:

        for t in target:
            f.write(f'{t}\n')
            
    target = [t[1] for t in target]    
    
    
    contexts = []
    for pos in tqdm.tqdm(range(len(reference)-11)):

        add = False
        for m in target:
            if m in reference[pos+1:pos+9]:
                contexts.append(reference[pos:pos+11])
                add = True
                break


        if add == False:
            if random() < args.err_rate:
                contexts.append(reference[pos:pos+11])
                
    contexts = set(contexts)
    print(f'Number of contexts in replicate {replicate}: {len(contexts)}')
    print()
    
    with open(f'{od}/test_{replicate}.fna', 'w') as f:


        cnt = 1

        for c in contexts:
            f.write(f'>context_{cnt}\n')
            f.write(f'{c}\n')
            cnt += 1