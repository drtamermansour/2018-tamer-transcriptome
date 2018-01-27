#! /usr/bin/env python
"""
construct cdbg like so:

bcalm -in gencode.v27.pc_transcripts.fa.gz -kmer-size 31 -abundance-min 1 -out bcalm.gencode

"""
import sys
import argparse
import screed
import khmer
import pickle
import shelve
from collections import defaultdict
import bbhash
import numpy


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('unitigs')
    parser.add_argument('transcriptomes', nargs='+')
    parser.add_argument('-k', '--ksize', type=int, default=31)
    parser.add_argument('-o', '--output')
    args = parser.parse_args()

    assert args.output

    kh = khmer.Nodetable(args.ksize, 1, 1)

    ## store all the kmers of the compacted DBG into bbhash index
    all_kmers = []
    for n, record in enumerate(screed.open(args.unitigs)):
        if n % 10000 == 0:
            print('... cdbg', n)
            if n > 20000 and 0:
                break

        all_kmers.extend(kh.get_kmer_hashes(record.sequence))

    print('building MPHF for {} k-mers in {} nodes.'.format(len(all_kmers), n))
    x = bbhash.PyMPHF(all_kmers, len(all_kmers), 4, 1.0)

    ###

    ## create arrays of kmers and nodes of cdbg using the mphf as an array index
    mphf_to_kmer = numpy.zeros(len(all_kmers), numpy.uint64)
    mphf_to_cdbg = numpy.zeros(len(all_kmers), numpy.uint32)

    for n, record in enumerate(screed.open(args.unitigs)):
        if n % 10000 == 0:
            print('... cdbg', n)
            if n > 20000 and 0:
                break

        cdbg_id = int(record.name.split(' ')[0])
        kmers = kh.get_kmer_hashes(record.sequence)

        for kmer in kmers:
            mphf = x.lookup(kmer)
            mphf_to_kmer[mphf] = kmer
            mphf_to_cdbg[mphf] = cdbg_id
        
    ###

    print('walking the transcriptome')

    family_ids = {}
    family_counter = 0

    cdbg_to_family_id = defaultdict(set)

    n = 0
    for tr_filename in args.transcriptomes:
        for record in screed.open(tr_filename):
            n += 1
            if n % 1000 == 0:
                print('...', tr_filename, n)
                if n > 5000 and 0:
                    break
            
            if len(record.sequence) < args.ksize:
                    break
   
            # get the family name
            family_name = record.name.split('|')[1]

            # convert to family ID, generating a new one if we need one
            family_id = family_ids.get(family_name)
            if family_id is None:
                family_id = family_counter
                family_counter += 1
                family_ids[family_name] = family_id

            # for all k-mers, 
            hashvals = kh.get_kmer_hashes(record.sequence)
            for hashval in hashvals:

                # find cDBG ID
                mphf = x.lookup(hashval)
                #if mphf is None:
                #    continue
                assert mphf is not None
                cdbg_id = mphf_to_cdbg[mphf]

                # link cDBG ID to family ID
                cdbg_to_family_id[cdbg_id].add(family_id)

    mphf_filename = args.output + '.mphf'
    array_filename = args.output + '.arr'
    x.save(mphf_filename)

    with open(array_filename, 'wb') as fp:
        pickle.dump((mphf_to_kmer, mphf_to_cdbg, family_ids, cdbg_to_family_id),
                    fp)


if __name__ == '__main__':
    main()
