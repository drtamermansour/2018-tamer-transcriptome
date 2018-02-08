#! /usr/bin/env python
import sys
import argparse
import screed
import khmer
from khmer.utils import write_record
import pickle
import numpy
import bbhash


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('database')
    parser.add_argument('reads')
    parser.add_argument('-k', '--ksize', type=int, default=31)
    args = parser.parse_args()

    #if args.reads == '-':
    #    args.reads = sys.stdin

    kh = khmer.Nodetable(args.ksize, 1, 1)

    mphf_filename = args.database + '.mphf'
    array_filename = args.database + '.arr'
    print('loading database {}'.format(args.database))
    
    with open(array_filename, 'rb') as fp:
        mphf_to_kmer, mphf_to_cdbg, family_ids, cdbg_to_family_id = pickle.load(fp)
    mphf = bbhash.load_mphf(mphf_filename)

    print('done!')

    def get_kmer_to_family_ids(hashval):
        mphf_hash = mphf.lookup(hashval)
        if mphf_hash is None:
            return set()
        
        kmer_hash = mphf_to_kmer[mphf_hash]
        if kmer_hash != hashval:
            return set()

        cdbg_id = mphf_to_cdbg[mphf_hash]
        id_list = cdbg_to_family_id[cdbg_id]
        return id_list

    n_same = 0
    n_fusion = 0
    n_unmatched = 0
    n_amb_same = 0
    n_amb_fusion = 0
    n_others = 0

    fusion_filename = args.database + '_fusion.fa'
    amb_fusion_filename = args.database + '_ambFusion.fa'

    fusion_fp = open(fusion_filename, 'w')
    ambfusion_fp = open(amb_fusion_filename, 'w')

    n = 0
    for record in screed.open(args.reads):
        n += 1
        if n % 1000 == 0:
            print('...', n)
            #if n > 5000:
            #    break

        hashvals = kh.get_kmer_hashes(record.sequence)
        if len(hashvals) <= 1:
            continue

        # find a matching k-mer at the beginning of the read
        first = hashvals[0]
        first_ids = get_kmer_to_family_ids(first)
        idx = 1
        while idx < len(hashvals) and len(first_ids) == 0:
            first = hashvals[idx]
            first_ids = get_kmer_to_family_ids(first)
            idx += 1

        if len(first_ids) == 0:
            print('no single match')
            n_unmatched += 1
        elif idx == len(hashvals):
            print('same, only last kmer matched ')
            n_same += 1
        else:
            # find a matching k-mer at the end of the read
            last = hashvals[-1]
            last_ids = get_kmer_to_family_ids(last)
            idy = len(hashvals) - 2
            while idy > idx and len(last_ids) == 0:
                last = hashvals[idy]
                last_ids = get_kmer_to_family_ids(last)
                idy -= 1
            
            if len(last_ids) == 0:
                print('same, only first kmer matched ')
                n_same += 1
            elif len(first_ids) == 1 and len(last_ids) == 1:
                if first_ids == last_ids: 
                    n_same += 1
                else:
                    print('fusion class A {} {} at kmers {} {}'.format(first_ids, last_ids, idx, idy+2))
                    n_fusion += 1
                    write_record(record, fusion_fp)
            else: 
                intersect_ids = first_ids.intersection(last_ids)
                if len(intersect_ids) > 0:
                    n_amb_same += 1
                else:
                    ## try to resolve ambiguity 
                    # find the least unambiguous interval at the beginning of the read
                    while idx <= idy:#  and len(first_ids) > 1:
                        first = hashvals[idx]
                        temp_ids = get_kmer_to_family_ids(first)
                        if len(temp_ids) == 0:
                            idx += 1
                        else:
                            intersect_ids = first_ids.intersection(temp_ids)
                            if len(intersect_ids) == 0:
                                break
                            else:
                                first_ids = intersect_ids
                                idx += 1

                    # find the last unambiguous interval at the end of the read
                    other = 0
                    while idy >= idx:#  and len(last_ids) > 1:
                        last = hashvals[idy]
                        temp_ids = get_kmer_to_family_ids(last)
                        if len(temp_ids) == 0:
                            idy -= 1
                        else:
                            intersect_ids = last_ids.intersection(temp_ids)
                            if len(intersect_ids) == 0:
                                other = 1
                                break
                            else:
                                last_ids = intersect_ids
                                idy -= 1

                    if other:
                        print('other {} {} {}'.format(first_ids, temp_ids, last_ids))
                        n_others += 1
                    elif len(first_ids) == 1 and len(last_ids) == 1:
                        print('fusion class B {} {} at kmers {} {}'.format(first_ids, last_ids, idx, idy+2))
                        n_fusion += 1
                        write_record(record, fusion_fp)
                    else: 
                        print('ambiguous fusion {} {} at kmers {} {}'.format(first_ids, last_ids, idx, idy+2))
                        n_amb_fusion += 1
                        ambfusion_fp
                        write_record(record, ambfusion_fp)

    print('same:', n_same)
    print('fusion:', n_fusion)
    print('unmatched:', n_unmatched)
    print('ambiguous same:', n_amb_same)
    print('ambiguous fusion:', n_amb_fusion)
    print('others:', n_others)


if __name__ == '__main__':
    main()
