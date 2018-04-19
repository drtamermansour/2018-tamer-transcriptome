#! /usr/bin/env python
"""
Example usage
cat read_sample_n10000.fastq | python process_bimode.py -k 31 --force_single human_cds /dev/stdin > human_cds.process_2PE_single.log
cat pair_sample_n10000.fasta | python process_bimode.py -k 31 --paired human_cds /dev/stdin > human_cds.process_2PE_paired.log
cat pair_sample_n10000.fasta | python process_bimode.py -k 31 --paired -u read_sample_n10000.fastq human_cds /dev/stdin > human_cds.process_2PE_both.log

"""
import sys
import argparse
import screed
import khmer
from khmer.utils import (write_record, broken_paired_reader,
                         clean_input_reads)
from khmer.kfile import check_valid_file_exists
from khmer.khmer_logger import log_error
from contextlib import contextmanager
import pickle
import numpy
import bbhash

n_unmatched = 0
n_same = 0
n_amb_same = 0
n_clear_fusion = 0
n_ambig_fusion = 0
n_mutli_fusion = 0

@contextmanager
def catch_io_errors(ifile, out1, out2, out3, out4, out5, out6, force, corrupt_files):
    """Context manager to do boilerplate handling of IOErrors."""
    try:
        yield
    except (IOError, OSError, ValueError) as error:
        log_error('** ERROR: {error}', error=str(error))
        log_error('** Failed on {name}: ', name=ifile)
        if not force:
            log_error('** Exiting!')

            sys.exit(1)
        else:
            log_error('*** Skipping error file, moving on...')
            corrupt_files.append(ifile)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('database')
    parser.add_argument('input_filenames', metavar='input_sequence_filename',
                        help='Input FAST[AQ] sequence filename.', nargs='+')
    parser.add_argument('-k', '--ksize', type=int, default=31)
    parser.add_argument('-p', '--paired', action='store_true',
                        help='require that all sequences be properly paired')
    parser.add_argument('--force_single', dest='force_single',
                        action='store_true',
                        help='treat all sequences as single-ended/unpaired')
    parser.add_argument('-u', '--unpaired-reads',
                        metavar="unpaired_reads_filename",
                        help='include a file of unpaired reads to which '
                        '-p/--paired does not apply.')
    parser.add_argument('-f', '--force', dest='force',
                        help='continue past file reading errors',
                        action='store_true')
    args = parser.parse_args()

    force_single = args.force_single


    #if args.reads == '-':
    #    args.reads = sys.stdin

    # check that input files exist 
    check_valid_file_exists(args.input_filenames)

    filenames = []
    for pathfilename in args.input_filenames:
        filenames.append(pathfilename)

    # make a list of all filenames and if they're paired or not;
    # if we don't know if they're paired, default to allowing but not
    # forcing pairing.
    files = []
    for element in filenames:
        files.append([element, args.paired])
    if args.unpaired_reads:
        files.append([args.unpaired_reads, False])

    # create object of Nodetable in Khmer to use its 
    kh = khmer.Nodetable(args.ksize, 1, 1)

    # load database
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


    def readFusion(read):
        global n_unmatched,n_same,n_amb_same,n_clear_fusion,n_ambig_fusion,n_mutli_fusion
        flag = None
        lf_ids = set()
        rt_ids = set()
        families = []
        shared_kmers = []
        gaps = []

        hashvals = kh.get_kmer_hashes(read.sequence)

        # find a matching k-mer at the beginning of the read
        lf = hashvals[0]
        lf_ids = get_kmer_to_family_ids(lf)
        idx = 1
        while idx < len(hashvals) and len(lf_ids) == 0:
            lf = hashvals[idx]
            lf_ids = get_kmer_to_family_ids(lf)
            idx += 1

        if len(lf_ids) == 0:
            #print('no single match')
            n_unmatched += 1
            flag = "unmatched"
        elif idx == len(hashvals):
            #print('same, only last kmer matched')
            families.append(lf_ids)
            if len(lf_ids) == 1:
                n_same += 1 
                flag = "unique"
            else:
                n_amb_same += 1
                flag = "ambiguous"
        else: # len(lf_ids) > 0 & idx < len(hashvals)
            # find a matching k-mer at the end of the read
            rt = hashvals[-1]
            rt_ids = get_kmer_to_family_ids(rt)
            idy = len(hashvals) - 2
            while idy >= idx and len(rt_ids) == 0:
                rt = hashvals[idy]
                rt_ids = get_kmer_to_family_ids(rt)
                idy -= 1

            if len(rt_ids) == 0:
                #print('same, only one non-last kmer matched ')
                families.append(lf_ids)
                if len(lf_ids) == 1:
                    n_same += 1  
                    flag = "unique"
                else:
                    n_amb_same += 1
                    flag = "ambiguous"
            else:
                intersect_ids = lf_ids.intersection(rt_ids)
                if len(intersect_ids) > 0:
                    families.append(intersect_ids)
                    if len(intersect_ids) == 1:
                        n_same += 1
                        flag = "unique"
                    else:
                        n_amb_same += 1
                        flag = "ambiguous"
                else:# fusion to be resolved
                    shared_kmer = 1
                    gap_size = 0
                    gap = False
                    while idx <= idy+1:
                        temp = hashvals[idx]
                        temp_ids = get_kmer_to_family_ids(temp)
                        if len(temp_ids) > 0:
                            intersect_ids = lf_ids.intersection(temp_ids)
                            if len(intersect_ids) > 0:
                                lf_ids = intersect_ids
                                shared_kmer += 1
                                gap_size = 0
                            else: # len(intersect_ids) == 0
                                families.append(lf_ids)
                                shared_kmers.append(shared_kmer)
                                lf_ids = temp_ids
                                shared_kmer = 1
                                gaps.append(gap_size)
                                gap_size = 0 
                        else:
                            gap_size += 1
                        idx += 1

                    families.append(lf_ids)
                    shared_kmers.append(shared_kmer)

                    assert len(families) > 1
                    if len(families) == 2:
                        if len(families[0]) == 1 and len(families[1]) == 1:
                            n_clear_fusion += 1
                            flag = "clear_fusion"
                        else:
                            n_ambig_fusion += 1
                            flag = "ambig_fusion"
                    else: # len(families) > 2
                        n_mutli_fusion += 1
                        flag = "multi_fusion"

        #if len(families) == 0:
        #    families = "-"

        #if len(shared_kmers) == 0:
        #    shared_kmers = "-"

        return flag, families, shared_kmers, gaps

    
    fusion_filename = args.database + '_fusion.fa'
    fusion_fp = open(fusion_filename, 'w')
    fusionInfo_filename = args.database + '_fusion.info'
    fusionInfo_fp = open(fusionInfo_filename, 'w')
    print("fileName","recordIndex","whichInPair","align_class","gene_families","shared_kmers", "gaps", file=fusionInfo_fp, sep='\t')
    fusionCalc_filename = args.database + '_fusion.calc'
    fusionCalc_fp = open(fusionCalc_filename, 'w')
    print("fileName","recordIndex","whichInPair","align_class","familiy_A", "familiy_B", "no_families", "len_families", "shared_kmers", "gaps", "sorted_keys",
          file=fusionCalc_fp, sep='\t')


    fusionPairs_filename = args.database + '_fusionPairs.fa'
    fusPair_fp = open(fusionPairs_filename, 'w')
    fusionPairsInfo_filename = args.database + '_fusionPairs.info'
    fusPairInfo_fp = open(fusionPairsInfo_filename, 'w')
    print("fileName","recordIndex","fusion_class","R1_family","R2_family",file=fusPairInfo_fp, sep='\t')
    fusionPairsCalc_filename = args.database + '_fusionPairs.calc'
    fusPairCalc_fp = open(fusionPairsCalc_filename, 'w')
    print("fileName","recordIndex","fusion_class","familiy_A","familiy_B","len_families","sorted_keys", file=fusPairCalc_fp, sep='\t')

    corrupt_files = []
    family_names = dict(zip(family_ids.values(),family_ids.keys()))
    n = 0
    n_paired_fusion = 0
    sameRef = ("unique","ambiguous")
    fusion = ("clear_fusion","ambig_fusion","multi_fusion")
    for filename, require_paired in files:
        with catch_io_errors(filename, fusion_fp, fusionInfo_fp, fusionCalc_fp, fusPair_fp, fusPairInfo_fp, fusPairCalc_fp,
                             args.force, corrupt_files):
            screed_iter = clean_input_reads(screed.open(filename))
            reader = broken_paired_reader(screed_iter, min_length=args.ksize,
                                          force_single=force_single,
                                          require_paired=require_paired)

           
            for r_index, is_paired, read0, read1 in reader:
                n += 1
                if n % 10000 == 0:
                    print('...', n)
                    #if n > 5000:
                    #    break
        
                flag0, families0, shared_kmers0, gaps0 = readFusion(read0)
 
                if not is_paired  and flag0 in fusion:
                    #families_names0 = []
                    #for gp in families0:
                    #    gp_names = []
                    #    for family_id in gp: 
                    #        family_name = family_names[family_id]
                    #        gp_names.append(family_name)
                    #    families_names0.append(gp_names)

                    print(filename, r_index, "single", flag0, families0, shared_kmers0, gaps0, file=fusionInfo_fp, sep='\t')
                    write_record(read0, fusion_fp)

                    #i = 1  
                    #while i < len(families0):
                    #    for g1 in families0[i-1]:
                    #        for g2 in families0[i]:
                    #            print(filename, r_index, "single", flag0, sorted([g1,g2]), len(families0), len(families0[i-1]), len(families0[i]),
                    #                  shared_kmers0, gaps0, file=fusionCalc_fp, sep='\t')
                    #    i += 1

                    i = len(families0)-1
                    for g1 in families0[0]:
                        g1_name = family_names[g1]
                        for g2 in families0[i]:
                            g2_name = family_names[g2]
                            print(filename, r_index, "single", flag0, '{}:{}'.format(g1,g1_name), '{}:{}'.format(g2,g2_name), len(families0),
                                  [len(f) for f in families0], shared_kmers0, gaps0, sorted([g1,g2]), file=fusionCalc_fp, sep='\t')

                if is_paired:
                    flag1, families1, shared_kmers1, gaps1  = readFusion(read1)

                    if flag0 in fusion or flag1 in fusion:
                        print(filename, r_index, "Read_1", flag0, families0, shared_kmers0, gaps0, file=fusionInfo_fp, sep='\t')
                        write_record(read0, fusion_fp)
                        print(filename, r_index, "Read_2", flag1, families1, shared_kmers1, gaps1, file=fusionInfo_fp, sep='\t')
                        write_record(read1, fusion_fp)

                        if flag0 in fusion:
                            i = len(families0)-1
                            for g1 in families0[0]:
                                g1_name = family_names[g1]
                                for g2 in families0[i]:
                                    g2_name = family_names[g2]
                                    print(filename, r_index, "Read_1", flag0, '{}:{}'.format(g1,g1_name), '{}:{}'.format(g2,g2_name), len(families0),
                                          [len(f) for f in families0], shared_kmers0, gaps0, sorted([g1,g2]), file=fusionCalc_fp, sep='\t')

                        if flag1 in fusion:
                            i = len(families1)-1
                            for g1 in families1[0]:
                                g1_name = family_names[g1]
                                for g2 in families1[i]:
                                    g2_name = family_names[g2]
                                    print(filename, r_index, "Read_2", flag1, '{}:{}'.format(g1,g1_name), '{}:{}'.format(g2,g2_name), len(families1),
                                          [len(f) for f in families1], shared_kmers1, gaps1, sorted([g1,g2]), file=fusionCalc_fp, sep='\t')


                    elif flag0 in sameRef and flag1 in sameRef:
                        if len(families0[0].intersection(families1[0])) == 0:
                            n_paired_fusion += 1

                            if flag0 == "unique" and flag1 == "unique":
                                fusion_class = "clear_fusion"
                            else:
                                fusion_class = "ambig_fusion"
                       
                            print(filename, r_index, fusion_class, families0, families1, file=fusPairInfo_fp, sep='\t')
                            write_record(read0, fusPair_fp)
                            write_record(read1, fusPair_fp)

                            for g1 in families0[0]:
                                g1_name = family_names[g1]
                                for g2 in families1[0]:
                                    g2_name = family_names[g2]
                                    print(filename, r_index, fusion_class, '{}:{}'.format(g1,g1_name), '{}:{}'.format(g2,g2_name),
                                          [len(f) for f in (families0[0],families1[0])], sorted([g1,g2]), file=fusPairCalc_fp, sep='\t')


    print('No of input fragments: ', n)
    print('unmatched:', n_unmatched)
    print('Unique:', n_same)
    print('Ambiguous:', n_amb_same)
    print('Single read clear fusion:', n_clear_fusion)
    print('Single read ambiguous fusion:', n_ambig_fusion)
    print('Single read multi fusion:', n_mutli_fusion)
    print('paired read fusion:', n_paired_fusion)


if __name__ == '__main__':
    main()
