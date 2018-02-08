#! /usr/bin/env python
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
n_fusion = 0
n_amb_fusion = 0
n_mutli_fusion = 0

@contextmanager
def catch_io_errors(ifile, out1, out2, force, corrupt_files):
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


    def readFusion(read):#,whichInPair,record_index,ifile):
        global n_unmatched,n_same,n_amb_same,n_fusion,n_amb_fusion,n_mutli_fusion
        flag = -1
        first_ids = set()
        temp_ids = set()
        last_ids = set()
  
        hashvals = kh.get_kmer_hashes(read.sequence)
        #if len(hashvals) <= 1:
        #    return flag,first_ids,temp_ids,last_ids#continue

        # find a matching k-mer at the beginning of the read
        first = hashvals[0]
        first_ids = get_kmer_to_family_ids(first)
        idx = 1
        while idx < len(hashvals) and len(first_ids) == 0:
            first = hashvals[idx]
            first_ids = get_kmer_to_family_ids(first)
            idx += 1

        idy = len(hashvals) - 2
        if len(first_ids) == 0:
            print('no single match')
            n_unmatched += 1
        elif idx == len(hashvals):
            print('same, only last kmer matched')
            n_same += 1  ## Actaully this can be n_amb_same. We can differentiate by len(first_ids) 
            flag = 0
        else:
            # find a matching k-mer at the end of the read
            last = hashvals[-1]
            last_ids = get_kmer_to_family_ids(last)
            #idy = len(hashvals) - 2
            while idy > idx and len(last_ids) == 0:
                last = hashvals[idy]
                last_ids = get_kmer_to_family_ids(last)
                idy -= 1

            if len(last_ids) == 0:
                print('same, only one non-last kmer matched ')
                n_same += 1
                flag = 0
            elif len(first_ids) == 1 and len(last_ids) == 1:
                if first_ids == last_ids:
                    n_same += 1
                    flag = 0
                else:
                    #print('{} {} {} unique {} {} {} {}'.format(ifile, record_index, whichInPair,
                    #      first_ids, last_ids, idx, idy+2), file=fusionInfo_fp)
                    n_fusion += 1
                    flag = "unique"
                    #write_record(read, fusion_fp)
            else:
                intersect_ids = first_ids.intersection(last_ids)
                if len(intersect_ids) > 0:
                    n_amb_same += 1
                    flag = 0
                else:
                    ## try to resolve ambiguity 
                    # find the least unambiguous interval at the beginning of the read
                    tdx = idx
                    while tdx < idy:#  and len(first_ids) > 1:
                        first = hashvals[tdx]
                        temp_ids = get_kmer_to_family_ids(first)
                        if len(temp_ids) == 0:
                            tdx += 1  ## skip unmatched kmers
                        else:
                            intersect_ids = first_ids.intersection(temp_ids)
                            if len(intersect_ids) != 0:
                                first_ids = intersect_ids ## select the shared family ids
                                idx += 1
                                temp_ids = set()
                                tdx += 1
                            else:
                                break    ## stop if you reach a likely fusion point

                    # find the least unambiguous interval at the end of the read
                    multi = 0
                    tdy = idy
                    while tdy >= tdx:#  and len(last_ids) > 1:
                        last = hashvals[tdy]
                        temp_ids = get_kmer_to_family_ids(last)
                        if len(temp_ids) == 0:
                            tdy -= 1  ## skip unmatched kmers
                        else:
                            intersect_ids = last_ids.intersection(temp_ids)
                            if len(intersect_ids) != 0:
                                last_ids = intersect_ids ## select the shared family ids
                                idy -= 1
                                temp_ids = set()
                                tdy -= 1
                            else:
                                multi = 1
                                break   ## stop if you reach a likely second fusion point in the same read

                    if multi:
                        #print('{} {} {} multi {} {} {} {} {} {}'.format(ifile, record_index, whichInPair,
                        #      first_ids, temp_ids, last_ids, idx, idy+1, idy+2), file=multifusInfo_fp)
                        n_mutli_fusion += 1
                        flag = "multi"
                        #write_record(read, multifusInfo_fp)
                    elif len(first_ids) == 1 and len(last_ids) == 1:
                        #print('{} {} {} resolved {} {} {} {}'.format(ifile, record_index, whichInPair,
                        #      first_ids, last_ids, idx, idy+2), file=fusionInfo_fp)
                        n_fusion += 1
                        flag = "resolved"
                        #write_record(read, fusion_fp)
                    else:
                        #print('{} {} {} ambiguous {} {} {} {}'.format(ifile, record_index, whichInPair,
                        #      first_ids, last_ids, idx, idy+2), file=fusionInfo_fp)
                        n_amb_fusion += 1
                        flag = "ambiguous"
                        #write_record(read, ambfusion_fp)

        if len(last_ids) == 0:
            last_ids = first_ids

        if len(temp_ids) == 0:
            temp_ids = "-"

        return flag, first_ids, temp_ids, last_ids, idx, idy+2


    fusion_filename = args.database + '_fusion.fa'
    fusion_fp = open(fusion_filename, 'w')
    fusionInfo_filename = args.database + '_fusion.info'
    fusionInfo_fp = open(fusionInfo_filename, 'w')
    print("fileName","recordIndex","whichInPair","fusion_class","family_A","family_mid","family_B","kmer_pos_A","kmer_pos_B",
          file=fusionInfo_fp, sep='\t')

    #ambfusion_filename = args.database + '_ambFusion.fa'
    #ambfusion_fp = open(ambfusion_filename, 'w')
    #ambfusionInfo_filename = args.database + '_ambFusion.info'
    #ambfusionInfo_fp = open(ambfusionInfo_filename, 'w')
    #print('fileName recordIndex whichInPair fusion_class family_A family_B kmer_pos_A kmer_pos_B',
    #      file=ambfusionInfo_fp)

    #multifus_filename = args.database + '_multiFusion.fa'
    #multifus_fp = open(multifus_filename, 'w')
    #multifusInfo_filename = args.database + '_multiFusion.info'
    #multifusInfo_fp = open(multifusInfo_filename, 'w')
    #print('fileName recordIndex whichInPair fusion_class family_A family_B family_C kmer_pos_A kmer_pos_B kmer_pos_C',
    #      file=multifusInfo_fp)

    fusionPairs_filename = args.database + '_fusionPairs.fa'
    fusPair_fp = open(fusionPairs_filename, 'w')
    fusionPairsInfo_filename = args.database + '_fusionPairs.info'
    fusPairInfo_fp = open(fusionPairsInfo_filename, 'w')
    print("fileName","recordIndex","fusion_class","R1_family","R2_family","R1_kmer_pos","R2_kmer_pos",
          file=fusPairInfo_fp, sep='\t')

    corrupt_files = []
    n = 0
    n_paired_fusion = 0
    for filename, require_paired in files:
        with catch_io_errors(filename, fusion_fp, fusPair_fp,
                             args.force, corrupt_files):
            screed_iter = clean_input_reads(screed.open(filename))
            reader = broken_paired_reader(screed_iter, min_length=args.ksize,
                                          force_single=force_single,
                                          require_paired=require_paired)

           
            for r_index, is_paired, read0, read1 in reader:
                n += 1
                if n % 1000 == 0:
                    print('...', n)
                    #if n > 5000:
                    #    break
        
                flag0, first_ids0, temp_ids0, last_ids0, Pos_A0, Pos_B0 = readFusion(read0)#,1,r_index,filename)
                if not is_paired  and flag0 not in (0,-1):
                    print(filename, r_index, "single", flag0,
                          first_ids0, temp_ids0, last_ids0, Pos_A0, Pos_B0, file=fusionInfo_fp, sep='\t')
                    write_record(read0, fusion_fp)

                if is_paired:
                    flag1, first_ids1, temp_ids1, last_ids1, Pos_A1, Pos_B1 = readFusion(read1)#,2,r_index,filename)
                    
                    if flag0 not in (0,-1) or flag1 not in (0,-1):
                        print(filename, r_index, "Read_1", flag0,
                              first_ids0, temp_ids0, last_ids0, Pos_A0, Pos_B0, file=fusionInfo_fp, sep='\t')
                        write_record(read0, fusion_fp)
                        print(filename, r_index, "Read_2", flag1,
                              first_ids1, temp_ids1, last_ids1, Pos_A1, Pos_B1, file=fusionInfo_fp, sep='\t')
                        write_record(read1, fusion_fp)

                    elif flag0 == 0 and flag1 == 0:
                        R1_family = first_ids0.intersection(last_ids0)
                        R2_family = first_ids1.intersection(last_ids1)
                        if len(R1_family.intersection(R2_family)) == 0:
                            n_paired_fusion += 1
                            if len(R1_family) == 1 and len(R2_family) == 1:
                                fusion_class = "unique"
                            else:
                                fusion_class = "ambiguous"
                       
                            print(filename, r_index, fusion_class, R1_family, R2_family,
                                  "{},{}".format(Pos_A0, Pos_B0), "{},{}".format(Pos_A1, Pos_B1), file=fusPairInfo_fp, sep='\t')
                            write_record(read0, fusPair_fp)
                            write_record(read1, fusPair_fp)



    print('same:', n_same)
    print('fusion:', n_fusion)
    print('unmatched:', n_unmatched)
    print('ambiguous same:', n_amb_same)
    print('ambiguous fusion:', n_amb_fusion)
    print('multi fusion:', n_mutli_fusion)
    print('paired fusion:', n_paired_fusion)


if __name__ == '__main__':
    main()
