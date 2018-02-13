#! /usr/bin/env python
# https://www.daniweb.com/programming/software-development/code/216636/multiple-word-replace-in-text-python
# https://stackoverflow.com/questions/15658187/replace-all-words-from-word-list-with-another-string-in-python
import sys
import argparse
import screed
import re

def multiwordReplace(text, wordDic):
    """
    take a text and replace words that match a key in a dictionary with
    the associated value, return the changed text
    """
    rc = re.compile(r'\b%s\b' % r'\b|\b'.join(map(re.escape, wordDic)))
    #rc = re.compile('|'.join(map(re.escape, wordDic)))
    def translate(match):
        return wordDic[match.group(0)]
    return rc.sub(translate, text)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('transcriptome', help='Input reference transcriptome with headers having ens gene and common names in 2nd and 6th positions.')
    parser.add_argument('input_filenames', help='input file to update', metavar='Text file with Ensembl IDs', nargs='+')
    parser.add_argument('-o', '--output')
    args = parser.parse_args()
    assert args.output

    family_ids = {}
    n = 0
    for record in screed.open(args.transcriptome):
        n += 1
        if n % 1000 == 0:
            print('...', args.transcriptome, n)
            #if n > 5000 and 0:
            #    break

        # fill the dectionary
        ens_name = record.name.split('|')[1]
        com_name = record.name.split('|')[5]
        #family_ids[ens_name] = com_name
        family_ids[ens_name] = '{}|{}'.format(ens_name,com_name)


    for in_filename in args.input_filenames:
        with open(in_filename, 'r') as inp, open(args.output, 'w') as outp:
            old_text=inp.read()
            new_text = multiwordReplace(old_text, family_ids)
            print(new_text.rstrip('\n'), file=outp) 


if __name__ == '__main__':
    main()

