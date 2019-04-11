from collections import defaultdict, Counter
import pickle
from os.path import join, exists, splitext
import time
import os
import sys
import argparse
sys.path.insert(0, os.path.abspath(".."))
sys.path.insert(0, os.path.abspath("../.."))
from CM122_starter_code.helpers import read_reads, read_reference, pretty_print_aligned_reads_with_ref


def hash_end(end, genome_ht):
    """
    Uses hashing to identify the set of locations spanned by
    a read.

    :param end: A single end of a read
    :param genome_ht: A hash of the genome with uniform key length
    :return:
    """
    key_length = max(map(len, genome_ht))
    end_pieces = [end[i * key_length: (i + 1) * key_length]
                  for i in range(len(end) // key_length)]

    hashed_read_locations = [genome_ht[read_piece]
                             for read_piece in end_pieces]
    start_positions = [[x - i * key_length for x in hashed_read_locations[i]]
                       for i in range(len(hashed_read_locations))]
    start_counter = Counter()

    for position_list in start_positions:
        start_counter.update(position_list)

    if not start_counter:
        return -1, 0
    else:
        best_alignment_location, best_alignment_count = \
            start_counter.most_common(1)[0]

    if best_alignment_count < 2:
        return -1, best_alignment_count
    else:
        return best_alignment_location, best_alignment_count


def hash_read(read, genome_ht):
    """
    Uses hashing to identify the set of locations spanned by
    a read.

    :param read: A single read
    :param genome_ht: A hash of the genome with uniform key length
    :return:
    """

    oriented_reads = [(read[0][::i], read[1][::j]) for i, j in ((1, -1), (-1, 1))]
    ## Either one end is forward and the other end is reversed, or vice versa.

    best_score = -1
    best_alignment_locations = (-1, -1)
    best_oriented_read = ('', '')
    for oriented_read in oriented_reads:
        hash_results = [hash_end(_, genome_ht) for _ in oriented_read]
        hash_locations = [_[0] for _ in hash_results]
        hash_score = sum([_[1] for _ in hash_results])
        if hash_score > best_score:
            best_alignment_locations = hash_locations
            best_oriented_read = oriented_read
    return best_oriented_read, best_alignment_locations


def make_genome_hash(reference, key_length):
    """

    :param reference: The reference as a string stored
    :param key_length: The length of keys to use.
    :return:
    """
    genome_hash = defaultdict(list)
    for i in range(len(reference) - key_length):
        ref_piece = reference[i: i + key_length]
        genome_hash[ref_piece].append(i)
    return genome_hash


def build_hash_and_pickle(ref_fn, key_length, force_rebuild=False):
    reference_hash_pkl_fn = '{}_hash_keylength_{}.pkl'.format(splitext(ref_fn)[0], key_length)
    if exists(reference_hash_pkl_fn) and not force_rebuild:
        ref_genome_hash = pickle.load(open(reference_hash_pkl_fn, 'rb'))
        if max(map(len, ref_genome_hash)) == key_length:
            return ref_genome_hash
        else:
            pass
    else:
        pass
    reference = read_reference(ref_fn)
    ref_genome_hash = make_genome_hash(reference, key_length)
    pickle.dump(ref_genome_hash, open(reference_hash_pkl_fn, 'wb'))
    return ref_genome_hash


def hashing_algorithm(paired_end_reads, genome_ht):
    """

    :param paired_end_reads:
    :param genome_ht:
    :return:
    """
    alignments = []
    genome_aligned_reads = []
    count = 0

    for read in paired_end_reads:
        alignment, genome_aligned_read = hash_read(read, genome_ht)
        alignments.append(alignment)
        genome_aligned_reads.append(genome_aligned_read)
        count += 1
        if count % 100 == 0:
            time_passed = (time.process_time())/60
            print('{} reads aligned'.format(count), 'in {:.3} minutes'.format(time_passed))
            remaining_time = time_passed/count*(len(paired_end_reads)-count)
            print('Approximately {:.3} minutes remaining'.format(remaining_time))
    return alignments, genome_aligned_reads

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='basic_hasher.py takes in data for homework assignment 2 consisting '
                                     'of a genome and a set of reads and aligns the reads to the reference genome.')
    parser.add_argument('-g', '--referenceGenome', required=True, dest='reference_file',
                        help='File containing a reference genome.')
    parser.add_argument('-r', '--reads', required=True, dest='reads_file',
                        help='File containg sequencing reads.')
    parser.add_argument('-o', '--outputFile', required=True, dest='output_file',
                        help='Output file name.')

    args = parser.parse_args()
    reference_fn = args.reference_file
    reads_fn = args.reads_file
    output_fn = args.output_file

    key_length = 7
    reads = read_reads(reads_fn)
    # If you want to speed it up, cut down the number of reads by
    # changing the line to reads = read_reads(reads_fn)[:<x>] where <x>
    # is the number of reads you want to work with.
    genome_hash_table = build_hash_and_pickle(reference_fn, key_length)
    ref = read_reference(reference_fn)
    genome_aligned_reads, alignments = hashing_algorithm(reads, genome_hash_table)
    output_str = pretty_print_aligned_reads_with_ref(genome_aligned_reads, alignments, ref)
    with(open(output_fn, 'w')) as output_file:
        output_file.write(output_str)
