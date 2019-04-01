#!/usr/bin/env python3
import collections, sys
from Bio import Seq, SeqIO, SeqRecord
import os
import argparse

class PyDBGAssembler:
    def __init__(
            self, fa_or_fq_file_path_list, output_path=None, kmer_len=31, low_abund_filter_threshold=0):
        self.cwd = os.path.dirname(os.path.realpath(__file__))
        if output_path is not None:
            self.output_path = output_path
            if not os.path.exists(os.path.dirname(self.output_path)):
                os.makedirs(os.path.dirname(self.output_path))
        else:
            self.output_path = os.path.join(self.cwd, 'assembled_contigs.fasta')
        self.kmer_len = kmer_len
        self.completed_contigs_as_fasta = []
        # The minimum abundance of a kmer for it not to be discarded
        self.abund_filter_cutoff = low_abund_filter_threshold
        self.fa_or_fq_file_path_list = fa_or_fq_file_path_list
        self.contig_fwd_list = []
        self.contig_rev_list = []
        self.kmer_count_dict = collections.defaultdict(int)
        self._init_kmer_count_dict()
        # A set that holds all of the kmers that have already been used to make a contig
        self.used_kmers_set = set()
        # A list that holds all of the contigs assembled so far
        self.list_of_completed_contigs = []
        self.g_dict = {}
        # attributes that will be updated for each file or kmer
        self.reads_of_file = None
        self.kmer_list_of_contig_fwd = None
        self.kmer_list_of_contig_rev = None
        self.kmer_list_of_contig = None
        # The current kmer from which we will build a contig from in the fwd and rev_comp
        self.current_build_kmer_fwd = None
        self.current_build_kmer_rev = None

    def do_assembly(self):
        print("Starting assembly")
        for kmer_key in self.kmer_count_dict:
            if kmer_key not in self.used_kmers_set:
                self.current_build_kmer_fwd = kmer_key
                self.current_build_kmer_rev = self._rev_comp(self.current_build_kmer_fwd)
                contig_as_string, c = self._get_contig()
                for y in c:
                    self.used_kmers_set.add(y)
                    self.used_kmers_set.add(self._rev_comp(y))
                self.list_of_completed_contigs.append(contig_as_string)
        print("Done.")

        self._contig_list_to_fasta()

        self._write_out_contig_fasta()
        print("Assembly complete")

        # # TODO I don't really know what this self.g_dict business is about.
        # # I think its for creating the de Bruijn graphs as a visual output.
        # # I'll debug it through and see where it gets us, but otherwise, I can just delete
        # self.g_dict = {}
        # heads = {}
        # tails = {}
        # for i, x in enumerate(self.list_of_completed_contigs):
        #     self.g_dict[i] = ([], [])
        #     heads[x[:self.kmer_len]] = (i, '+')
        #     tails[self._rev_comp(x[-self.kmer_len:])] = (i, '-')
        #
        # for i in self.g_dict:
        #     x = self.list_of_completed_contigs[i]
        #     for y in self._fwd_seq_generator(x[-self.kmer_len:]):
        #         if y in heads:
        #             self.g_dict[i][0].append(heads[y])
        #         if y in tails:
        #             self.g_dict[i][0].append(tails[y])
        #     for z in self._fwd_seq_generator(self._rev_comp(x[:self.kmer_len])):
        #         if z in heads:
        #             self.g_dict[i][1].append(heads[z])
        #         if z in tails:
        #             self.g_dict[i][1].append(tails[z])
        #
        # return self.g_dict, self.list_of_completed_contigs

    def _write_out_contig_fasta(self):
        print(f"Writing out contig fasta to {self.output_path}")
        with open(self.output_path, 'w') as f:
            for line in self.completed_contigs_as_fasta:
                f.write(f'{line}\n')
        print("Done.")

    def _contig_list_to_fasta(self):
        for i, contig_seq in enumerate(self.list_of_completed_contigs):
            self.completed_contigs_as_fasta.extend([f'>{i}', f'{contig_seq}'])

    def _get_contig(self):
        self.kmer_list_of_contig_fwd = self._get_contig_forward(self.current_build_kmer_fwd)

        self.kmer_list_of_contig_rev = self._get_contig_forward(self.current_build_kmer_rev)

        # I think perhaps it is looking for cases where we have created a circular piece of DNA.
        # In this case there is no need to search for the reverse build as well
        if self._contig_is_circular(self.kmer_list_of_contig_fwd[-1]):
            self.kmer_list_of_contig = self.kmer_list_of_contig_fwd
        else:
            self.kmer_list_of_contig = [self._rev_comp(x) for x in self.kmer_list_of_contig_rev[-1:0:-1]] + self.kmer_list_of_contig_fwd
        return self._contig_to_string(), self.kmer_list_of_contig

    def _contig_is_circular(self, kmer):
        return self.current_build_kmer_fwd in self._fwd_seq_generator(kmer)

    def _get_contig_forward(self, kmer):
        c_fw = [kmer]

        while True:
            # If more than one of the candidate kmers exists in the dictionary
            # Abandon further contig extension
            if sum(candidate in self.kmer_count_dict for candidate in self._fwd_seq_generator(c_fw[-1])) != 1:
                break

            # If exactly one of the candidates kmers exists in the dictionary,
            # this is the candidate we will use to build forwards
            candidate = [x for x in self._fwd_seq_generator(c_fw[-1]) if x in self.kmer_count_dict][0]

            if candidate == kmer or candidate == self._rev_comp(kmer):
                break  # break out of cycles or mobius contigs

            if candidate == self._rev_comp(c_fw[-1]):
                break  # break out of hairpins

            if sum(x in self.kmer_count_dict for x in self._rev_seq_generator(candidate)) != 1:
                break

            c_fw.append(candidate)

        return c_fw

    def _init_kmer_count_dict(self):
        """Create the dictionary that will hold the kmer counts"""
        print('Creating kmer count dict')
        for file_path in self.fa_or_fq_file_path_list:
            self._get_list_of_reads_from_fq_or_fa_file(file_path)
            self._break_reads_into_kmers_and_pop_dict()
        self._remove_low_abund_kmers()
        print('Done.')

    def _break_reads_into_kmers_and_pop_dict(self):
        for read in self.reads_of_file:
            seq_str = self._get_read_sequence_as_str(read)
            list_of_sub_seqs = self._split_str_into_multiple_reads_by_n(seq_str)
            for seq in list_of_sub_seqs:
                for km in self._create_kmers_in_sequence_generator(seq):
                    self.kmer_count_dict[km] += 1
                rev_seq = self._rev_comp(seq)
                for km in self._create_kmers_in_sequence_generator(rev_seq):
                    self.kmer_count_dict[km] += 1

    @staticmethod
    def _split_str_into_multiple_reads_by_n(seq_s):
        return seq_s.split('N')

    @staticmethod
    def _get_read_sequence_as_str(read):
        return str(read.seq)

    def _remove_low_abund_kmers(self):
        low_abund_kmers = [x for x in self.kmer_count_dict if self.kmer_count_dict[x] < self.abund_filter_cutoff]
        for kmer in low_abund_kmers:
            del self.kmer_count_dict[kmer]

    def _get_list_of_reads_from_fq_or_fa_file(self, file_path):
        file_ext = os.path.splitext(file_path)[1]
        if file_ext == '.fq' or file_ext == '.fastq':
            self.reads_of_file = SeqIO.parse(file_path, 'fastq')
        elif file_ext == '.fasta' or file_ext == '.fa':
            self.reads_of_file = SeqIO.parse(file_path, 'fasta')
        else:
            raise RuntimeError(f'unknown extension of file {file_path}')

    def _create_kmers_in_sequence_generator(self, seq):
        for i in range(len(seq) - self.kmer_len + 1):
            yield seq[i:i + self.kmer_len]

    @staticmethod
    def _rev_comp(km):
        return Seq.reverse_complement(km)

    @staticmethod
    def _fwd_seq_generator(km):
        """A generator that yields the four kmers to search for (by extending in the 3' direction)
        that would allow us to build the current contig forwards
        """
        for x in 'ACGT':
            yield km[1:]+x
    @staticmethod
    def _rev_seq_generator(km):
        """A generator that yields the four kmers that the current candidate kmer could have been extended from
        """
        for x in 'ACGT':
            yield x + km[:-1]


    def _contig_to_string(self):
        return self.kmer_list_of_contig[0] + ''.join(kmer[-1] for kmer in self.kmer_list_of_contig[1:])


def process_args():
    global args
    default_output_dir = os.path.join(os.path.dirname(os.path.realpath(__file__)), 'assembled_contigs.fasta')
    parser = argparse.ArgumentParser(
        description='Simple Python de Brujin graph-based assembler adapted from https://github.com/pmelsted/dbg')
    parser.add_argument("files", nargs='*', help="The seq files to generate contigs from")
    parser.add_argument("-k", "--kmer_length", type=int, help="The length of kmer to use", default=20)
    parser.add_argument("-o", "--output_path", help="Full path to which the output fasta should be written",
                        default=default_output_dir)
    return parser.parse_args()


if __name__ == "__main__":
    args = process_args()
    dbga = PyDBGAssembler(fa_or_fq_file_path_list=args.files, kmer_len=args.kmer_length, output_path=args.output_path)
    dbga.do_assembly()
