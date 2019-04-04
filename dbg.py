#!/usr/bin/env python3
import collections, sys
from Bio import Seq, SeqIO, SeqRecord
import os
import argparse
import sys

class PyDBGAssembler:
    def __init__(
            self, fa_or_fq_file_path_list, kmer_len, is_relaxed, is_pick_random, verbose, remove,
            low_abund_filter_threshold, assess_directionality, output_path=None,
            ):
        self.cwd = os.path.dirname(os.path.realpath(__file__))
        if output_path is not None:
            self.output_path = output_path
            if not os.path.exists(os.path.dirname(self.output_path)):
                os.makedirs(os.path.dirname(self.output_path))
        else:
            self.output_path = os.path.join(self.cwd, 'assembled_contigs.fasta')
        self.is_relaxed = is_relaxed
        if kmer_len is None:
            if self.is_relaxed:
                self.kmer_len = 20
            else:
                self.kmer_len = 31
        else:
            self.kmer_len = kmer_len
        self.completed_contigs_as_fasta = []
        # The minimum abundance of a kmer for it not to be discarded
        self.abund_filter_cutoff = int(low_abund_filter_threshold)
        self.fa_or_fq_file_path_list = fa_or_fq_file_path_list
        self.contig_fwd_list = []
        self.contig_rev_list = []
        # the total number of kmers e.g. [AACT, AACT, AACT, AGAT] == 4 not 2.
        self.total_number_of_kmer_seqs = 0
        self.kmer_count_dict = collections.defaultdict(int)
        self._init_kmer_count_dict()
        if not is_pick_random:
            self.kmer_list_for_contigs = [a[0] for a in
                                           sorted(self.kmer_count_dict.items(), key=lambda x:x[1], reverse=True)]
        else:
            self.kmer_list_for_contigs = list(self.kmer_count_dict.keys())
        self.num_kmers = len(self.kmer_list_for_contigs)
        # A set that holds all of the kmers that have already been used to seed a contig extension
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
        self.verbose = verbose
        self.remove = remove
        self.assess_directionality = assess_directionality
        # Attributes used in directionality assessment
        self.lhs_seed_kmer = None
        self.rhs_seed_kmer = None
        self.forward_contig_list = None
        self.reverse_contig_list = None
        # The set of kmers that have been used in building the current contig
        self.temp_used_build_kmers_set = set()



    def do_assembly(self):
        print("Starting assembly...")
        sys.stdout.write('.')
        count = 0
        one_percent = int(self.num_kmers / 100)
        for kmer_key in self.kmer_list_for_contigs:
            if count%one_percent == 0:
                sys.stdout.write('.')
                sys.stdout.flush()
            count += 1
            if kmer_key not in self.used_kmers_set:
                self.current_build_kmer_fwd = kmer_key
                self.current_build_kmer_rev = self._rev_comp(self.current_build_kmer_fwd)
                contig_as_string, contig_kmer_list = self._get_contig()
                for used_kmer in contig_kmer_list:
                    self.used_kmers_set.add(used_kmer)
                    self.used_kmers_set.add(self._rev_comp(used_kmer))
                self.list_of_completed_contigs.append(contig_as_string)
        print("\nDone.")

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
        """This is where we go from a seed kmer and return a contig as both a string and as a list of kmers.
        If we are not assessing directionality (--assess_directionality) then this takes the form of extending
        the kmer forwards and then back and then returning it. If we are assessing directionality then we extend
        the kmer forward and back as before but then we take the first and last kmer of the assembled contig
        and do a new assembly in the forward and rever direction using these two kmers as seeds for each respectively.
        We then check to see if the two contigs created agree. If they do, great, if the don't then we return
        (i.e. choose) the contig that has the highest mean support of its kmers. If we are working with the --remmove
        flag, this means taht we are strictly only allowing each kmer to be use once across the entire assembly
        of all contigs. Rather than directly delete the used kmers when in remove mode we put them into a used
        list and see if kmers are in the used list then at the end of get contig creation we remove the kmers
        from the dict. This way we can easily reinstate the set of kmers available when doing directionality
        assessment"""

        self.kmer_list_of_contig_fwd = self._get_contig_forward(self.current_build_kmer_fwd)

        self.kmer_list_of_contig_rev = self._get_contig_forward(self.current_build_kmer_rev)

        # I think perhaps it is looking for cases where we have created a circular piece of DNA.
        # In this case there is no need to search for the reverse build as well
        if self._contig_is_circular(self.kmer_list_of_contig_fwd[-1]):
            self.kmer_list_of_contig = self.kmer_list_of_contig_fwd
        else:
            self.kmer_list_of_contig = [self._rev_comp(x) for x in self.kmer_list_of_contig_rev[-1:0:-1]] + self.kmer_list_of_contig_fwd
            if self.assess_directionality and self.kmer_list_for_contigs:
                return self._assess_directionality()
        if self.remove:  # delete the kmers in the temp used set from the dict before exiting
            self._del_temp_kmers_from_count_dict()
        return self._contig_to_string(), self.kmer_list_of_contig

    def _assess_directionality(self):
        self._make_forward_extension_contig()
        self._make_reverse_extension_contig()
        if not self.reverse_contig_list or not self.forward_contig_list:
            if self.reverse_contig_list:
                self._del_temp_kmers_from_count_dict(self.reverse_contig_list)
                return self._contig_to_string(self.reverse_contig_list), self.reverse_contig_list
            elif self.forward_contig_list:
                self._del_temp_kmers_from_count_dict(self.forward_contig_list)
                return self._contig_to_string(self.forward_contig_list), self.forward_contig_list
            else:  # neither list is good
                pass


        return self._if_agree_return_fwd_contig_else_return_contig_with_most_support()

    def _if_agree_return_fwd_contig_else_return_contig_with_most_support(self):
        if self.forward_contig_list == self.reverse_contig_list:  # If no disagreement then contig is good
            if self.remove:
                self._del_temp_kmers_from_count_dict()
            return self._contig_to_string(self.forward_contig_list), self.forward_contig_list
        else:  # contigs disagree then take the contig that has the highest average support for its kmers
            return self._return_contig_with_highest_average_kmer_support()

    def _make_reverse_extension_contig(self):
        self.rhs_seed_kmer = self._rev_comp(self.kmer_list_of_contig[-1])
        # reset the used kmer set
        self._reset_temp_used_kmer_set()
        self.reverse_contig_list = [self._rev_comp(x) for x in self._get_contig_forward(self.rhs_seed_kmer)[-1:0:-1]]

    def _make_forward_extension_contig(self):

        self.lhs_seed_kmer = self.kmer_list_of_contig[0]
        # reset the used kmer set
        self._reset_temp_used_kmer_set()
        self.forward_contig_list = self._get_contig_forward(self.lhs_seed_kmer)

    def _return_contig_with_highest_average_kmer_support(self):
        forward_mean_kmer_support = self._calc_mean_kmer_support(self.forward_contig_list)
        reverse_mean_kmer_support = self._calc_mean_kmer_support(self.reverse_contig_list)
        if forward_mean_kmer_support > reverse_mean_kmer_support:
            if self.remove:
                self._del_temp_kmers_from_count_dict(self.forward_contig_list)
            return self._contig_to_string(self.forward_contig_list), self.forward_contig_list
        else:
            if self.remove:
                self._del_temp_kmers_from_count_dict(self.reverse_contig_list)
            return self._contig_to_string(self.reverse_contig_list), self.reverse_contig_list

    def _reset_temp_used_kmer_set(self):
        self.temp_used_build_kmers_set = set()

    def _calc_mean_kmer_support(self, kmer_list):
        return sum([self.kmer_count_dict[kmer] for kmer in kmer_list])/len(kmer_list)

    def _del_temp_kmers_from_count_dict(self, kmer_set_of_list=None):
        """ When we are working with --remove flag then we are strictly only allowing each kmer to be used once
        in building the kmers. We keep track of which kmers have been used through the contig creations using
        the self.temp_used_build_kmers_set. At then end of the contig building we need to make sure that these kmers
        are deleted from the dict and the kmer set is reset."""

        if kmer_set_of_list is None:
            for kmer in self.temp_used_build_kmers_set:
                del self.kmer_count_dict[kmer]
                if self._rev_comp(kmer) != kmer: # cannot delete a seq if rev comp is same as fwd comp.
                    del self.kmer_count_dict[self._rev_comp(kmer)]
        else:
            for kmer in kmer_set_of_list:
                del self.kmer_count_dict[kmer]
                if self._rev_comp(kmer) != kmer: # cannot delete a seq if rev comp is same as fwd comp.
                    del self.kmer_count_dict[self._rev_comp(kmer)]
        self.temp_used_build_kmers_set = set()


    def _contig_is_circular(self, kmer):
        return self.current_build_kmer_fwd in self._fwd_seq_generator(kmer)

    def _get_contig_forward(self, kmer):
        c_fw = [kmer]

        while True:
            # cand is abbrev. for candidate
            cand_kmers = [(cand, self.kmer_count_dict[cand]) for cand in
                          self._fwd_seq_generator(c_fw[-1]) if
                          cand in self.kmer_count_dict if
                          cand not in self.temp_used_build_kmers_set]

            if not self.is_relaxed: # operating in strict mode
                # If more than one of the candidate kmers exists in the dictionary
                # Abandon further contig extension
                if len(cand_kmers) != 1:
                    break
                    # If exactly one of the candidates kmers exists in the dictionary,
                    # this is the candidate we will use to build forwards
                candidate = cand_kmers[0][0]
            else: # operating in relaxed mode
                # allow there to be more than one suitable candidate and work with the candidate that is more abundant
                if len(cand_kmers) != 0:
                    candidate = sorted(cand_kmers, key=lambda x:x[1], reverse=True)[0][0]
                else: # no candidates
                    break

            if candidate == kmer or candidate == self._rev_comp(kmer):
                break  # break out of cycles or mobius contigs

            if candidate == self._rev_comp(c_fw[-1]):
                break  # break out of hairpins

            if not self.is_relaxed:
                # This check will only be made if we are operating in --strict mode
                if sum(x in self.kmer_count_dict for x in self._rev_seq_generator(candidate)) != 1:
                    break
            else:
                if candidate == c_fw[-1]:
                    if self.verbose:
                        print(f'breaking contig assembly due to homopolymer {candidate}')
                    break

            c_fw.append(candidate)
            if self.remove:
                self.temp_used_build_kmers_set.add(candidate)
                self.temp_used_build_kmers_set.add(self._rev_comp(candidate))

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
                    self.total_number_of_kmer_seqs += 1
                rev_seq = self._rev_comp(seq)
                for km in self._create_kmers_in_sequence_generator(rev_seq):
                    self.kmer_count_dict[km] += 1
                    self.total_number_of_kmer_seqs += 1
        print(f'{self.total_number_of_kmer_seqs} total kmers found; represented as '
              f'{len(self.kmer_count_dict)} unique kmers')

    @staticmethod
    def _split_str_into_multiple_reads_by_n(seq_s):
        return seq_s.split('N')

    @staticmethod
    def _get_read_sequence_as_str(read):
        return str(read.seq)

    def _remove_low_abund_kmers(self):
        low_abund_kmers = [x for x in self.kmer_count_dict if
                           self.kmer_count_dict[x] < self.abund_filter_cutoff]
        print(f'{len(low_abund_kmers)} unique kmers to be removed at cutoff threshold of {self.abund_filter_cutoff}')
        print(f'removing low abund kmers from count dict...')
        for kmer in low_abund_kmers:
            del self.kmer_count_dict[kmer]

        print(f'{sum(self.kmer_count_dict.values())} total kmers and {len(self.kmer_count_dict)} unique kmers '
              f'now remaining in count dict')
        print('Done.')
        print('Count dictionary population complete.')

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


    def _contig_to_string(self, contig_list=None):
        if contig_list is None:
            return self.kmer_list_of_contig[0] + ''.join(kmer[-1] for kmer in self.kmer_list_of_contig[1:])
        else:
            return contig_list[0] + ''.join(kmer[-1] for kmer in contig_list[1:])


def process_args():
    global args
    default_output_dir = os.path.join(os.path.dirname(os.path.realpath(__file__)), 'assembled_contigs.fasta')
    parser = argparse.ArgumentParser(
        description='Simple Python de Brujin graph-based assembler adapted from https://github.com/pmelsted/dbg')
    parser.add_argument("files", nargs='+', help="The seq files to generate contigs from")
    parser.add_argument("-k", "--kmer_length", type=int, help="The length of kmer to use")
    parser.add_argument("-o", "--output_path", help="Full path to which the output fasta should be written",
                        default=default_output_dir)
    parser.add_argument("-t", "--abund_threshold", type=int,
                        help="The minimum abundance at which a kmer must be found to be used in the assembly. "
                             "Given as an int. [0]", default=0)
    is_relaxed = parser.add_mutually_exclusive_group(required=False)
    is_relaxed.add_argument("--relaxed",
                        help="If relax is true then the assembler will allow contigs to continue to be build "
                             "if there are more than 1 candidate extension kmers. In such cases, it will move "
                             "the extension forwards with the most abundant kmer. [False]", action="store_true",
                        default=False)
    is_relaxed.add_argument("--strict",
                        help="If strict is true then the assembler will not allow contigs to coninue to be built"
                             "if there are more than 1 candidate extension kmers. If no flag is passed"
                             "kmers will be used in order of abundance with more abundant kmers used first"
                             " [True]", action="store_true",
                        default=True)
    parser.add_argument("--pick_by_random",
                        help="If this flag is passed, the seed kmers from which the contigs will be assembled"
                             "will be picked randomly from the count dict. [False]", action="store_true",
                        default=False)
    parser.add_argument("-v", "--verbosity", action="store_true", help="Enable a more verbose output [False]",
                        default=False)
    parser.add_argument("--remove",
                        help="When set, this flag will mean that every kmer in the dictionary may stricly only be used"
                             "once.", action="store_true", default=False)
    parser.add_argument("--assess_directionality",
                        help="If set to true, after a contig is generated from a seed kmer, rather than directly "
                             "returning that contig, the first and last kmer of this contig are used as seed kmers"
                             " to make two additional contigs, one representing extension in the forward direction "
                             "and one representing extension in the reverse direction. These two contigs are assessed "
                             "for similarity. If they are identical, one of the contigs will be returned. If there are"
                             " differences between the contigs, the contig with the highest average per kmer support"
                             " will be returned.", action="store_true", default=False)
    return parser.parse_args()


if __name__ == "__main__":
    args = process_args()
    dbga = PyDBGAssembler(
        fa_or_fq_file_path_list=args.files, kmer_len=args.kmer_length,
        output_path=args.output_path, is_relaxed=args.relaxed, is_pick_random=args.pick_by_random,
        verbose=args.verbosity, remove=args.remove, low_abund_filter_threshold=args.abund_threshold,
        assess_directionality=args.assess_directionality)
    dbga.do_assembly()
