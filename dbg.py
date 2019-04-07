#!/usr/bin/env python3
import collections
from Bio import Seq, SeqIO
import os
import argparse
import sys

class PyDBGAssembler:
    def __init__(
            self, fa_or_fq_file_path_list, kmer_len, is_relaxed, is_pick_random, verbose, remove,
            low_abund_filter_threshold, assess_directionality, no_rev_comp, no_kmer_list_output, output_dir):
        self.cwd = os.path.dirname(os.path.realpath(__file__))
        self.verbose = verbose
        self.no_rev_comp = no_rev_comp
        self.no_kmer_list_output = no_kmer_list_output
        self.is_relaxed = is_relaxed
        if kmer_len is None:
            if self.is_relaxed:
                self.kmer_len = 20
            else:
                self.kmer_len = 31
        else:
            self.kmer_len = kmer_len
        self.fa_or_fq_file_path_list = fa_or_fq_file_path_list
        if output_dir is not None:
            self.output_dir = output_dir
            self.kmer_contig_output_path = os.path.join(self.output_dir, 'assembled_contigs.fasta')
            self.kmer_list_output_path = os.path.join(self.output_dir, 'kmer_list_{}.tsv'.format(self.kmer_len))
            if not os.path.exists(self.output_dir):
                os.makedirs(self.output_dir)
        else:
            self.output_dir = os.path.abspath(os.path.dirname(fa_or_fq_file_path_list[0]))
            self.kmer_contig_output_path = os.path.join(self.output_dir, 'assembled_contigs.fasta')
            self.kmer_list_output_path = os.path.join(self.output_dir, 'kmer_list_{}.tsv'.format(self.kmer_len))
        self.completed_contigs_as_fasta = []
        # The minimum abundance of a kmer for it not to be discarded
        self.abund_filter_cutoff = int(low_abund_filter_threshold)
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
        self.current_seed_kmer = None
        self.remove = remove
        self.assess_directionality = assess_directionality
        # Attributes used in directionality assessment
        self.lhs_seed_kmer = None
        self.rhs_seed_kmer = None
        self.forward_contig_list = None
        self.reverse_contig_list = None
        # The set of kmers that have been used in building the current contig
        self.temp_used_build_kmers_set = set()
        self.global_contig_count = 0

    def do_assembly(self):
        print("Starting assembly...\n\n\n")
        if not self.verbose:
            sys.stdout.write('.')
            count = 0
        one_percent = int(self.num_kmers / 100)
        for kmer_key in self.kmer_list_for_contigs:
            if not self.verbose:
                if count%one_percent == 0:
                    sys.stdout.write('.')
                    sys.stdout.flush()
                count += 1
            if kmer_key not in self.used_kmers_set:
                self.current_seed_kmer = kmer_key
                contig_as_string, contig_kmer_list = self._get_contig()
                for used_kmer in contig_kmer_list:
                    self.used_kmers_set.add(used_kmer)
                    if not self.no_rev_comp:
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
        print("Writing out contig fasta to {}".format(self.kmer_contig_output_path))
        with open(self.kmer_contig_output_path, 'w') as f:
            for line in self.completed_contigs_as_fasta:
                f.write('{}\n'.format(line))
        print("Done.")
        if self.verbose:
            print('Kmer lists was written to {}'.format(self.kmer_list_output_path))

    def _contig_list_to_fasta(self):
        for i, contig_seq in enumerate(self.list_of_completed_contigs):
            self.completed_contigs_as_fasta.extend(['>{}'.format(i), '{}'.format(contig_seq)])

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
        self.global_contig_count += 1
        if self.verbose:
            print('>seed_kmer_{}'.format(self.global_contig_count))
            print(self.current_seed_kmer)
        self.kmer_list_of_contig_fwd = self._extend_contig(self.current_seed_kmer, dir='fwd')
        if self.verbose:
            print('>forward_extension_{}'.format(self.global_contig_count))
            print(self._contig_to_string(self.kmer_list_of_contig_fwd))
        self.kmer_list_of_contig_rev = self._extend_contig(self.current_seed_kmer, dir='rev')
        if self.verbose:
            print('>reverse_extension_{}'.format(self.global_contig_count))
            print(self._contig_to_string(self.kmer_list_of_contig_rev))

        # I think perhaps it is looking for cases where we have created a circular piece of DNA.
        # In this case there is no need to search for the reverse build as well
        if self._contig_is_circular(self.kmer_list_of_contig_fwd[-1]):
            self.kmer_list_of_contig = self.kmer_list_of_contig_fwd
            if self.verbose:
                print('>assembled_contig_{}'.format(self.global_contig_count))
                print(self._contig_to_string(self.kmer_list_of_contig))
                print('resultant kmer is circular. no directionality assessment\n\n\n')
        else:
            self.kmer_list_of_contig = self.kmer_list_of_contig_rev[:-1] + self.kmer_list_of_contig_fwd
            if self.verbose:
                print('>assembled_contig_{}'.format(self.global_contig_count))
                print(self._contig_to_string(self.kmer_list_of_contig))
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
                if self.verbose:
                    print('{}: no fwd directionality extension. '
                          'Returning rev directionality extension.\n\n\n'.format(self.global_contig_count))
                return self._contig_to_string(self.reverse_contig_list), self.reverse_contig_list
            elif self.forward_contig_list:
                self._del_temp_kmers_from_count_dict(self.forward_contig_list)
                if self.verbose:
                    print('{}: no rev directionality extension. '
                          'Returning fwd directionality extension.\n\n\n'.format(self.global_contig_count))
                return self._contig_to_string(self.forward_contig_list), self.forward_contig_list
            else:  # neither list is good
                raise RuntimeError('Ooops Neither directionality extension exists')
        return self._if_agree_return_fwd_contig_else_return_contig_with_most_support()

    def _if_agree_return_fwd_contig_else_return_contig_with_most_support(self):
        if self.forward_contig_list == self.reverse_contig_list:  # If no disagreement then contig is good
            if self.verbose:
                print('{}: fwd and rev directionality extensions are identical. '
                      'Returning fwd extension.\n\n\n'.format(self.global_contig_count))
            if self.remove:
                self._del_temp_kmers_from_count_dict()
            return self._contig_to_string(self.forward_contig_list), self.forward_contig_list
        else:  # contigs disagree then take the contig that has the highest average support for its kmers
            if self.verbose:
                print('{}: fwd and rev directionality extensions are NOT identical. '
                      'Comparing extensions.'.format(self.global_contig_count))
            return self._return_contig_with_highest_average_kmer_support()

    def _make_reverse_extension_contig(self):
        self.rhs_seed_kmer = self.kmer_list_of_contig[-1]
        if self.verbose:
            print('>rev_directionality_extension_seed_{}'.format(self.global_contig_count))
            print(self.rhs_seed_kmer)
        # reset the used kmer set
        self._reset_temp_used_kmer_set()
        self.reverse_contig_list = self._extend_contig(self.rhs_seed_kmer, dir='rev')
        if self.verbose:
            print('>rev_directionality_extension_{}'.format(self.global_contig_count))
            print(self._contig_to_string(self.forward_contig_list))

    def _make_forward_extension_contig(self):
        self.lhs_seed_kmer = self.kmer_list_of_contig[0]
        if self.verbose:
            print('>fwd_directionality_extension_seed_{}'.format(self.global_contig_count))
            print(self.lhs_seed_kmer)
        # reset the used kmer set
        self._reset_temp_used_kmer_set()
        self.forward_contig_list = self._extend_contig(self.lhs_seed_kmer, dir='fwd')
        if self.verbose:
            print('>fwd_directionality_extension_{}'.format(self.global_contig_count))
            print(self._contig_to_string(self.forward_contig_list))

    def _return_contig_with_highest_average_kmer_support(self):
        forward_mean_kmer_support = self._calc_mean_kmer_support(self.forward_contig_list)
        reverse_mean_kmer_support = self._calc_mean_kmer_support(self.reverse_contig_list)
        if self.verbose:
            print('{}: fwd directionality extension average kmer support: {}'.format(self.global_contig_count, forward_mean_kmer_support))
            print('{}: rev directionality extension average kmer support: {}'.format(self.global_contig_count, reverse_mean_kmer_support))
        if forward_mean_kmer_support > reverse_mean_kmer_support:
            if self.remove:
                self._del_temp_kmers_from_count_dict(self.forward_contig_list)
            if self.verbose:
                print('{}: returning fwd directionality extension as contig\n\n\n'.format(self.global_contig_count))
            return self._contig_to_string(self.forward_contig_list), self.forward_contig_list
        else:
            if self.remove:
                self._del_temp_kmers_from_count_dict(self.reverse_contig_list)
            if self.verbose:
                print('{}: returning rev directionality extension as contig\n\n\n'.format(self.global_contig_count))
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
                # no need to delete the rev comp as this has already been added to the temp_used...
        else:
            for kmer in kmer_set_of_list:
                del self.kmer_count_dict[kmer]
                if not self.no_rev_comp:
                    if self._rev_comp(kmer) != kmer: # cannot delete a seq if rev comp is same as fwd comp.
                        del self.kmer_count_dict[self._rev_comp(kmer)]
        self.temp_used_build_kmers_set = set()

    def _contig_is_circular(self, kmer):
        return self.current_seed_kmer in self._fwd_seq_generator(kmer)

    def _extend_contig(self, kmer, dir):
        """Dir is the direction that we are extending. Previously we would work irrespective of direction by simply
        rev complementing the kmer from the forward extension. But we now want to be able to work in a system that
        doesn't automatically use the rev comp reads in the count dict. The previous approch of extending the rev
        comp of the seed kmer relied on all kmers being in the count dict as both their fwd and rev complement. By
        specificng the direction we are working in we can make ourselves independent of needing to have the rev comp
        of kmers in the count dict."""
        kmers_of_contig = [kmer]

        while True:
            # cand is abbrev. for candidate
            if dir == 'fwd':
                cand_kmers = [(cand, self.kmer_count_dict[cand]) for cand in
                              self._fwd_seq_generator(kmers_of_contig[-1]) if
                              cand in self.kmer_count_dict if
                              cand not in self.temp_used_build_kmers_set]
            else:  # dir == 'rev'
                cand_kmers = [(cand, self.kmer_count_dict[cand]) for cand in
                              self._rev_seq_generator(kmers_of_contig[0]) if
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

            if candidate == self._rev_comp(kmers_of_contig[-1]):
                break  # break out of hairpins

            if not self.is_relaxed:
                # This check will only be made if we are operating in --strict mode
                if sum(x in self.kmer_count_dict for x in self._rev_seq_generator(candidate)) != 1:
                    break
            else:
                if candidate == kmers_of_contig[-1]:
                    if self.verbose:
                        print('Breaking contig assembly due to homopolymer {}'.format(candidate))
                    break
            self._put_candidate_in_kmer_list_for_contig(candidate, dir, kmers_of_contig)

            if self.remove:
                self.temp_used_build_kmers_set.add(candidate)
                if not self.no_rev_comp:
                    if self._rev_comp(candidate) != candidate:
                        self.temp_used_build_kmers_set.add(self._rev_comp(candidate))

        return kmers_of_contig

    def _put_candidate_in_kmer_list_for_contig(self, candidate, dir, kmers_of_contig):
        """If we are building in the fwd direction we can simply append. If we are working in the rev direction
        then we should put the contig at the begining of the list. This way there is no reordering of the list
        required before appending the two lists"""
        if dir == 'fwd':
            kmers_of_contig.append(candidate)
        else:  # dir == 'rev':
            kmers_of_contig.insert(0, candidate)

    def _init_kmer_count_dict(self):
        """Create the dictionary that will hold the kmer counts"""
        print('Creating kmer count dict')
        for file_path in self.fa_or_fq_file_path_list:
            self._get_list_of_reads_from_fq_or_fa_file(file_path)
            self._break_reads_into_kmers_and_pop_dict()
        self._remove_low_abund_kmers()
        print('Done.')
        self._write_out_kmer_lists()

    def _write_out_kmer_lists(self):
        """Write out the kmers and their abundance as kmer sequence \t kmer abundance
        with one kmer abundance pair per line. If we have verbose set then we should print out the writing feedback
        at the end of the assembly rather than here else this will be obscured by all the output to stdout.
        """
        if not self.no_kmer_list_output:
            if not self.verbose:
                print('Writing kmer lists to {}'.format(self.kmer_list_output_path))
            with open(self.kmer_list_output_path, 'w') as f:
                for kmer, abund in sorted(self.kmer_count_dict.items(), key=lambda x:x[1], reverse=True):
                    f.write('{}\t{}\n'.format(kmer, abund))
            if not self.verbose:
                print('Done')

    def _break_reads_into_kmers_and_pop_dict(self):
        for read in self.reads_of_file:
            seq_str = self._get_read_sequence_as_str(read)
            list_of_sub_seqs = self._split_str_into_multiple_reads_by_n(seq_str)
            for seq in list_of_sub_seqs:
                for km in self._create_kmers_in_sequence_generator(seq):
                    self.kmer_count_dict[km] += 1
                    self.total_number_of_kmer_seqs += 1
                if not self.no_rev_comp:
                    rev_seq = self._rev_comp(seq)
                    for km in self._create_kmers_in_sequence_generator(rev_seq):
                        self.kmer_count_dict[km] += 1
                        self.total_number_of_kmer_seqs += 1
        print('{} total kmers found; '
              'represented as {} unique kmers'.format(self.total_number_of_kmer_seqs, len(self.kmer_count_dict)))

    @staticmethod
    def _split_str_into_multiple_reads_by_n(seq_s):
        return seq_s.split('N')

    @staticmethod
    def _get_read_sequence_as_str(read):
        return str(read.seq)

    def _remove_low_abund_kmers(self):
        low_abund_kmers = [x for x in self.kmer_count_dict if
                           self.kmer_count_dict[x] < self.abund_filter_cutoff]
        print('{} unique kmers to be removed at cutoff threshold of {}'.format(len(low_abund_kmers), self.abund_filter_cutoff))
        print('removing low abund kmers from count dict...')
        for kmer in low_abund_kmers:
            del self.kmer_count_dict[kmer]

        print('{} total kmers and {} unique kmers now remaining '
              'in count dict'.format(sum(self.kmer_count_dict.values()), len(self.kmer_count_dict)))
        print('Done.')
        print('Count dictionary population complete.')

    def _get_list_of_reads_from_fq_or_fa_file(self, file_path):
        file_ext = os.path.splitext(file_path)[1]
        if file_ext == '.fq' or file_ext == '.fastq':
            self.reads_of_file = SeqIO.parse(file_path, 'fastq')
        elif file_ext == '.fasta' or file_ext == '.fa':
            self.reads_of_file = SeqIO.parse(file_path, 'fasta')
        else:
            raise RuntimeError('Unknown extension of file {}'.format(file_path))

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
    parser = argparse.ArgumentParser(
        description='Simple Python de Brujin graph-based assembler adapted from https://github.com/pmelsted/dbg')
    parser.add_argument("files", nargs='+', help="The seq files to generate contigs from")
    parser.add_argument("-k", "--kmer_length", type=int, help="The length of kmer to use")
    parser.add_argument("-o", "--output_dir",
                        help="Directory to which the output fasta and kmer lists should be written")
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
    parser.add_argument("-v", "--verbose", action="store_true", help="Enable a more verbose output [False]",
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
    parser.add_argument("--no_rev_comp",
                        help="When this flag is passed, the kmer dict will only be made from the sequences that are "
                             "contined in the fasta and fastq files passed to the assembler. i.e. no reverse"
                             " complements of the kmers will be added to the kmer dictionary. Reverse complements"
                             " of kmers will therefore not be available for making contigs. [False]", default=False,
                        action="store_true")
    parser.add_argument("--no_kmer_list_output", help="When this flag is passed. No kmer lists will be printed out. "
                                                      "[False]", action="store_true", default=False)
    # return parser.parse_args(["--no_rev_comp", "--remove", "--relaxed", "--verbose",
    # "--assess_directionality", "../data/test_fastq.fastq"])
    return parser.parse_args()


if __name__ == "__main__":
    args = process_args()
    dbga = PyDBGAssembler(
        fa_or_fq_file_path_list=args.files, kmer_len=args.kmer_length,
        output_dir=args.output_dir, is_relaxed=args.relaxed, is_pick_random=args.pick_by_random,
        verbose=args.verbose, remove=args.remove, low_abund_filter_threshold=args.abund_threshold,
        assess_directionality=args.assess_directionality, no_rev_comp=args.no_rev_comp,
        no_kmer_list_output=args.no_kmer_list_output)
    dbga.do_assembly()
