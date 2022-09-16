import os
import re
import sys
import itertools
from Bio.Seq import Seq

from minorg.functions import splitlines
from minorg.fasta import fasta_to_dict, dict_to_fasta
from minorg.minweight_sc import (
    enum_approx_order_SC,
    approx_min_SC,
    Set,
    SetOfSets
)

# # ## test
# import minorg
# from minorg import grna, minimum_set_v2, minweight_sc
# x = grna.gRNAHits()
# x.parse_from_mapping("/mnt/chaelab/rachelle/zzOtherzz/Jinge/nrg_adr/nrg1ab_vdw_minlen1500/nrg1ab_vdw/nrg1ab_vdw/nrg1ab_vdw_gRNA_all.map")
# x2 = x.filter_hits(all_checks = True)
# x2 = x2.filter_seqs(all_checks = True)
# y2 = minorg.minimum_set_v2.collapse_gRNAHits(x2)
# z2 = minorg.minimum_set_v2.limited_minweight_SC(y2, 10)
# z2_2 = [x for x in z2 if len(x) == 2] ## get set cover solutions w/ 2 sets
# [[a for a in sos] for sos in z2_2] ## show CollapsedgRNA
# z2_2[1].get('015')[0].most_5prime().relative_distance_to_5prime(z2_2[0].get('008')[0].most_5prime()) ## compare

impossible_set_message_default = "If you used the full programme, consider adjusting --minlen, --minid, and --merge-within, and/or raising --check-reciprocal to restrict candidate target sequences to true biological replicates. Also consider using --domain to restrict the search."

def manual_check_prompt(grnas, set_num = None):
    """
    Prints prompt for manual check.
    
    Arguments:
        grnas (list): list of :class:`~minorg.grna.gRNASeq` or :class:`~minorg.minimum_set.gRNA` objects.
            Should already be sorted in printing order.
    
    Returns
    -------
    str
        Of user input
    """
    if set_num is not None:
        print(f"\n\tID\tsequence (Set {set_num})")
    else:
        print("\n\tID\tsequence")
    for grna in grnas:
        print(f"\t{grna.id}\t{grna.seq}")
    ## obtain user input
    usr_input = input(("Hit 'x' to continue if you are satisfied with these sequences."
                       " Otherwise, enter the sequence ID (case-sensitive) or sequence of"
                       " an undesirable gRNA and hit the return key to update this list: "))
    return usr_input

class gRNA(Set):
    def __init__(self, name, gRNASeq_obj):
        self.gRNASeq = gRNASeq_obj
        target_ids = set(grnahit.target_id for grnahit in gRNASeq_obj.hits)
        ## weight is equivalent to number of targets covered
        super().__init__(name = name, weight = len(target_ids), elements = target_ids)
    def __repr__(self):
        return (f"gRNA(name='{self.name}', seq='{self.seq}',"
                f" num_hits={len(self.hits)}, coverage={self.weight})")
    @property
    def id(self) -> str: return self.gRNASeq.id
    @property
    def seq(self) -> str: return self.gRNASeq.seq
    @property
    def targets(self) -> list: return self.gRNASeq.targets
    @property
    def coverage(self) -> int: return self.weight
    @property
    def hits(self) -> list: return list(self.gRNASeq.hits)
    @property
    def relative_5prime_pos(self) -> float:
        """
        Returns a value calculated from all hits that is
        ONLY VALID FOR COMPARISON BETWEEN gRNA with the SAME COVERAGE.
        The smaller the value, the closer to the 5' end.
        """
        return self._relative_5prime_pos(self.hits)
    def _relative_5prime_pos(self, gRNAHit_objs) -> float:
        """
        Returns a value calculated from hits in 'gRNAHit_objs' that is
        ONLY VALID FOR COMPARISON BETWEEN gRNA with the SAME COVERAGE.
        The smaller the value, the closer to the 5' end.
        """
        return sum((hit.range[0] if hit.target.sense != '-'
                    else (hit.target_len - hit.range[1]))
                   for hit in gRNAHit_objs)/len(gRNAHit_objs)
    def subset_hits_by_target(self, target_ids) -> list:
        """
        Returns list of :class:`~minorg.grna.gRNAHit` objects where the hit is
        to a target with the same name as at lest one entry in 'target_ids'.
        """
        target_ids = set(target_ids)
        return [hit for hit in self.hits if hit.target_id in target_ids]
    def relative_distance_to_5prime(self, other) -> float:
        """
        Takes another gRNA object and returns self._relative_5prime_pos(<hits in shared targets>) - other._relative_5prime_pos(<hits in shared targets>).
        If self, is closer to the 5' end, the value returned will be negative.
        
        Returns
        -------
        float
        
        Raises
        ------
        Exception
            If no common targets between the two gRNA objects (cannot compare across different targets)
        """
        common_targets = self.common_targets(other)
        if not common_targets:
            raise Exception("Error: No common targets")
        self_pos = self._relative_5prime_pos(self.subset_hits_by_target(common_targets))
        other_pos = other._relative_5prime_pos(other.subset_hits_by_target(common_targets))
        return self_pos - other_pos
    def closer_to_5prime(self, other) -> bool:
        """
        Takes another gRNA object and returns True if self is closer to 5' end
        (of sense strand, if sense information exists) by
        comparing relative positions between shared targets.
        If tie, return True.
        
        Returns
        -------
        bool
        
        Raises
        ------
        Exception
            If no common targets between the two gRNA objects (cannot compare across different targets)
        """
        return self.relative_distance_to_5prime(other) <= 0
    def common_targets(self, other):
        return self.intersection(other)

class CollapsedgRNA(Set):
    def __init__(self, name, gRNA_objs):
        self.gRNAs = set(gRNA_objs)
        ## weight is equivalent to number of targets covered
        super().__init__(name = name, weight = len(gRNA_objs[0]), elements = gRNA_objs[0])
    def __repr__(self):
        return (f"CollapsedgRNA(name='{self.name}', num_gRNA={len(self.gRNAs)}, coverage={self.weight})")
    def most_5prime(self):
        """
        Returns :class:`~minorg.minimum_set.gRNA` object closest to 5'
        
        Returns
        -------
        :class:`~minorg.minimum_set.gRNA`
        """
        return min(self.gRNAs, key = lambda grna:grna.relative_5prime_pos)
    def least_5prime(self):
        """
        Returns :class:`~minorg.minimum_set.gRNA` object furthest from 5'
        (Usually but not necessarily closest to 3')

        Returns
        -------
        :class:`~minorg.minimum_set.gRNA`
        """
        return max(self.gRNAs, key = lambda grna:grna.relative_5prime_pos)
    def remove(self, gRNA_obj) -> None:
        """
        Remove a single :class:`~minorg.minimum_set.gRNA` object from self.
        
        Arguments:
            gRNA_obj (:class:`~minorg.minimum_set.gRNA`): gRNA object
        """
        if gRNA_obj in self.gRNAs:
            self.gRNAs.remove(gRNA_obj)
    def add(self, gRNA_obj) -> None:
        """
        """
        if set(gRNA_obj) != set(self):
            raise Exception("Error: Cannot add gRNA to CollapsedgRNA. Coverage not identical.")
        self.gRNAs.add(gRNA_obj)
    def copy(self):
        return self.__class__(self.name, list(self.gRNAs))
    def num_grna(self):
        return len(self.gRNAs)
    def is_empty(self):
        return self.num_grna() == 0

class SetOfCollapsedgRNA(SetOfSets):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
    def __repr__(self):
        return (f"SetOfCollapsedgRNA(num_CollapsedgRNA={len(self)}, weight={self.w})")
    @property
    def max_coverage(self): return len(max(self, key = lambda cg:len(cg)))
    @property
    def min_coverage(self): return len(min(self, key = lambda cg:len(cg)))
    def copy(self):
        return self.__class__(*self.sets)
    def write(self, fout):
        with open(fout, 'w+') as f:
            f.write('\t'.join(["coverage group", "coverage", "gRNA id",
                               "gRNA sequence", "relative pos"]) + '\n')
            for collapsed_grna in sorted(self, key = lambda cg:-len(cg)):
                for grna in sorted(collapsed_grna.gRNAs, key = lambda g:g.relative_5prime_pos):
                    f.write('\t'.join(map(str,
                        [collapsed_grna.name, len(collapsed_grna), grna.id,
                         grna.seq, round(grna.relative_5prime_pos, 2)])) + '\n')
        return
    def remove_empty(self):
        self.remove(*self.empty_collapsed_grna())
    def all_not_empty(self):
        return all((not x.is_empty()) for x in self)
    def empty_collapsed_grna(self):
        return [x for x in self if x.is_empty()]
    def remove_grna(self, *gRNAs):
        for cg in self:
            for grna in gRNAs:
                cg.remove(grna)
        return
    def add_grna(self, *gRNAs):
        added = set()
        for grna in gRNAs:
            for cg in self:
                if cg == grna:
                    cg.add(grna)
                    added.update({grna})
        if set(gRNAs) - added != set():
            print(("The following gRNA object(s) was/were not added due to incompatible coverage:"
                   f" {','.join(grna.name for grna in (set(gRNAs) - added))}"))
        return
    def _generate_grna_set(self, collapsed_grna, prioritise_3prime = False, resort = True):
        """
        Generates a single gRNA set from list of CollapsedgRNA objects in 'collapsed_grna'.
        
        Arguments:
            collapsed_grna (iter): of :class:`~minorg.minimum_set.CollapsedgRNA` objects
            prioritise_3prime (bool): tie-break with proximity to 3' end instead of 5' (default=False)
            resort (bool): sort 'collapsed_grna' by descending order of coverage size,
                tie-breaks by sorted target names
        
        Returns:
            list: Of :class:`~minorg.minimum_set.gRNA`, one from each :class:`~minorg.minimum_set.CollapsedgRNA`
        """
        if resort:
            sorted_collapsed_grna = sorted(collapsed_grna, key = lambda cg:(-len(cg), sorted(cg)))
        else:
            sorted_collapsed_grna = list(collapsed_grna)
        ## copy the CollapsedgRNA objects so we can modify them without changing the original
        if not all(len(cg.gRNAs) > 0 for cg in sorted_collapsed_grna):
            raise Exception("Error: Unable to make gRNA set. Some CollapsedgRNA do not have associated gRNA.")
        return [cg.most_5prime() if not prioritise_3prime else cg.least_5prime()
                for cg in sorted_collapsed_grna]
    def generate_grna_set(self, prioritise_3prime = False, consume = False):
        """
        Outputs list of 1 gRNA from each CollapsedgRNA in self. 
        By default, the gRNA closest to the 5' (sense, if information exists) strand are selected from each
        CollapsedgRNA object.
        
        Arguments:
            prioritise_3prime (bool): select gRNA from CollapsedgRNA using proximity to 3' end
                instead of 5' end (default=False)
            consume (bool): remove output gRNA from self's CollapsedgRNA objects permanently
        
        Returns:
            list: Of :class:`~minorg.minimum_set.gRNA` objects
        """
        grnas = self._generate_grna_set(self, resort = True, prioritise_3prime = prioritise_3prime)
        if consume:
            self.remove_grna(*grnas)
            # print(self.all_not_empty())
        return grnas
    def generate_grna_sets(self, prioritise_3prime = False, set_num = 1, max_set_num = 1, manual = True):
        """
        Yields lists of 1 gRNA from each CollapsedgRNA in self.
        
        Arguments:
            prioritise_3prime (bool): tie-break with proximity to 3' end instead of 5' (default=False)
            set_num (int): number to print for manual check. Also used with 'max_set_num' to set limit
                on number of sets to generate. (default=1)
            max_set_num (int): maximum set number. Maximum of sets generated will be
                '<max_set_num> - <set_num> + 1' (default=1)
            manual (bool): enable manual screening of each gRNA set for approval at interactive terminal
                (default=False)
        
        Returns
        -------
        generator
            Of [<1 :class:`~minorg.minimum_set.gRNA` object from each CollapsedgRNA in self>]
        """
        collapsed_grna = sorted(self, key = lambda cg:(-len(cg), sorted(cg)))
        ## copy the CollapsedgRNA objects so we can modify them without changing the original
        collapsed_grna = [cg.copy() for cg in collapsed_grna]
        output = []
        while all(len(cg.gRNAs) > 0 for cg in collapsed_grna) and set_num <= max_set_num:
            selected_grna = self._generate_grna_set(collapsed_grna, resort = False,
                                                    prioritise_3prime = prioritise_3prime)
            if manual:
                sorted_grna = sorted(selected_grna, key = lambda grna:grna.id)
                usr_input = manual_check_prompt(sorted_grna, set_num = set_num)
                if usr_input.upper() != 'X':
                    id_list = tuple(grna.id for grna in selected_grna)
                    seq_list = tuple(str(grna.seq).upper() for grna in selected_grna)
                    if usr_input in id_list:
                        index = id_list.index(usr_input)
                        collapsed_grna[index].remove(selected_grna[index])
                    elif usr_input.upper() in seq_list:
                        index = seq_list.index(usr_input)
                        collapsed_grna[index].remove(selected_grna[index])
                    else:
                        print("Invalid input.")
                    continue ## skip removal of selected gRNA and yielding of selected gRNA set
            ## remove selected gRNA from candidate pool
            for i, grna in enumerate(selected_grna):
                collapsed_grna[i].remove(grna)
            set_num += 1
            yield selected_grna
        # return output
        return

def collapse_gRNAHits(gRNAHits):
    """
    Create :class:`SetOfCollapsedgRNA` from :class:`~minorg.grna.gRNAHits` object.
    """
    grna_seqs = gRNAHits.gRNAseqs.values()
    equivalents = {}
    for grna_seq in grna_seqs:
        gRNA_obj = gRNA(grna_seq.id, grna_seq)
        cov = tuple(sorted(gRNA_obj))
        equivalents[cov] = equivalents.get(cov, []) + [gRNA_obj]
    ## make CollapsedgRNA objects
    sorted_equivalents = sorted(equivalents.items(), key = lambda x:(-len(x[0]), x[0]))
    ## name collapsed gRNA according to sorted order
    max_digits = len(str(len(sorted_equivalents)))
    return SetOfCollapsedgRNA(*[CollapsedgRNA("cg_"+str(i+1).zfill(max_digits), dat[1])
                                for i, dat in enumerate(sorted_equivalents)])

def limited_minweight_SC(collapsed_grnas, num_sets, targets = None,
                         num_lengths_to_track = None,
                         low_coverage_penalty = 0):
    """
    Executes :func:`~minorg.minweight_sc.enum_approx_order_SC` for a capped number of iterations
    (max(20, 2*num_sets)) while seeding each run with a different CollapsedgRNA, starting with the
    CollapsedgRNA with the highest coverage to the lowest. CollapsedgRNA will be removed from candidate
    list once all CollapsedgRNA with equivalent coverage have been used as seed. Stops when coverage of
    the next CollapsedgRNA to be seeded has a coverage of less than <total targets>/<num CollapsedgRNA
    in <num_lengths_to_track>th smallest set cover solution>.
    
    Arguments:
        collapsed_grnas (CollapsedgRNA): :class:`CollapsedgRNA`CollapsedgRNA object
        num_sets (int): desired number of sets. Used to inform maximum number of iterations.
        targets (list or set or tuple): targets IDs (str) of targets to be coverd
        num_lengths_to_track (int): length of <num_lengths_to_track>th smallest set cover solution
            will be used to determine whether to terminate search (see function description)
        low_coverage_penalty (float): multiplier for value calculated by :func:`~minorg.minweight_sc`,
            which will then be multiplied by <number of remaining targets that are not covered by gRNA> 
            and then added to the output of that function. Effectively penalises large sets of many
            small coverage gRNA. This might make the set less redundant, but will likely reduce set size.
    
    Returns
    -------
    list
        Of set cover solutions (:class:`CollapsedgRNA`)
    """
    num_iter = max(20, 2*num_sets)
    if targets is None:
        targets = collapsed_grnas.elements
    else:
        targets = set(targets)
    num_targets = len(targets)
    if num_lengths_to_track is None:
        num_lengths_to_track = max(10, 2*num_sets)
    sc_solutions = []
    min_lens = [num_targets]
    ## group gRNA with identical coverage
    cov_grna = {}
    for collapsed_grna in collapsed_grnas:
        cov = len(collapsed_grna)
        cov_grna[cov] = cov_grna.get(cov, []) + [collapsed_grna]
    ## execute enum_approx_order_SC by seeding with gRNA obj from highest coverage to lowest
    ## once all gRNA with given coverage size have been seeded, remove all from collapsed_grnas
    ## if cov < num_targets/len(<num_sets>th smallest solution), stop
    min_len = lambda *args: min_lens[-1]
    terminate_loop = lambda cov: cov < (1 if len(min_lens) < num_lengths_to_track
                                        else num_targets/min_len())
    retain_min_len = lambda Cs: [C for C in Cs if len(C) <= min_len()]
    def update_min_lens(Cs):
        ml = min_len()
        candidates = [len(C) for C in Cs if len(C) <= ml]
        output = sorted(min_lens + candidates)[:num_lengths_to_track]
        sample = set() if not Cs else Cs[0]
        # print(ml, type(sample), len(sample), [C for C in Cs if len(C) == 1][:2],
        #       len(candidates), len(min_lens), output)
        return output
    ## make a copy of collapsed_grnahits that we can modify
    collapsed_grnas = collapsed_grnas.copy()
    ## remove collapsed_grna if it covers ALL targets
    if num_targets in cov_grna:
        ## there should only be one entry in cov_grna[num_targets] since gRNA w/ same coverage are collapsed
        ## and there's only possible coverage that well covers all targets but for the sake of consistency
        ## we'll just go ahead and assume unknown length
        collapsed_grnas.remove(*cov_grna[num_targets])
        new_solutions = [SetOfCollapsedgRNA(collapsed_grna) for collapsed_grna in cov_grna[num_targets]]
        sc_solutions.extend(new_solutions)
        min_lens = update_min_lens(new_solutions)
        del cov_grna[num_targets]
    ## start set cover :)
    for cov, gRNAs in sorted(cov_grna.items(), key = lambda x:-x[0]):
        # print("coverage:", cov)
        if terminate_loop(cov): break
        curr_solutions = []
        # print(type(gRNAs[0]))
        for grna in gRNAs:
            solutions = enum_approx_order_SC(targets, collapsed_grnas,
                                             seed = grna, num_iter = num_iter,
                                             low_coverage_penalty = low_coverage_penalty)
            curr_solutions.extend(solutions)
        min_lens = update_min_lens(curr_solutions)
        sc_solutions = list(set(retain_min_len(sc_solutions) + \
                                retain_min_len([C for C in curr_solutions if C not in sc_solutions])))
        collapsed_grnas.remove(*gRNAs)
    return sc_solutions

## gets a single minimum set
def make_get_minimum_set(gRNA_hits, manual_check = True, exclude_seqs = set(), targets = None,
                         prioritise_nr = False, sc_algorithm = "LAR", num_sets = 1, tie_breaker = None,
                         low_coverage_penalty = 0.5, suppress_warning = False,
                         impossible_set_message = impossible_set_message_default):
    """
    Make function to generate minimum set of gRNA.
    
    Arguments:
        gRNA_hits (list): gRNAHit objects
        manual_check (bool): manually approve each gRNA set
        exclude_seqs (set/list): optional, gRNA sequences (str) to exclude
        targets (list): optional, target IDs (str)
        prioritise_nr (bool): prioritise non-redundancy.
            If used, 'sc_algorithm' will be ignored. (default=False)
        sc_algorithm (str): set cover algorithm when not prioritising non-redundancy.
            Only used if prioritise_nr=False. (default='LAR')
        num_sets (int): number of sets
        tie_breaker (func): tie-breaker function.
            Takes (1) 'gRNA_coverage' filtered for unselected gRNA seq,
            (2) unmodified 'gRNA_coverage',
            (3) list of IDs of targets covered by already selected gRNA
        impossible_set_message (str): message to print when gRNA cannot cover all targets
        suppress_warning (bool): suppress printing of warning when gRNA cannot cover all targets
    
    Returns
    -------
    func
        That returns list of minimum set of gRNA sequences (str)
    """
    # ## filter by excluded sequences
    if exclude_seqs:
        gRNA_hits.set_seqs_check("exclude", False, [str(s) for s in exclude_seqs])
        gRNA_hits = gRNA_hits.filter_seqs_all_checks_passed(quiet = suppress_warning)
        gRNA_hits = gRNA_hits.filter_hits_all_checks_passed(quiet = suppress_warning)
    if prioritise_nr:
        set_cover = make_set_cover_nr(gRNA_hits, num_sets = num_sets, target_ids = targets,
                                      num_lengths_to_track = None,
                                      low_coverage_penalty = low_coverage_penalty,
                                      suppress_warning = suppress_warning)
    else:
        ## tie breakers should return 2 values: <gRNASeq>, [<gRNAHits>]
        ## note: If antisense, tie break by minimum -end. Else, tie break by minimum start.
        ## note: tie-breaker uses AVERAGE distance of hits (to inferred N-terminus)
        if tie_breaker is None:
            tie_breaker = lambda *args: tuple(all_best_pos(*args).items())[0]
        set_cover = make_set_cover_pos(gRNA_hits, num_sets = num_sets, target_ids = targets,
                                       algorithm = sc_algorithm, id_key = lambda x:x.target_id,
                                       tie_breaker = tie_breaker,
                                       suppress_warning = suppress_warning)
    set_num = [0]
    def get_minimum_set():
        restore = []
        set_num[0] += 1
        while True:
            ## solve set_cover
            selected_grna = set_cover(restore = restore)
            restore = selected_grna
            ## if empty set, print message and break out of loop to exit and return the empty set
            if len(selected_grna) == 0:
                print(impossible_set_message)
                break
            ## if valid set AND manual check NOT requested, break out of loop to exit and return the valid set
            elif not manual_check: break
            ## if valid set AND manual check requested
            else:
                ## print gRNA sequences in selected_grna to screen for user to evaluate
                sorted_grna = sorted(selected_grna, key = lambda grna:grna.id)
                usr_input = manual_check_prompt(sorted_grna, set_num[0])
                if usr_input.upper() == 'X':
                    break
                else:
                    ## id_list and seq_list have same order as selected_grna
                    id_list = tuple(grna.id for grna in selected_grna)
                    seq_list = tuple(str(grna.seq).upper() for grna in selected_grna)
                    ## remove gRNA from list of gRNAs to restore upon next set cover generation
                    if usr_input in id_list:
                        index = id_list.index(usr_input)
                        restore.remove(selected_grna[index])
                    elif usr_input.upper() in seq_list:
                        index = seq_list.index(usr_input)
                        restore.remove(selected_grna[index])
                    else:
                        print("Invalid input.")
        return [str(grna.seq) for grna in selected_grna]
    return get_minimum_set

def make_set_cover_nr(gRNA_hits, num_sets = 1, target_ids = [], low_coverage_penalty = 0,
                      num_lengths_to_track = None, prioritise_3prime = False,
                      suppress_warning = False):
    """
    Create function that generates mutually exclusive gRNA sets with non-redundancy as a priority.
    
    Arguments:
        gRNA_hits (:class:`~minorg.grna.gRNAHits`): gRNAHits object
        num_sets (int): number of mutually exclusive gRNA sets to return (default=1)
        target_ids (list): list of target names/IDs (str) to be covered.
            If not provided, will be inferred from the set of target IDs covered by gRNA hits in 'gRNA_hits'
        manual (bool): manually approve each gRNA set for inclusion through interactive terminal (default=False)
        low_coverage_penalty (float): multiplier for value calculated by :func:`~minorg.minweight_sc`,
            which will then be multiplied by <number of remaining targets that are not covered by gRNA> 
            and then added to the output of that function. Effectively penalises large sets of many
            small coverage gRNA. This might make the set less redundant, but will likely reduce set size.
        prioritise_3prime (bool): tie-break with proximity to 3' end instead of 5' (default=False)
    
    Returns
    -------
    func
        That returns list (gRNA panel) of str (gRNA names)
    """
    collapsed_grnas = collapse_gRNAHits(gRNA_hits)
    if not target_ids:
        target_ids = set().union(*[set(cg) for cg in collapsed_grnas])
    else:
        target_ids = set(target_ids)
    ## function to regenerate set cover solutions from collapsed_grna object
    def generate_sc_solutions():
        ## sort in order of smallest set cover size, smallest redundancy, and size of largest set in set cover
        return sorted(limited_minweight_SC(collapsed_grnas, num_sets, targets = target_ids,
                                           low_coverage_penalty = low_coverage_penalty,
                                           num_lengths_to_track = num_lengths_to_track),
                      key = lambda C:(len(C), C.redundancy, -C.max_coverage))
    sc_solutions = []
    sc_solutions.extend(generate_sc_solutions())
    eliminated_grna = []
    ## function to generate set covers
    def make_set_cover(restore = []):
        ## restore only works if gRNA belonged in the current set cover
        curr_sc = sc_solutions[0]
        for grna in restore:
            curr_sc.add_grna(grna)
            eliminated_grna.remove(grna)
        ## if current set cover solution has at least one CollapsedgRNA with no gRNA left
        while not curr_sc.all_not_empty():
            sink = sc_solutions.pop(0) ## remove set cover solution
            ## generate more possible gRNA sets if no pre-generated set covers are left
            if not sc_solutions:
                collapsed_grnas.remove_grna(*eliminated_grna)
                collapsed_grnas.remove_empty()
                sc_solutions.extend(generate_sc_solutions())
                if not sc_solutions:
                    if not suppress_warning:
                        print(("\nError: The provided gRNA sequences cannot cover all"
                               " target sequences at least once.\n"))
                    return []
            ## select next solution
            curr_sc = sc_solutions[0]
        ## consume=True -> remove selected gRNA from CollapsedgRNA
        output = curr_sc.generate_grna_set(prioritise_3prime = prioritise_3prime, consume = True)
        eliminated_grna.extend(output)
        return output
    return make_set_cover

## tie-breaker functions for make_set_cover_pos
def all_best_nr(potential_coverage, all_coverage, covered):
    """
    Get all gRNA with equivalent non-redundnacy.
    
    This function prioritises gRNA with fewest target overlap with already covered targets.
    
    Arguments:
        potential_coverage (dict): dictionary of {'<gRNA seq>': [<list of gRNAHit obj>]}
            where gRNAHits in <list of gRNAHit obj> only contain hits to targets
            NOT already covered; AND
            where ONLY as yet unchosen gRNASeq obj are included
        all_coverage (dict): dictionary of {'<gRNA seq>': [<list of gRNAHit obj>]}
            where gRNAHits in <list of gRNAHit obj> only contain hits to all targets
            REGARDLESS of whether they've already been covered; AND
            where all gRNA's gRNASeq obj are included REGARDLESS of whether they've
            already been chosen
        covered (set): set of IDs of targets already covered
    
    Returns
    -------
    dict
        Of {'<gRNA seq>': [<list of gRNAHit obj>]} subset of 'potential_coverage'
        with equivalent non-redundancy
    """
    ## get redundancy count
    potential_redundancy = {grna_seq: len(set(hit.target_id for hit in hits
                                              if hit.target_id in covered))
                            for grna_seq, hits in all_coverage.items()
                            if grna_seq in potential_coverage}
    best_redundancy = min(potential_redundancy.values())
    return {grna_seq: potential_coverage[grna_seq]
            for grna_seq, redundancy in potential_redundancy.items()
            if redundancy == best_redundancy}

def all_best_pos(potential_coverage, all_coverage, covered):
    """
    Get all gRNA with equivalent closeness to 5'.
    
    This function prioritises gRNA closest to 5' end.
    
    Arguments:
        potential_coverage (dict): dictionary of {'<gRNA seq>': [<list of gRNAHit obj>]}
            where gRNAHits in <list of gRNAHit obj> only contain hits to targets
            NOT already covered; AND
            where ONLY as yet unchosen gRNASeq obj are included
        all_coverage (dict): dictionary of {'<gRNA seq>': [<list of gRNAHit obj>]}
            where gRNAHits in <list of gRNAHit obj> only contain hits to all targets
            REGARDLESS of whether they've already been covered; AND
            where all gRNA's gRNASeq obj are included REGARDLESS of whether they've
            already been chosen
        covered (set): set of IDs of targets already covered
    
    Returns
    -------
    dict
        Of {'<gRNA seq>': [<list of gRNAHit obj>]} subset of 'potential_coverage'
        with equivalent closeness to 5'
    """
    ## get closeness to 5'
    proximity = {grna_seq: sum((hit.target_len - hit.range[1] if
                                hit.target.sense == '-' else
                                hit.range[0])
                               for hit in hits)/len(hits)
                 for grna_seq, hits in potential_coverage.items()}
    best_proximity = min(proximity.values())
    return {grna_seq: potential_coverage[grna_seq]
            for grna_seq, prox in proximity.items()
            if prox == best_proximity}

def tie_break_first(cov, all_cov, coverage):
    """
    Arbitrarily returns the first gRNA sequence and its associated gRNAHit objects.
    
    Arguments:
        potential_coverage (dict): dictionary of {'<gRNA seq>': [<list of gRNAHit obj>]}
            where gRNAHits in <list of gRNAHit obj> only contain hits to targets
            NOT already covered; AND
            where ONLY as yet unchosen gRNASeq obj are included
        all_coverage (dict): dictionary of {'<gRNA seq>': [<list of gRNAHit obj>]}
            where gRNAHits in <list of gRNAHit obj> only contain hits to all targets
            REGARDLESS of whether they've already been covered; AND
            where all gRNA's gRNASeq obj are included REGARDLESS of whether they've
            already been chosen
        covered (set): set of IDs of targets already covered
    
    Returns
    -------
    str
        Of gRNA sequence
    list
        Of gRNAHit objects (for as yet uncovered targets) associated with the above gRNA sequence
    """
    return tuple(cov.items())[0]

## note that tie_breaker function should work on dictionaries of {gRNA_seq: {gRNAHit objects}} and return a tuple or list of two values: (gRNA_seq, {gRNAHit objects})
def make_set_cover_pos(gRNA_hits, num_sets = 1, target_ids = [], algorithm = "LAR",
                       id_key = lambda x: x, tie_breaker = tie_break_first, suppress_warning = False):
    """
    Create function that generates mutually exclusive minimum gRNA sets with position as a priority.
    
    Arguments:
        gRNA_hits (list): gRNAHit objects
        target_ids (list): list of target names/IDs (str) to be covered.
            If not provided, will be inferred from the set of target IDs covered by gRNA hits in 'gRNA_hits'
        algorithm (str): set cover algorithm
        exclude_seqs (set/list): gRNA sequences (str) to exclude
        id_key (func): function to extract target ID from gRNAHit obj
        tie_breaker (func): tie-breaker function.
            Takes (1) 'gRNA_coverage' filtered for unselected gRNA seq,
            (2) unmodified 'gRNA_coverage',
            (3) list of IDs of targets covered by already selected gRNA
        suppress_warning (bool): suppress printing of warning when gRNA cannot cover all targets
    
    Returns
    -------
    list
        Minimum set of gRNA (:class:`gRNA`)
    """
    # exclude_seqs = set(str(s).upper() for s in exclude_seqs)
    # gRNA_coverage = {seq: hits for seq, hits in gRNA_hits.hits.items()
    #                  if str(seq).upper() not in exclude_seqs}
    gRNA_coverage = gRNA_hits.hits
    eliminated_gRNA = {}
    ## prepare target ids
    if not target_ids:
        target_ids = set(hit.target_id for hit in gRNA_hits.flatten_hits())
    else:
        target_ids = set(target_ids)
    ## selected set cover algorithm
    if algorithm in ("LAR", "greedy"):
        set_cover_algo = set_cover_LAR if algorithm == "LAR" else set_cover_greedy
    else:
        raise Exception(f"Invalid algorithm name: '{algorithm}'")
    def coverage_possible():
        return set(id_key(hit) for hits in gRNA_coverage.values() for hit in hits) >= set(target_ids)
    ## function to generate set covers
    def make_set_cover(restore = []):
        for grna in restore:
            gRNA_coverage[grna.seq] = eliminated_gRNA[grna.seq]
        if not coverage_possible():
            if not suppress_warning:
                print(("\nError: The provided gRNA sequences cannot cover all target sequences"
                       " at least once.\n"))
            return []
        selected_grnas = set_cover_algo(gRNA_coverage, target_ids, id_key = id_key, tie_breaker = tie_breaker)
        ## remove selected gRNA from candidates, and covert to gRNA object
        output = []
        for grna_seq in selected_grnas:
            ## remove
            eliminated_gRNA[grna_seq] = gRNA_coverage[grna_seq]
            del gRNA_coverage[grna_seq]
            ## convert gRNA sequences to gRNA object
            grna_seq_obj = gRNA_hits.get_gRNAseq_by_seq(grna_seq)
            output.append(gRNA(grna_seq_obj.id, grna_seq_obj))
        return output
    return make_set_cover

## LAR algorithm
def set_cover_LAR(gRNA_coverage, target_ids,
                  id_key = lambda x: x, tie_breaker = tie_break_first):
    """
    Set cover algorithm LAR.
    
    Algorithm described in: Yang, Q., Nofsinger, A., Mcpeek, J., Phinney, J. and Knuesel, R. (2015). A Complete Solution to the Set Covering Problem. In International Conference on Scientific Computing (CSC) pp. 36–41
    
    Arguments:
        gRNA_coverage (dict): {'<gRNA seq>': [<list of gRNAHit obj associated w/ that gRNA seq>]}
        target_ids (list): IDs of targets to cover
        id_key (func): function to extract target ID from gRNAHit obj
        tie_breaker (func): tie-breaker function.
            Takes (1) 'gRNA_coverage' filtered for unselected gRNA seq,
            (2) unmodified 'gRNA_coverage',
            (3) list of IDs of targets covered by already selected gRNA
    
    Returns
    -------
    set
        Minimum set of gRNA sequences (str)
    """
    main_sorting_key = lambda k, uncovered_count: len(set(id_key(x) for x in uncovered_count[k]))
    result_cover, covered = {}, set()
    from copy import deepcopy
    uncovered_count = deepcopy(gRNA_coverage)
    ## get set cover
    while {len(v) for v in uncovered_count.values()} != {0}:
        sorting_key = lambda k: main_sorting_key(k, uncovered_count)
        max_val = sorting_key(max(uncovered_count, key = sorting_key))
        max_items = {seq: hits for seq, hits in uncovered_count.items()
                     if sorting_key(seq) == max_val}
        seq, coverage = tie_breaker(max_items, gRNA_coverage, covered)
        coverage = set(id_key(target) for target in coverage)
        covered |= coverage
        result_cover[seq] = coverage
        ## update coverage of unchosen seqs
        uncovered_count = {seq: set(hit for hit in hits if id_key(hit) not in covered)
                           for seq, hits in uncovered_count.items()}
    ## remove redundant sequences
    for seq_id, coverage in result_cover.items():
        coverage_remaining = set().union(*[v for k, v in result_cover.items() if k != seq_id])
        if coverage.issubset(coverage_remaining):
            del result_cover[seq_id]
    return set(result_cover.keys())

## greedy algorithm
def set_cover_greedy(gRNA_coverage, target_ids,
                     id_key = lambda x: x, tie_breaker = tie_break_first):
    """
    Greedy set cover algorithm.
    
    Arguments:
        gRNA_coverage (dict): {'<gRNA seq>': [<list of gRNAHit obj associated w/ that gRNA seq>]}
        target_ids (list): IDs of targets to cover
        id_key (func): function to extract target ID from gRNAHit obj
        tie_breaker (func): tie-breaker function.
            Takes (1) 'gRNA_coverage' filtered for unselected gRNA seq,
            (2) unmodified 'gRNA_coverage',
            (3) list of IDs of targets covered by already selected gRNA
    
    Returns
    -------
    set
        Minimum set of gRNA sequences (str)
    """
    main_sorting_key = lambda k, covered: len(set(id_key(x) for x in gRNA_coverage[k]) - set(covered))
    coverage_set = lambda k: set(id_key(x) for x in gRNA_coverage[k])
    covered, desired = set(), []
    while set(covered) != set(target_ids):
        sorting_key = lambda k: main_sorting_key(k, covered)
        max_val = sorting_key(max(gRNA_coverage, key = sorting_key))
        max_items = {k: v for k, v in gRNA_coverage.items() if sorting_key(k) == max_val}
        subset = tie_breaker(max_items, gRNA_coverage, covered)[0]
        covered |= coverage_set(subset)
        desired.append(subset)
    return set(desired)
