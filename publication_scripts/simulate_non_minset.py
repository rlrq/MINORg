import sys
sys.path.append("/mnt/chaelab/rachelle/src")

import itertools
from data_manip import splitlines, parse_get_data

proj_dir = "/mnt/chaelab/rachelle/minorg/"

f_dat_map = proj_dir + "/data/TN3_vdw_nbs_cas12a_23nt/TN3_vdw_gRNA_all.wo1925_2.map"

get, dat = parse_get_data(f_dat_map)

set_pos_wo_replacement = prioritise_pos_wo_replacement(dat, get, 5)
set_pos_w_replacement = prioritise_pos_w_replacement(dat, get, 5)
# set_pos_w_replacement_overlap = report_set_overlap(set_pos_w_replacement)
reps = 100
sets_random_wo_replacement = [random_grna_wo_replacement(dat, get, 5) for i in range(reps)]
# sets_random_wo_replacement_size = [{k: len(v) for k, v in rep.items()} for rep in sets_random_wo_replacement]
sets_random_w_replacement = [random_grna_w_replacement(dat, get, 5) for i in range(reps)]
# sets_random_w_replacement_overlap = [report_set_overlap(rep) for rep in sets_random_w_replacement]
# sets_random_w_replacement_size = [{k: len(v) for k, v in rep.items()} for rep in sets_random_w_replacement]

fout_sim = proj_dir + "/results/TN3_vdw_nbs_cas12a_23nt_simulate_non_minset_wcoverageinfo_rep100.tsv"
with open(fout_sim, "w+") as f:
    f.write('\t'.join(["simulation", "replicate", "group", "coverage",
                       "num gRNA", "shared gRNA", "gRNA id"]) + '\n')
    simulations = {"5prime_w_replacement": [set_pos_w_replacement],
                   "5prime_wo_replacement": [set_pos_wo_replacement],
                   "random_w_replacement": sets_random_w_replacement,
                   "random_wo_replacement": sets_random_wo_replacement}
    for simulation, replicates in simulations.items():
        for i, replicate in enumerate(replicates):
            rep, cov = replicate
            group_overlap = report_set_overlap(rep)
            for group in sorted(rep.keys()):
                members = rep[group]
                coverage = cov[group]
                f.write('\t'.join(map(str, [simulation, i, group, coverage, len(members),
                                            group_overlap[group], ','.join(sorted(members))])) + '\n')


def group_by_field(dat, get, field):
    output = {}
    for entry in dat:
        field_value = get(entry, field)
        output[field_value] = output.get(field_value, []) + [entry]
    return output

def prioritise_pos_wo_replacement(dat, get, sets = 1):
    selected_grna = {}
    selected_grna_coverage = {}
    used_grna = set()
    ## filter for passing gRNA
    dat = [entry for entry in dat if "fail" not in entry[8:]]
    d_dat = group_by_field(dat, get, "target id")
    targets = set(d_dat.keys())
    ## sort grna by position (per target)
    sorted_dat = {target: (sorted(entries, key = lambda entry: int(get(entry, "start")))
                           if get(entries[0], "target sense") == "sense"
                           else sorted(entries, key = lambda entry: -int(get(entry, "end"))))
                  for target, entries in d_dat.items()}
    ## generate groups
    for group in range(1, sets+1):
        grna_in_group = set()
        group_coverage = set()
        for target, entries in sorted_dat.items():
            for i, entry in enumerate(entries):
                grna_id = get(entry, "gRNA id")
                if grna_id not in used_grna:
                    grna_in_group = grna_in_group.union({grna_id})
                    group_coverage = group_coverage.union({target})
                    sorted_dat[target] = entries[i+1:]
                    break
        selected_grna[group] = grna_in_group
        selected_grna_coverage[group] = group_coverage
        used_grna = used_grna.union(grna_in_group)
    ## print warning if not all targets have enough gRNA
    groups_wo_enough_coverage = {k: len(v) for k, v in selected_grna_coverage.items()
                                 if v != targets}
    if groups_wo_enough_coverage:
        print(f"The following groups cannot cover all targets.")
        print("group\tnumber of targets")
        for group, num_targets in groups_wo_enough_coverage.items():
            print(f"{group}\t{num_targets}")
    return selected_grna, {k: len(cov) for k, cov in selected_grna_coverage.items()}

def prioritise_pos_w_replacement(dat, get, sets = 1):
    selected_grna = {}
    ## filter for passing gRNA
    dat = [entry for entry in dat if "fail" not in entry[8:]]
    d_dat = group_by_field(dat, get, "target id")
    targets_wo_enough_grna = {} ## target_id: number of gRNA
    ## get top n gRNA (from 5' end)
    for target_id, entries in d_dat.items():
        if len(entries) < sets:
            targets_wo_enough_grna[target_id] = len(entries)
        sense = get(entries[0], "target sense") == "sense"
        if sense:
            selected_grna[target_id] = sorted(entries, key = lambda entry: int(get(entry, "start")))[:sets]
        else:
            selected_grna[target_id] = sorted(entries, key = lambda entry: -int(get(entry, "end")))[:sets]
    ## make "sets"
    output = {i+1: set(grna_id for grna_id in
                       ((None if len(v) < i+1 else get(v[i], "gRNA id"))
                        for v in selected_grna.values())
                       if grna_id is not None) for i in range(sets)}
    ## print warning if not all targets have enough gRNA
    if targets_wo_enough_grna:
        print(f"The following targets did not have enough gRNA to populate {sets} sets.")
        print("target id\tnumber of gRNA")
        for target_id, num_grna in targets_wo_enough_grna.items():
            print(f"{target_id}\t{num_grna}")
    ## get set coverage
    d_dat_bygrna = group_by_field(dat, get, "gRNA id")
    coverage = {s: len(set(itertools.chain(*[get(d_dat_bygrna[grna_id], "target id", return_list = True)
                                             for grna_id in grnas]))) for s, grnas in output.items()}
    return output, coverage

def random_grna_wo_replacement(dat, get, sets = 1):
    selected_grna = {}
    selected_grna_coverage = {}
    dat = [entry for entry in dat if "fail" not in entry[8:]]
    d_dat = group_by_field(dat, get, "gRNA id")
    targets = set(group_by_field(dat, get, "target id").keys())
    ## randomly order gRNA, and we'll group them according to non-overlapping N from start that cover all targets
    import random
    shuffled = list(d_dat.keys())
    random.shuffle(shuffled)
    ## start assigning groups
    group = 1
    for grna in shuffled:
        entries = d_dat[grna]
        coverage = set(get(entries, "target id"))
        selected_grna[group] = selected_grna.get(group, []) + [grna]
        selected_grna_coverage[group] = selected_grna_coverage.get(group, set()).union(coverage)
        if selected_grna_coverage[group] == targets:
            group += 1
            if group > sets: break
    ## determine if target number of sets reached
    if len(selected_grna) < sets:
        print(f"Unable to generate {sets} sets.")
    ## determine if all groups can cover all targets
    groups_wo_complete_coverage = {k: len(v) for k, v in selected_grna_coverage.items()
                                   if v != targets}
    if groups_wo_complete_coverage:
        print("The following group(s) could not cover all target:")
        print(f"group\tnumber of targets")
        for group, num_targets in groups_wo_complete_coverage.items():
            print(f"{group}\t{num_targets}")
    return selected_grna, {k: len(cov) for k, cov in selected_grna_coverage.items()}

def random_grna_w_replacement(dat, get, sets = 1):
    selected_grna = {}
    selected_grna_coverage = {}
    dat = [entry for entry in dat if "fail" not in entry[8:]]
    d_dat = group_by_field(dat, get, "gRNA id")
    targets = set(group_by_field(dat, get, "target id").keys())
    ## randomly order gRNA, and we'll group them according to non-overlapping N from start that cover all targets
    import random
    for group in range(1, sets+1):
        shuffled = list(d_dat.keys())
        random.shuffle(shuffled)
        for grna in shuffled:
            entries = d_dat[grna]
            coverage = set(get(entries, "target id"))
            selected_grna[group] = selected_grna.get(group, []) + [grna]
            selected_grna_coverage[group] = selected_grna_coverage.get(group, set()).union(coverage)
            if selected_grna_coverage[group] == targets: break
    ## determine if all groups can cover all targets
    groups_wo_complete_coverage = {k: len(v) for k, v in selected_grna_coverage.items()
                                   if v != targets}
    if groups_wo_complete_coverage:
        print("The following group(s) could not cover all target:")
        print(f"group\tnumber of targets")
        for group, num_targets in groups_wo_complete_coverage.items():
            print(f"{group}\t{num_targets}")
    return selected_grna, {k: len(cov) for k, cov in selected_grna_coverage.items()}

def report_set_overlap(sets):
    output = {}
    for group, members in sets.items():
        output[group] = len([member for member in members
                             if sum(member in v for v in sets.values()) > 1])
    return output
