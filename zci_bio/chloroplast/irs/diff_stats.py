#!/usr/bin/env python3

from common_utils.properties_db import PropertiesDB

# _shorted_char = dict(equal='=', replace='*', insert='+', delete='-')

# Sto moze biti problematicno, upitno:
#  - Suljina IR-ova. Ako je <20k, onda nes ne stima
#  - Predug data['diff']
#  - Cudan data['diff']
#  - Velika razlika izmedju duljina IR-ova po alatima

# Sto znam da moze biti problem:
#  - male rupe, neznam koliko odvojene; Chloroplot, Blast metode
#  - velike razlike; GeSeq jer trzi krajeve

# Stat:
#  - max indel, max replace
#  - min udaljenost izmedju razlika


def diff_stats(project, args):
    seq_idents = args.seq_ident if args.seq_ident else []
    if args.table_step:
        step = project.read_step(args.table_step, check_data_type='table')
        seq_idents.extend(step.get_column_values_by_type('seq_ident'))
    if args.sequences_step:
        step = project.read_step(args.sequences_step, check_data_type='sequences')
        seq_idents.extend(step.all_sequences())
    #
    if args.first_n:
        seq_idents = seq_idents[:args.first_n]
    #
    diff_stats_idents(seq_idents, args.merge_distance)


def diff_stats_idents(seq_idents, merge_distance):
    properties_db = PropertiesDB()
    # ncbi?
    methods_ann = [f'annotation {a}' for a in ('airpg', 'chloe', 'chloroplot', 'ge_seq', 'org_annotate', 'pga', 'plann')]
    for seq_ident in sorted(seq_idents):
        print(seq_ident)
        annots = properties_db.get_properties_keys2(seq_ident, methods_ann)
        max_data = None
        # Find the best one
        for method, data in annots.items():
            method = method.split()[1]
            if 'ira' in data:
                # Skip problematic IRs
                if data.get('type') == '???':
                    continue
                if method in ('airpg', 'ncbi') and len(data.get('diff', [])) > 1000:
                    print(f'Warning: sequence {seq_ident}, {method} diff length {len(data["diff"])}!')
                    continue
                if len(data.get('diff', [])) > 100:
                    print(f'Info: sequence {seq_ident}, {method} diff length {len(data["diff"])}!')

                # Probati po redu pa sredjivati jednu po jednu!!!

                if 'diff' in data:
                    analyse_diff(data)  # , merge_distance)
                # if not max_data or _is_better(data, max_data):
                #     max_data = data
        # print(seq_ident, list(annots.keys()))


def analyse_diff(data):
    diff = data['diff']
    if len(diff[0]) == 4:
        # ToDo: prebaciti u diff module ili vidjeti zasto uopce postoji!!!
        # First element can be tuple of type ('equal', length, 0, 0)?!?!
        d0 = diff[0]
        assert d0[0] == 'equal', d0
        assert d0[2] == d0[3] == 0, d0
        diff[0] = ('equal', 0, d0[1], 0, d0[1])

    diff = [('I' if d[0] in ('insert', 'delete') else d[0][0].upper(), max(d[2] - d[1], d[4] - d[3])) for d in diff]
    equal_lengths = sorted(d[1] for d in diff if d[0] == 'E')

    return dict(max_indel=max((d[1] for d in diff if d[0] == 'I'), default=0),
                max_replace=max((d[1] for d in diff if d[0] == 'R'), default=0),
                equal_lengths=equal_lengths,
                diff_simple=diff)


def _is_better(data, max_data):
    data['ir_lengths']
    pass


if __name__ == '__main__':
    import sys
    diff_stats_idents(sys.argv[2:], int(sys.argv[1]))

    # import argparse
    # parser = argparse.ArgumentParser(description="Mounts file in specified directory")
    # parser.add_argument('-i', '--seq-ident', action='append', help='Sequence identifier')
    # params = parser.parse_args()

    # diff_stats_idents(params.seq_ident)
