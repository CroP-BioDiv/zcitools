from zcitools.utils.import_methods import import_CAI

# Check: https://github.com/Benjamin-Lee/CodonAdaptationIndex


def calcualte_rscu(step, calc_filename):
    RSCU = import_CAI('RSCU')
    step_idents, cds_seqs = zip(*step.extract_all_cds())
    r = RSCU(cds_seqs)
    print(r)
