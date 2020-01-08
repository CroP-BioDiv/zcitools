from zci_bio.annotations.steps import AnnotationsStep
from common_utils.file_utils import write_fasta  # copy_file, link_file

_instructions = """
Open web page http://www.herbalgenomics.org/cpgavas/
Probably one of mirrors:
 Mirror 1: Central China  : http://47.96.249.172:16019/analyzer/home
 Mirror 2: East Coast USA : http://47.90.241.85:16019/analyzer/home  (more stable)

For each sequence (fas file) do:
 * Upload file: sequence.fas
 * Specify project name, species name if needed, and email address for notification.
 * Leave other data on default
 * Submit job
 * When job is finished:
 - download Global multi-GenBank file into job directory ({abspath})
 - run zcit command: zcit.py finish {step_name}


The paper describing CPGAVAS2 can be found here:
https://academic.oup.com/nar/advance-article/doi/10.1093/nar/gkz345/5486746

"""


def create_cpgavas_data(step_data, sequences_step):
    step = AnnotationsStep(sequences_step.project, step_data, remove_data=True)

    # Store sequence
    for seq_ident in sequences_step.all_sequences():
        seq = sequences_step.get_sequence(seq_ident)
        seq = seq.replace('N', '')
        # ToDo: napraviti mapiranje
        write_fasta(step.step_file(seq_ident + '.fas'), [(seq_ident, seq)])

    # Store instructions
    with open(step.step_file('INSTRUCTIONS.txt'), 'w') as out:
        out.write(_instructions.format(abspath=step.absolute_path(), step_name=step_data['step_name']))

    #
    step.set_sequences(sequences_step.all_sequences())
    step.save(completed=False)
    return step


def finish_cpgavas_data(step_obj):
    prnt("ToDo: ...")
    # # Check file named: GeSeqJob-<num>-<num>_GLOBAL_multi-GenBank.gbff
    # for f in step_obj.step_files():
    #     if f.startswith('GeSeqJob') and f.endswith('_GLOBAL_multi-GenBank.gbff'):
    #         filename = f
    #         break
    # else:
    #     print("Warning: can't find GeSeq output file!")
    #     return

    # # Leave original file
    # # ToDo: repair and filter data???
    # # ToDo: inverted_region 126081..1 !!! To_ind > from_ind!!!
    # copy_file(step_obj.step_file(filename), step_obj.get_all_annotation_filename())
    # step_obj._check_data()
    # step_obj.save()
