from step_project.base_commands import CreateStepFromStepCommand
from .standardization.commands import *  # Import commands


class ChloroplastAlign(CreateStepFromStepCommand):
    _COMMAND = 'align_chloroplast'
    _HELP = "Align chloroplast genomes, by specifics."
    _COMMAND_GROUP = 'Chloroplast'
    _INPUT_STEP_DATA_TYPE = ('sequences', 'annotations')  # If there are not IRs, at least hwole genomes will be aligned

    @staticmethod
    def set_arguments(parser):
        from ..alignments.common_methods import ALIGN_PROGRAMS
        CreateStepFromStepCommand.set_arguments(parser)
        parser.add_argument('sequences', nargs='+', help='Sequence identifiers (at least two)')
        parser.add_argument('-p', '--alignment-program', default='mafft', choices=ALIGN_PROGRAMS,
                            help='Alignmnet program to use.')
        parser.add_argument('-a', '--to-align', default='whole',
                            help='What to align. Options: whole, parts, trnH-GUG, offset. ' +
                                 'Format: "opt1:opt2:...". Only first character is needed.')
        parser.add_argument('-r', '--run', action='store_true', help='Run locally')
        parser.add_argument('-o', '--keep-offset', type=int, help='Offset to keep.')

    def _to_align(self):
        return sorted(x[0].lower() for x in self.args.to_align.split(';'))

    def step_base_name(self):
        return f'AlignChloroplasts_{"".join(self._to_align())}_{"_".join(sorted(self.args.sequences))}'

    def run(self, step_data):
        from .utils import chloroplast_alignment
        args = self.args
        return chloroplast_alignment(
            step_data, self._input_step(no_data_check=True),
            args.sequences, self._to_align(), args.run, args.alignment_program, args.keep_offset)
