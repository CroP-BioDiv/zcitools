from common_utils.misc import find_registered
from .workflows.chloroplast_normalization import ChloroplastNormalization
registered_commands, registered_steps = find_registered(__file__, 'zci_bio', exclude_dirs=['workflows'])
registered_workflows = [ChloroplastNormalization]
