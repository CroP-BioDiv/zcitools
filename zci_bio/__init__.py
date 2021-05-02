from common_utils.misc import find_registered
from .workflows.phylogenetic_analysis import PhylogeneticAnalysis
from .workflows.chloroplast_normalization import ChloroplastNormalization
from .workflows.irs_statistics import IRsStatistics
registered_commands, registered_steps = find_registered(__file__, 'zci_bio', exclude_dirs=['workflows'])
registered_workflows = [PhylogeneticAnalysis, ChloroplastNormalization, IRsStatistics]
