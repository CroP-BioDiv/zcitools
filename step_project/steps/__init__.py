# Step classes. Used for filling _type_2_step_cls dict
from .table import TableStep
from .images import ImagesStep
from .alignments import AlignmentStep, AlignmentsStep
from .phylogenetics import RAxMLStep, RAxMLSteps, MrBayesStep, MrBayesSteps

registered_steps = [cls for cls in locals().values() if getattr(cls, '_STEP_TYPE', None)]
