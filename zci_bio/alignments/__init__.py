from .commands import ClustalO, mVISTA
from .steps import AlignmentStep, AlignmentsStep, mVISTAStep

registered_commands = [ClustalO, mVISTA]
registered_steps = [AlignmentStep, AlignmentsStep, mVISTAStep]
