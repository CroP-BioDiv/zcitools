# Step classes. Used for filling _type_2_step_cls dict
from .images import ImagesStep

registered_steps = [cls for cls in locals().values() if getattr(cls, '_STEP_TYPE', None)]
