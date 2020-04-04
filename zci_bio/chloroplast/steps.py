from step_project.base_step import SimpleStep


class ChloroplastSSCBlast(SimpleStep):
    "Stores data for Blasting SSC ends in sequences"
    _STEP_TYPE = 'chloroplast_ssc_blast'


class ChloroplastOrientateStep(SimpleStep):
    "Stores info about chloroplast orientation"
    _STEP_TYPE = 'chloroplast_orientate'
