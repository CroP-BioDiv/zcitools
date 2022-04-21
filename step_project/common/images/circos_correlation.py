import os.path
import subprocess
from .steps import ImagesStep
from common_utils.exceptions import ZCItoolsValueError
from common_utils.data_types.correlation_matrix import CorrelationMatrix
from common_utils.file_utils import ensure_directory, write_str_in_file, get_settings

_circos_conf = """
<colors>
{colors}
</colors>

# Groups
karyotype = data/karyotype.txt

<ideogram>
<spacing>
default = 0.020r
</spacing>

thickness        = 40p
stroke_thickness = 0
stroke_color     = vdgrey
fill             = yes
fill_color       = black

# fractional radius position of chromosome ideogram within image
radius         = 0.90r
show_label     = yes
label_font     = bold
label_radius   = dims(image,radius) - 100p
label_size     = 50
label_parallel = yes

show_bands     = no
</ideogram>

# 1-correlation group parts
<highlights>
z          = 0
<highlight>
file       = data/tiles.txt
r0         = 0.999r-30p
r1         = 0.999r-5p
stroke_thickness = 0
</highlight>
</highlights>

# Correlations
<links>
<link>
ribbon           = yes
flat             = yes
file             = data/links.txt
bezier_radius    = 0.0r
radius           = 0.999r-30p
thickness        = 10
color            = grey
stroke_color     = dgrey
stroke_thickness = 1

<rules>
<rule>
condition  = var(dist) <= 1.5
bezier_radius    = 0.3r
</rule>
</rules>
</link>

</links>

<image>
<<include etc/image.conf>>
</image>

<<include etc/colors_fonts_patterns.conf>>
<<include etc/housekeeping.conf>>
"""


def create_circos_correlation(project, step_data, params):
    # Read correlation data
    cm = None
    if params.input_filename:
        cm = CorrelationMatrix.from_file(params.input_filename)

    if not cm:
        raise ZCItoolsValueError('No correlation input data!')
    num_c = cm.num_columns()
    if num_c < 2:
        raise ZCItoolsValueError('Not much of a matrix!')

    step = ImagesStep(project, step_data, remove_data=True)
    one_width = params.one_width
    gap_correlations = params.gap_correlations
    ow_2 = one_width // 2
    one_plus_gap = one_width + gap_correlations

    # Note: column lowercase names are used as column identifiers
    data_dir = step.step_file('data')
    etc_dir = step.step_file('etc')
    ensure_directory(data_dir)
    ensure_directory(etc_dir)

    colors = dict((lc, 'green') for lc in cm._columns_lower)  # ToDo: some defaults
    colors['plus_'] = 'blue'
    colors['minus_'] = 'red'
    for col_def in params.group_color or []:
        col_fields = col_def.split(',', 1)
        if len(col_fields) == 2 and cm.check_column(col_fields[0]):
            colors[cm.check_column(col_fields[0])] = col_fields[1]
        else:
            print(f"Warning: '{col_def}' is not column color definition!")

    # data directory
    # karyotype.txt: defines groups (as chromosomes)
    # chr - <name> <label> <start> <end> <color>
    # ...
    gl = (num_c - 1) * one_width + (num_c - 2) * gap_correlations  # group length
    write_str_in_file(os.path.join(data_dir, 'karyotype.txt'),
                      '\n'.join(f"chr - {lc} {c} 0 {gl} color_{lc}"
                                for lc, c in zip(cm._columns_lower, cm._columns)))

    # tiles.txt: defines abs(correlation) == 1 interval, as tiles
    # <name> <start> <end> [options]
    with open(os.path.join(data_dir, 'tiles.txt'), 'w') as out:
        for idx1, c1 in enumerate(cm._columns_lower):
            for idx2, c2 in enumerate(cm._columns_lower):
                if idx1 == idx2:
                    continue
                pos = (idx1 - idx2 - 1) if idx1 > idx2 else (idx1 - idx2 + (num_c - 1))
                start = pos * one_plus_gap
                out.write(f"{c1} {start} {start + one_width} fill_color=color_{c2}\n")

    # cells.txt: defines correlations as links
    # <cell_idx> <group_1> <start_1> <end_1> color=color_{plus|minus}_,dist={int}
    # <cell_idx> <group_2> <start_2> <end_2> color=color_{plus|minus}_,dist={int}
    # ...
    with open(os.path.join(data_dir, 'links.txt'), 'w') as out:
        cell_idx = 0
        for idx1, c1 in enumerate(cm._columns_lower):
            rest_c = cm._columns_lower[idx1 + 1:]
            for idx2, c2 in enumerate(rest_c):
                corr = cm.get(c1, c2)
                if corr is not None:
                    w = round(abs(corr) * one_width)
                    w_1 = w // 2
                    w_2 = w - w_1  # - 1?
                    centar = ow_2 + idx2 * one_plus_gap
                    color = 'plus_' if corr >= 0 else 'minus_'
                    dist = min(idx2 + 1, idx1 + (len(rest_c) - idx2))
                    atts = f"color=color_{color},dist={dist}"
                    out.write(f"cell_{cell_idx} {c1} {gl - centar - w_2} {gl - centar + w_1} {atts}\n")
                    out.write(f"cell_{cell_idx} {c2} {centar - w_1} {centar + w_2} {atts}\n")
                    cell_idx += 1

    # etc directory
    write_str_in_file(os.path.join(etc_dir, 'circos.conf'), _circos_conf.format(
        colors='\n'.join(f"color_{lc} = {c}" for lc, c in colors.items())
    ))

    subprocess.run(['circos', '-conf', 'etc/circos.conf'], cwd=step.directory)

    # View it
    if params.show_image:
        image_viewer = get_settings().get('image_viewer')
        if image_viewer:
            subprocess.Popen([image_viewer, step.step_file('circos.png')])

    #
    # # step.set_table_data(data, columns)
    # step.save()
    # return step
