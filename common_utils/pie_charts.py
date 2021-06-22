import math
from collections import namedtuple
from .import_method import import_matplotlib_pyplot


class Value(namedtuple('Value', 'label, value, color')):
    def __new__(cls, label, value, color=None):
        return super(Value, cls).__new__(cls, label, value, color)


def _values_2_lists(values, normalize=True, revert=False):
    labels = [v.label for v in values if v.value]
    ratios = [v.value for v in values if v.value]
    colors = [v.color for v in values if v.value and v.color]
    if revert:
        labels = labels[::-1]
        ratios = ratios[::-1]
        colors = colors[::-1]
    #
    if len(colors) != len(ratios):
        if colors:
            print('Warning: not all colors are specified!')
        colors = None

    if normalize:
        s_ratios = sum(ratios)
        if abs(s_ratios - 1) > 0.01:
            ratios = [x / s_ratios for x in ratios]
    return ratios, labels, colors


def _save_chart(pyplot, output_filename, dpi=150):
    if output_filename:
        pyplot.savefig(f'{output_filename}.svg')
        pyplot.savefig(f'{output_filename}.png', dpi=dpi, pad_inches=0)
        pyplot.savefig(f'{output_filename}.eps', format='eps', dpi=dpi)

    pyplot.show()


def pie_chart_bar_of_pie(pie_data, explode=None, output_filename=None):
    # https://matplotlib.org/stable/gallery/pie_and_polar_charts/bar_of_pie.html
    # pie_data is dict with attrs:
    #  - pie : list _Value objects
    #  - bar : list _Value objects
    #  - title : string
    pyplot = import_matplotlib_pyplot(use_arial=True)
    pyplot.rc('font', size=5)          # controls default text sizes
    pyplot.rc('legend', fontsize=5)    # legend fontsize

    # make figure and assign axis objects
    fig, (ax1, ax2) = pyplot.subplots(1, 2, figsize=(2.74, 1.52), dpi=450)
    fig.subplots_adjust(wspace=0, left=0.08, right=.89, top=1., bottom=0.05)

    _pie_chart_bar_of_pie(ax1, ax2, pie_data, explode=explode)
    _save_chart(pyplot, output_filename, dpi=450)


def pie_charts_bar_of_pie(pies_data, explode=None, output_filename=None):
    pyplot = import_matplotlib_pyplot(use_arial=True)
    num_pies = len(pies_data)
    pyplot.rc('font', size=16)          # controls default text sizes
    pyplot.rc('legend', fontsize=16)    # legend fontsize

    # make figure and assign axis objects
    # fig, axis = pyplot.subplots(1, 3 * num_pies - 1, figsize=(9 * 2 * num_pies, 5))
    # fig.subplots_adjust(wspace=0, left=0.05, right=1)  # , left=-1, right=1, top=0.9, bottom=0.05)
    # for i in range(1, num_pies):
    #     axis[i * 3 - 1].axis('off')
    fig, axis = pyplot.subplots(1, 2 * num_pies, figsize=(9 * 2 * num_pies, 5))
    fig.subplots_adjust(wspace=0, left=0.02, right=1.02, bottom=0.05, top=0.92)  # , top=1.02)

    for i, pie_data in enumerate(pies_data):
        # _pie_chart_bar_of_pie(axis[i * 3], axis[i * 3 + 1], pie_data, explode=explode)
        _pie_chart_bar_of_pie(axis[i * 2], axis[i * 2 + 1], pie_data, explode=explode)

    _save_chart(pyplot, output_filename)


def _pie_chart_bar_of_pie(ax1, ax2, pie_data, explode=None):
    from matplotlib.patches import ConnectionPatch

    # pie chart parameters
    pie_ratios, pie_labels, pie_colors = _values_2_lists(pie_data['pie'])
    num_parts = len(pie_ratios)

    # Explode
    if explode:
        _explode = [0.1 if explode is True else explode] + [0] * (num_parts - 1)
    else:
        _explode = None

    # rotate so that first wedge is split by the x-axis
    angle = -180 * pie_ratios[0]
    # wedges, texts, autotexts =
    ax1.pie(pie_ratios, autopct='%1.1f%%', startangle=angle, labels=pie_labels, explode=_explode, colors=pie_colors)
    if title := pie_data.get('title'):
        ax1.set_title(title, loc='left')

    # bar chart parameters
    xpos = 0
    bottom = 0
    bar_ratios, bar_labels, bar_colors = _values_2_lists(pie_data['bar'], revert=True)
    num_bars = len(bar_ratios)
    bar_width = pie_data.get('bar_width', .2)
    if not bar_colors:
        d = (0.9 - 0.3) / (num_bars - 1)
        bar_colors = [[0.1, 0.3, 0.3 + d * x] for x in range(num_bars)]

    x_label = bar_width * 0.6
    for j, height in enumerate(bar_ratios):
        ax2.bar(xpos, height, bar_width, bottom=bottom, color=bar_colors[j])
        ypos = bottom + ax2.patches[j].get_height() / 2
        bottom += height
        ax2.text(x_label, ypos, bar_labels[j], ha='left', va='center')

    # ax2.legend(('50-65', 'Over 65', '35-49', 'Under 35'))
    ax2.axis('off')
    ax2.set_xlim(- 2.5 * bar_width, 2.5 * bar_width)

    # use ConnectionPatch to draw lines between the two plots
    # get the wedge data
    theta1, theta2 = ax1.patches[0].theta1, ax1.patches[0].theta2
    center, r = ax1.patches[0].center, ax1.patches[0].r
    bar_height = sum([item.get_height() for item in ax2.patches])

    # draw top connecting line
    x = r * math.cos(math.pi / 180 * theta2) + center[0]
    y = r * math.sin(math.pi / 180 * theta2) + center[1]
    con = ConnectionPatch(xyA=(-bar_width / 2, bar_height), coordsA=ax2.transData,
                          xyB=(x, y), coordsB=ax1.transData)
    con.set_color('black')
    con.set_linewidth(0.25)
    con.set_linestyle((0, (1, 4)))
    ax2.add_artist(con)

    # draw bottom connecting line
    x = r * math.cos(math.pi / 180 * theta1) + center[0]
    y = r * math.sin(math.pi / 180 * theta1) + center[1]
    con = ConnectionPatch(xyA=(-bar_width / 2, 0), coordsA=ax2.transData,
                          xyB=(x, y), coordsB=ax1.transData)
    con.set_color('black')
    ax2.add_artist(con)
    con.set_linewidth(0.25)
    con.set_linestyle((0, (1, 4)))
