import os
from ..utils.import_methods import import_pygraphviz
from ..utils.file_utils import get_settings


def create_graph(zcit):
    pygraphviz = import_pygraphviz()
    graph = pygraphviz.AGraph(strict=True, directed=True)
    # graph.graph_attr['rankdir'] = 'LR'  # Orientation left-right
    # graph.node_attr['shape'] = 'plaintext'
    # graph.node_attr['shape'] = 'record'
    # graph.edge_attr['lblstyle'] = 'above, sloped'

    # Nodes
    edges = []
    for d in sorted(os.listdir('.')):
        if os.path.isdir(d) and os.path.isfile(os.path.join(d, 'description.yml')):
            step = zcit.read_step(d)
            node = step.directory
            label = node if step.is_completed() else '* ' + node
            graph.add_node(node, label=label)
            #
            prev_steps = step._step_data['prev_steps']
            if prev_steps:
                command = step._step_data['command']
                edges.extend((p, node, command) for p in prev_steps)

    # Edges
    for from_node, to_node, label in edges:
        graph.add_edge(from_node, to_node, label=label)

    # Export
    output_filename = 'graph'
    # graph.write(output_filename + '.dot')
    # prog=neato|dot|twopi|circo|fdp|nop
    ps_file = output_filename + '.ps'
    graph.draw(ps_file, prog='dot')

    ps_viewer = get_settings()['ps_viewer']
    if ps_viewer:
        os.system(f"{ps_viewer} {ps_file}")
