import os
from step_project.utils.import_methods import import_pygraphviz
from common_utils.file_utils import get_settings


def create_project_graph(project, output_filename='project_graph'):
    # Nodes
    nodes = []
    edges = []
    for d in sorted(os.listdir('.')):
        if os.path.isdir(d) and os.path.isfile(os.path.join(d, 'description.yml')):
            step = project.read_step(d, no_check=True)
            node = step.directory
            label = node if step.is_completed() else '* ' + node
            nodes.append((node, dict(label=label)))
            #
            prev_steps = step._step_data['prev_steps']
            if prev_steps:
                command = step._step_data['command']
                edges.extend((p, node, dict(label=command)) for p in prev_steps)

    create_graph_from_data(nodes, edges, output_filename)


def create_graph_from_data(nodes, edges, output_filename):
    pygraphviz = import_pygraphviz()
    attrs = dict(strict=True, directed=True)
    if len(nodes) > 10:
        attrs['rankdir'] = 'LR'
    graph = pygraphviz.AGraph(**attrs)
    # graph.graph_attr['rankdir'] = 'LR'  # Orientation left-right
    # graph.node_attr['shape'] = 'plaintext'
    # graph.node_attr['shape'] = 'record'
    # graph.edge_attr['lblstyle'] = 'above, sloped'

    for node, attrs in nodes:
        graph.add_node(node, **attrs)

    # Edges
    for from_node, to_node, attrs in edges:
        graph.add_edge(from_node, to_node, **attrs)

    # Export
    graph.write(output_filename + '.dot')
    # prog=neato|dot|twopi|circo|fdp|nop
    ps_file = output_filename + '.ps'
    graph.draw(ps_file, prog='dot')

    ps_viewer = get_settings()['ps_viewer']
    if ps_viewer:
        os.system(f"{ps_viewer} {ps_file}")
