import os
from ..steps import read_step
from ..utils.import_methods import import_pygraphviz


def create_graph():
    pygraphviz = import_pygraphviz()
    graph = pygraphviz.AGraph(strict=True, directed=True)

    # Nodes
    edges = []
    for d in sorted(os.listdir('.')):
        if os.path.isdir(d) and os.path.isfile(os.path.join(d, 'description.yml')):
            step = read_step(d)
            node = step.directory
            graph.add_node(node)
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
    graph.draw(output_filename + '.ps', prog='dot')

    # ?Show
