import os
import itertools
# import shlex
from collections import namedtuple
from common_utils.exceptions import ZCItoolsValueError
from common_utils.cache import cache
from common_utils.file_utils import remove_directory, merge_zip_files, write_str_in_file
from common_utils.value_data_types import sheets_2_excel
from .common.graph.project_graph import create_graph_from_data

_have_run_switch = ('align_genomes', 'mr_bayes', 'raxml')


class WfAction(namedtuple('WfAction', 'step_name, prev_steps, cmd, additional_reqs')):
    def __new__(cls, step_name, prev_steps, cmd, additional_reqs):
        # cmd is a string
        assert isinstance(step_name, str), step_name

        # prev_steps is None, str, or list of string
        assert isinstance(prev_steps, (list, tuple)), prev_steps
        assert all(isinstance(s, str) for s in prev_steps), prev_steps
        assert isinstance(additional_reqs, (list, tuple)), additional_reqs
        assert all(isinstance(s, str) for s in additional_reqs), additional_reqs

        # args is a list
        assert cmd and isinstance(cmd, (list, tuple)), cmd
        assert all(isinstance(s, str) for s in cmd), cmd

        return super(WfAction, cls).__new__(cls, step_name, prev_steps, cmd, additional_reqs)

    command = property(lambda self: self.cmd[0])
    all_prev_steps = property(lambda self: self.prev_steps + self.additional_reqs)

    def has_run_switch(self):
        return self.command in _have_run_switch


class BaseWorkflow:
    _WORKFLOW = None  # Name of workflow
    _COMMAND_METHODS = dict(run='cmd_run', graph='cmd_graph', summary='cmd_summary', actions='cmd_actions')

    def __init__(self, project, parameters):
        self.project = project
        self.parameters = parameters

    @classmethod
    def all_commands(cls):
        return list(cls._COMMAND_METHODS.keys())

    @staticmethod
    def required_parameters():
        raise NotImplementedError('')

    @staticmethod
    def format_parameters(params):
        return params

    @cache
    def actions(self):
        # Returns list of WfAction objects
        actions = []
        for x in self._actions():
            assert 2 <= len(x) <= 3, x
            sn, cmd = x[:2]
            a_prev = [] if len(x) == 2 else ([x[2]] if isinstance(x[2], str) else x[2])
            actions.append((sn, (cmd.split() if isinstance(cmd, str) else cmd), a_prev))
            # actions.append((sn, (shlex.split(cmd) if isinstance(cmd, str) else cmd), a_prev))

        assert all(isinstance(cmd, (list, tuple)) for _, cmd, _ in actions), actions
        assert all(isinstance(prevs, (list, tuple)) for _, _, prevs in actions), actions

        # Check commands
        commands_map = self.project.commands_map
        if (not_in := [cmd[0] for _, cmd, _ in actions if cmd[0] not in commands_map]):
            raise ZCItoolsValueError(f"Worflow actions, not existing command(s)! {', '.join(sorted(not_in))}")

        # Check step names
        step_names = set(sn for sn, _, _ in actions)
        return [WfAction(sn, [c for c in cmd if c in step_names], cmd, a_prevs) for sn, cmd, a_prevs in actions]

    def _actions(self):
        # List of tuples (step_name, cmd, additional_required_steps)
        raise NotImplementedError('')

    def get_summary(self):
        # Returns dict with data to store. Attributes can be: text, table
        raise NotImplementedError('')

    #
    @cache
    def all_step_names(self):
        return sorted(a.step_name for a in self.actions())

    def step_status(self, step_name):
        if not os.path.isdir(step_name):
            return 'not_in'
        step = self.project.read_step(step_name, no_check=True)
        if step.is_completed():
            return 'completed'
        if step.can_be_completed():
            return 'can_be_completed'
        return 'in_process'

    def steps_status(self):
        return dict((d, self.step_status(d)) for d in self.all_step_names())

    def _remove_steps_with_error(self):
        for sn in self.all_step_names():
            if os.path.isdir(sn) and not os.path.isfile(os.path.join(sn, 'description.yml')):
                print('Info: removing step with an error:', sn)
                remove_directory(sn)

    #
    def run_command(self, cmd, run_step):
        self._remove_steps_with_error()
        if cmd == 'run':
            self.cmd_run(run_step)
        else:
            getattr(self, self._COMMAND_METHODS[cmd])()

    def cmd_run(self, run_step):
        # Run all that can be run:
        #  * finish not finished step with output.zip in it.
        #  * run not run process action that has all dependencies satisfied

        # Finish not finished step with output.zip in it
        for s_obj in self.project.get_steps():
            if not s_obj.is_completed() and s_obj.can_be_completed():
                self.project.run_command_with_args('finish', s_obj.directory)

        # Run actions that have all dependencies satisfied
        to_finish = set()
        re_run = True
        something_done = False
        while re_run:
            # Note: after each calculation check actions, project structure can change!
            for action in self.actions():
                step_status = self.step_status(action.step_name)
                if step_status == 'not_in':
                    if all(self.step_status(s) == 'completed' for s in action.all_prev_steps):
                        cmd_args = list(action.cmd) + ['--step-name', action.step_name]
                        if run_step and action.has_run_switch():
                            cmd_args.append('-r')
                        self.project.run_command_with_args(*cmd_args)
                        something_done = True
                        break
                if step_status in ('in_process', 'can_be_completed'):
                    to_finish.add(action.step_name)
            else:
                re_run = False

        if to_finish:
            to_finish = '\n'.join(f' - {s}' for s in sorted(to_finish))
            print(f"""
There are steps to finish!
Check for INSTRUCTION.txt and calculate.zip in step(s):
{to_finish}
""")
        elif something_done:
            self.cmd_summary()

    def cmd_graph(self):
        node_styles = dict(not_in='dotted', completed='solid', in_process='dashed', can_be_completed='dashed,filled')
        edge_styles = dict(not_in='dashed', completed='solid', in_process='solid', can_be_completed='solid')
        # add_edge_styles = dict(not_in='dotted', completed='dotted', in_process='dotted', can_be_completed='dotted')
        status = self.steps_status()
        nodes = []
        edges = []
        for a in self.actions():
            sn = a.step_name
            nodes.append((sn, dict(label=sn, style=node_styles[status[sn]], shape='rectangle')))
            edges.extend((p, sn, dict(label=a.command, style=edge_styles[status[sn]])) for p in a.prev_steps)
            edges.extend((p, sn, dict(label=a.command, style='dotted')) for p in a.additional_reqs)

        create_graph_from_data(nodes, edges, 'workflow_graph')

    def cmd_summary(self):
        summary = self.get_summary()

        if text := summary.get('text'):
            print(text)
            write_str_in_file('workflow_summary.txt', text)

        if sheets := summary.get('sheets'):
            sheets_2_excel('workflow_summary.xlsx', sheets)

    def cmd_actions(self):
        for action in self.actions():
            print(f"{action.step_name:<22} {' '.join(action.cmd)}")
