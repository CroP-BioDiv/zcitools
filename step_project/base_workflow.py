import os
import itertools
from collections import namedtuple
from common_utils.exceptions import ZCItoolsValueError
from common_utils.cache import cache
from common_utils.file_utils import remove_directory, merge_zip_files, write_str_in_file
from .common.graph.project_graph import create_graph_from_data


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


class BaseWorkflow:
    _WORKFLOW = None  # Name of workflow
    _COMMAND_METHODS = dict(run='cmd_run', graph='cmd_graph', calcs='cmd_calcs', summary='cmd_summary')
    _SUMMARY_STEPS = None

    def __init__(self, project, parameters):
        self.project = project
        self.parameters = parameters

    @classmethod
    def all_commands(cls):
        return list(cls._COMMAND_METHODS.keys())

    @staticmethod
    def required_parameters():
        raise NotImplementedError('')

    @cache
    def actions(self):
        # Returns list of WfAction objects
        actions = []
        for x in self._actions():
            assert 2 <= len(x) <= 3, x
            sn, cmd = x[:2]
            a_prev = [] if len(x) == 2 else ([x[2]] if isinstance(x[2], str) else x[2])
            actions.append((sn, (cmd.split() if isinstance(cmd, str) else cmd), a_prev))

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

    def summary(self):
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
    def run_command(self, cmd):
        self._remove_steps_with_error()
        getattr(self, self._COMMAND_METHODS[cmd])()

    def cmd_run(self):
        # Run all that can be run:
        #  * finish not finished step with output.zip in it.
        #  * run not run process action that has all dependencies satisfied

        # Finish not finished step with output.zip in it
        for s_obj in self.project.get_steps():
            if not s_obj.is_completed() and s_obj.can_be_completed():
                self.project.run_command_with_args('finish', s_obj.directory)

        # Run not run process action that has all dependencies satisfied
        # ToDo: or not todo: topological sort, next to work, ...
        s_status = self.steps_status()
        re_run = True
        actions = self.actions()
        to_finish = set()
        while re_run:
            re_run = False
            for action in actions:
                sn = action.step_name
                if s_status[sn] == 'not_in':
                    if all(s_status[s] == 'completed' for s in action.all_prev_steps):
                        self.project.run_command_with_args(*action.cmd, '--step-name', sn)
                        re_run = True
                        # Update status
                        s_status[sn] = self.step_status(sn)
                if s_status[sn] in ('in_process', 'can_be_completed'):
                    to_finish.add(sn)
        if to_finish:
            print(f"""
There are steps to finish!
Check for INSTRUCTION.txt and calculate.zip in step(s):
{', '.join(sorted(to_finish))}
""")

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

    def cmd_calcs(self):
        # Collects data from calculate.zip's in unfinished steps
        to_zip = []
        for sn in self.all_step_names():
            if os.path.isdir(sn):
                step = self.project.read_step(sn, no_check=True)
                if not step.is_completed() and step.is_file('calculate.zip'):
                    to_zip.append(step.step_file('calculate.zip'))
        if not to_zip:
            print('No pending calculations!')
        else:
            merge_zip_files('calculate_pending.zip', to_zip, info=True)

    def cmd_summary(self):
        if not self._SUMMARY_STEPS:
            print('No summary steps!!!')
            return
        for s in self._SUMMARY_STEPS:
            pass

        # text = self.summary()
        # print(text)
        # write_str_in_file('workflow_summary.txt', text)
