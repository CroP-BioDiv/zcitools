import os
import itertools
from collections import namedtuple
from common_utils.exceptions import ZCItoolsValueError
from .common.graph.project_graph import create_graph_from_data


class WfAction(namedtuple('WfAction', 'step_name, prev_steps, cmd')):
    def __new__(cls, step_name, prev_steps, cmd):
        # cmd is a string
        assert isinstance(step_name, str), step_name

        # prev_steps is None, str, or list of string
        if isinstance(prev_steps, str):
            prev_steps = [prev_steps]
        elif not prev_steps:
            prev_steps = []
        assert isinstance(prev_steps, (list, tuple)), prev_steps
        assert all(isinstance(s, str) for s in prev_steps), prev_steps

        # args is a list
        assert cmd and isinstance(cmd, (list, tuple)), cmd
        assert all(isinstance(s, str) for s in cmd), cmd

        return super(WfAction, cls).__new__(cls, step_name, prev_steps, cmd)

    command = property(lambda self: self.cmd[0])


class BaseWorkflow:
    _WORKFLOW = None  # Name of workflow
    _COMMAND_METHODS = dict(run='cmd_run', graph='cmd_graph')

    def __init__(self, project, parameters):
        self.project = project
        self.parameters = parameters

    @classmethod
    def all_commands(cls):
        return list(cls._COMMAND_METHODS.keys())

    @staticmethod
    def required_parameters():
        raise NotImplementedError('')

    def actions(self):
        # Iterator of WfAction objects
        raise NotImplementedError('')

    def _check_actions(self, actions):
        # Check step names
        step_names = set(a.step_name for a in actions)
        prev_names = set(itertools.chain.from_iterable(a.prev_steps for a in actions))
        if (not_in := (prev_names - step_names)):
            raise ZCItoolsValueError(f"Worflow prev steps, not in project steps! {', '.join(sorted(not_in))}")
        # Check commands
        commands_map = self.project.commands_map
        if (not_in := [a.command for a in actions if a.command not in commands_map]):
            raise ZCItoolsValueError(f"Worflow actions, not existing command(s)! {', '.join(sorted(not_in))}")
        return actions

    #
    def all_steps(self):
        return sorted(a.step_name for a in self.actions())

    def steps_status(self):
        status = dict()
        for d in self.all_steps():
            if os.path.isdir(d):
                status[d] = 'completed' if self.project.read_step(d, no_check=True).is_completed() else 'in_process'
            else:
                status[d] = 'not_in'
        return status

    #
    def run_command(self, cmd):
        getattr(self, self._COMMAND_METHODS[cmd])()

    def cmd_run(self):
        # Run all that can be run:
        #  * finish not finished step with output.zip in it.
        #  * run not run process action that has all dependencies satisfied

        # Finish not finished step with output.zip in it
        for s_obj in self.project.get_steps():
            if not s_obj.is_completed() and s_obj.is_file('ouput.zip'):
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
                    if all(s_status[s] == 'completed' for s in action.prev_steps):
                        self.project.run_command_with_args(*action.cmd, '--fixed-name', sn)
                        re_run = True
                        # Update status
                        s_status[sn] = 'completed' \
                            if self.project.read_step(sn, no_check=True).is_completed() else 'in_process'
                if s_status[sn] == 'in_process':
                    to_finish.add(sn)
        if to_finish:
            print(f"""
There are steps to finish!
Check for INSTRUCTION.txt and calculate.zip in steps:
{', '.join(sorted(to_finish))}
""")

    def cmd_graph(self):
        styles = dict(not_in='dotted', completed='solid', in_process='dashedg')
        status = self.steps_status()
        nodes = []
        edges = []
        for a in self.actions():
            sn = a.step_name
            nodes.append((sn, dict(label=sn, style=styles[status[sn]])))
            edges.extend((p, sn, a.command) for p in a.prev_steps)

        create_graph_from_data(nodes, edges, 'workflow_graph')
