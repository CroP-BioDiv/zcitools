import os
import shutil
import datetime
import yaml
import multiprocessing
import subprocess
from zipfile import ZipFile, ZIP_DEFLATED

_TIMERUN_FILENAME = 'run_info.txt'


def lasted_str(lasted):
    days = lasted // (24 * 60 * 60)
    desc = f'{days}d:' if days else ''
    rest = lasted - days * 24 * 3600
    hours = rest // 3600
    desc += f'{hours:02}h:' if desc or hours else ''
    rest -= hours * 3600
    mins = rest // 60
    desc += f'{mins:02}m:' if desc or mins else ''
    secs = rest - mins * 60
    desc += f'{secs:02}s'
    return desc


class LogRun:
    def __init__(self, **data):
        started = datetime.datetime.now()
        self.started = started.timestamp()
        with open(_TIMERUN_FILENAME, 'w') as _out:
            _out.write(f'started: {started}\n')
            if data:
                _out.write('\n'.join(f'{k}: {v}' for k, v in data.items()))
                _out.write('\n')

    def finish(self, **data):
        ended = datetime.datetime.now()
        lasted = int(ended.timestamp() - self.started)  # Just seconds
        desc = lasted_str(lasted)
        print(f'Calculation lasted: {desc} ({lasted}s)')

        with open(_TIMERUN_FILENAME, 'a') as _out:
            if data:
                _out.write('\n'.join(f'{k}: {v}' for k, v in data.items()))
                _out.write('\n')
            _out.write(f'ended: {ended}\n')
            _out.write(f'lasted: {desc}\n')


def find_exe(default_exe, env_var, install_instructions, raise_desc):
    exe = os.getenv(env_var, default_exe)
    if not shutil.which(exe):
        print(install_instructions.format(exe=default_exe, env_var=env_var))
        if raise_desc:
            raise ValueError(f'No {raise_desc} installed! Tried {exe}')
        return
    return exe


def get_num_threads():
    return multiprocessing.cpu_count()


def get_num_logical_threads():
    try:
        import psutil
        return psutil.cpu_count(logical=True)
    except ImportError:
        return max(1, (get_num_threads() // 2))


def load_finish_yml():
    with open('finish.yml', 'r') as r:
        return yaml.load(r, Loader=yaml.CLoader)


def run_cmd(cmd, cwd=None, output_file=None):
    cmd = [str(x) for x in cmd]  # Fix arguments

    # Print
    print(f"Cmd: {f'cd {cwd}; ' if cwd else ''}{' '.join(cmd)}{f' > {output_file}' if output_file else ''}")

    # Run
    if output_file:
        with open(output_file, 'w') as _out:
            subprocess.run(cmd, cwd=cwd, stdout=_out)
    else:
        subprocess.run(cmd, cwd=cwd)


def zip_output(files, cwd=None, zip_filename='output.zip', skip_missing=False):
    if cwd:
        os.chdir(cwd)

    with ZipFile(zip_filename, 'w', compression=ZIP_DEFLATED) as output:
        for f_name in files:
            if os.path.isfile(f_name):
                output.write(f_name)
            elif not skip_missing:
                raise ValueError(f'File {f_name} is missing!')
        if os.path.isfile(_TIMERUN_FILENAME):
            output.write(_TIMERUN_FILENAME)
