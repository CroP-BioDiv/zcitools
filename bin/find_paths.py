#!/usr/bin/python3
import os


def find_paths(directory, shell_format=False):
    # Two layer directory structure: <dir>/<group>/<program>/current/bin
    bin_dirs = []
    for group in os.listdir(directory):
        g_dir = os.path.join(directory, group)
        if os.path.isdir(g_dir):
            print(g_dir)
            for program in os.listdir(g_dir):
                bin_dir = os.path.join(g_dir, program, 'current', 'bin')
                print('  ', bin_dir)
                if os.path.isdir(bin_dir):
                    bin_dirs.append(bin_dir)

    if shell_format:
        print(':'.join(os.path.abspath(d) for d in bin_dirs))


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description="Finds bin paths")

    # Input
    parser.add_argument('directory', help='Directory to look in')
    parser.add_argument('-s', '--shell-format', action='store_true', help="Returns paths delimited by ':'")

    params = parser.parse_args()
    if params.directory:
        find_paths(params.directory, shell_format=params.shell_format)
    elif not params.shell_format:
        print('Directory not set!')
