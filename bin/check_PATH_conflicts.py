#!/usr/bin/python3
import os
from collections import defaultdict

exe_2_path = defaultdict(list)
processed_paths = set()
for path in os.environ.get('PATH', '').split(os.pathsep):
    if path not in processed_paths:
        processed_paths.add(path)
        if os.path.isdir(path):
            for f in os.listdir(path):
                fp = os.path.join(path, f)
                if os.path.isfile(fp) and os.access(fp, os.X_OK):
                    exe_2_path[f].append(path)


for exe, paths in sorted(exe_2_path.items()):
    if len(paths) > 1:
        # Remove /usr/bin, /bin, /usr/local/bin, /usr/sbin, /sbin combination
        if all(p in ('/usr/bin', '/bin', '/usr/local/bin', '/usr/sbin', '/sbin') for p in paths):
            continue
        print(f"{exe:<20} : {', '.join(paths)}")
