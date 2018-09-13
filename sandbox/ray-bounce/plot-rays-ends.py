import argparse
import os
import subprocess
import sys
import re as re
import pprint

import matplotlib.pyplot as plt
import mpl_toolkits.mplot3d as plt3d

args_parser = argparse.ArgumentParser(description='Plots origins and targets for a ray-bounce test centered on {0, 0, 0}')
args_parser.add_argument('exec_path', type=str)
args_parser.add_argument('--origin_extents', type=str)
args_parser.add_argument('--target_extents', type=str)
args_parser.add_argument('--origin_subdiv', type=str)
args_parser.add_argument('--target_subdiv', type=str)
args_parser.add_argument('--origin_distance', type=str)
args = args_parser.parse_args()

exec_path = os.path.abspath(args.exec_path)

command = [exec_path, "--short-output", "--dry-run", "--target-offset", "0x0x0"]
if not args.origin_extents is None:
    command.extend(["--origin-extents", args.origin_extents])
if not args.target_extents is None:
    command.extend(["--target-extents", args.target_extents])
if not args.origin_subdiv is None:
    command.extend(["--origin-subdiv", args.origin_subdiv])
if not args.target_subdiv is None:
    command.extend(["--target-subdiv", args.target_subdiv])
if not args.origin_distance is None:
    command.extend(["--origin-distance", args.origin_distance])

command_out = subprocess.run(command, stdout=subprocess.PIPE, text=True)
print(command_out.stdout)

label_line_regexp = '[A-Za-z]+\n'
output_lists = [
    (label.strip(), [
        line.split(' ')
        for line in data.split('\n') if line != ''
    ])
    for label, data in zip(re.findall(label_line_regexp, command_out.stdout),
                           [
                               raw_data
                               for raw_data in re.split(label_line_regexp, command_out.stdout)
                               if raw_data != ''
                           ]
                           )
]

fig = plt.figure()
axes = fig.add_subplot(111, projection='3d')
for output in output_lists:
    axes.scatter([float(point[0])
                  for point in output[1]
                 ],
                 [float(point[1])
                  for point in output[1]
                 ],
                 [float(point[2])
                  for point in output[1]
                 ])

plt.show()
