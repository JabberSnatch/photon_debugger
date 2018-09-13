import argparse
import os
import subprocess
import sys

import matplotlib.pyplot as plt

class RunData:
    def __init__(self):
        self.run_count = 1
        self.target_offset = []
        self.primary_fail = []
        self.inward_fail = []
        self.outward_fail = []
        self.max_bound = []

methods = ["baseline", "exp_single", "exp_double", "exp_redecimal"]

args_parser = argparse.ArgumentParser(description='Benchmarking utility for the ray-bounce test case')
args_parser.add_argument('exec_path', type=str)
args_parser.add_argument('methods', type=str, nargs='+', choices=methods)
args_parser.add_argument('--range', '-r', type=int, nargs=2, default=[0, 25])
args_parser.add_argument('--subdiv', '-s', type=int, default=1)

print("begin parse")
args = args_parser.parse_args()
print("end parse")

exec_path = os.path.abspath(args.exec_path)
selected_methods = args.methods
test_range = args.range
subdiv_count = args.subdiv

figs = []
xticks = []

for method_index, method in enumerate(selected_methods):
    command = [exec_path, "--short-output", "--method", method, "--origin-extents", "90x45"]

    run_data = RunData()
    for i in range(test_range[0] * subdiv_count, test_range[1] * subdiv_count):
        tick = pow(2, i / subdiv_count)
        xticks.append(tick)

        offsets = tick
        if offsets.is_integer():
            offsets = int(offsets)
        offset_option = ["--target-offset", str(offsets)+"x"+str(offsets)+"x1.0"]
        print(command + offset_option)
        command_out = subprocess.run(command + offset_option, stdout=subprocess.PIPE, text=True)
        print(command_out.stdout)
        split_output = command_out.stdout.split(' ')

        run_data.run_count = int(split_output[0])
        run_data.target_offset.append(offsets)
        run_data.primary_fail.append(int(split_output[1]))
        run_data.inward_fail.append(int(split_output[2]))
        run_data.outward_fail.append(int(split_output[3]))
        run_data.max_bound.append(max([int(x) for x in split_output[4:8]]))

    print(run_data.primary_fail)
    print(run_data.inward_fail)

    figs.append(plt.figure(clear=True))
    axes = plt.axes()
    axes.plot(run_data.target_offset, run_data.primary_fail)
    axes.plot(run_data.target_offset, run_data.inward_fail)
    axes.plot(run_data.target_offset, run_data.outward_fail)
    axes.plot(run_data.target_offset, run_data.max_bound)
    axes.set_xscale('log', basex=2)

plt.show()
