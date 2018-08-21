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

exec_path = os.path.abspath(sys.argv[1])
print(exec_path)
methods = ["baseline", "exp_single", "exp_double", "exp_redecimal"]
selected_methods = [methods[1]]
figs = []
for method_index, method in enumerate(selected_methods):
    command = [exec_path, "--short-output", "--method", method]
    subdiv_count = 2

    run_data = RunData()
    for i in range(127 * subdiv_count):
        offsets = pow(2, i / subdiv_count)
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
