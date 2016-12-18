import matplotlib.pyplot as plt
import numpy as np


def plot(data, data1):
    plt.figure()
    plt.ion()
    for i in range(data.shape[0]):
        if i % 10 == 0:
            plt.plot(data[i], color='b', marker="o")
            plt.plot(data1[i], color='r', marker="v")
            plt.legend(("count", "real"))
            plt.title('Temperature_changes')
            plt.grid(True)
            plt.pause(0.01)
    plt.show(block=True)
    return

file = open("D:\\Documents\\sem5\\setki\\heat_conduction_1d_oop\\heat_conduction_1d_oop\\outputtest0.txt")
time_steps = int(file.readline())
#N = int(file.readline())

t = []

for i in range(time_steps):
    t.append(0.0)

dat1 = []
for i in range(time_steps):
    #stroka = file.readline()
    #t[i] = float(stroka)
    line = file.readline().split()
    dat1.append(np.array(line))
data1 = np.array(dat1)

file = open("D:\\Documents\\sem5\\setki\\heat_conduction_1d_oop\\heat_conduction_1d_oop\\output03.txt")
time_steps = int(file.readline())
#N = int(file.readline())

t = []

for i in range(time_steps):
    t.append(0.0)

dat = []
for i in range(time_steps):
    #stroka = file.readline()
    #t[i] = float(stroka)
    line = file.readline().split()
    dat.append(np.array(line))
data = np.array(dat)

plot(data, data1)
