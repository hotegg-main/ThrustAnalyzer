import matplotlib.pyplot as plt
import os
import numpy as np

# グラフの描画設定
plt.rcParams['font.family'] = 'Arial'
plt.rcParams['font.size']   = 12
plt.rcParams['figure.titlesize'] = 13
plt.rcParams["xtick.direction"]   = "in"
plt.rcParams["ytick.direction"]   = "in"
plt.rcParams["xtick.top"]         = True
plt.rcParams["ytick.right"]       = True
plt.rcParams["xtick.major.width"] = 1.5
plt.rcParams["ytick.major.width"] = 1.5
plt.rcParams["axes.linewidth"] = 1.5

color = ['#FF4B00',
         '#005AFF',
         '#03AF7A',
         '#4DC4FF',
         '#F6AA00',
         ]


def plot_thrust_curve(path, time, force):

    plt.figure(figsize=(8, 4))
    plt.plot(time, force, color=color[0])
    plt.xlim(left=0., right=time[-1])
    plt.xlabel('Time [sec]')
    plt.ylabel('Force [N]')
    plt.minorticks_on()
    plt.grid(linestyle='--')
    plt.savefig(path + os.sep + 'thrust.png')
    plt.close()

def compare_thrust_curve(path, time1, force1, time2, force2, time3, force3):

    plt.figure(figsize=(8, 4))
    plt.plot(time1, force1, color=color[0], label='Raw', linewidth=3)
    plt.plot(time2, force2, color=color[1], label='LPF', linewidth=3)
    # plt.plot(time3, force3, color=color[2], label='LPF+Ave', linestyle='--', linewidth=3)
    plt.plot(time3, force3, color=color[2], label='LPF+Ave', marker='o', markersize=4)
    plt.xlim(left=0., right=time2[-1])
    # plt.xlim(left=9., right=10.)
    plt.xlabel('Time [sec]')
    plt.ylabel('Force [N]')
    plt.minorticks_on()
    plt.legend()
    plt.grid(linestyle='--')
    plt.savefig(path + os.sep + 'thrust_compare.png')
    plt.close()

def plot_gradient(path, time, grad):

    grad_max = np.max(abs(grad))

    plt.figure(figsize=(8, 4))
    plt.plot(time, grad / grad_max, color=color[0])
    # plt.plot(time, grad, color=color[0])
    plt.xlim(left=0., right=time[-1])
    # plt.ylim(bottom=-1000., top=1000.)
    # plt.ylim(bottom=-5000., top=5000.)
    # plt.ylim(top=1000.)
    plt.xlabel('Time [sec]')
    plt.ylabel('Gradient [N/s]')
    plt.minorticks_on()
    plt.grid(linestyle='--')
    plt.savefig(path + os.sep + 'gradient.png')
    plt.close()

def plot_gradient2(path, time, grad, time2, grad2, time3, grad3):

    plt.figure(figsize=(8, 4))
    plt.plot(time, grad, color=color[0])
    plt.scatter(time2, grad2, color=color[1])
    plt.scatter(time3, grad3, color=color[2])
    plt.xlim(left=0., right=time[-1])
    # plt.ylim(bottom=-1000., top=1000.)
    # plt.ylim(top=1000.)
    # plt.ylim(bottom=-1000.)
    plt.xlabel('Time [sec]')
    plt.ylabel('Gradient [N/s]')
    plt.minorticks_on()
    plt.grid(linestyle='--')
    plt.savefig(path + os.sep + 'gradient2.png')
    plt.close()

def plot_std(path, time, std):

    plt.figure('Standard')
    plt.plot(time, std, color=color[0], marker='o')
    plt.xlim(left=0., right=time[-1])
    # plt.ylim(bottom=0., top=1500.)
    plt.xlabel('Time [sec]')
    plt.ylabel('Gradient [N/s]')
    plt.minorticks_on()
    plt.grid(linestyle='--')
    plt.savefig(path + os.sep + 'std.png')
    plt.close()

def plot_temp(x, y, xlabel='x', ylabel='y'):

    plt.plot(x, y, marker='o')
    plt.grid()
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.show()

def plot_temp2(x, y, point1, point2, xlabel='x', ylabel='y'):

    plt.plot(x, y, marker='o')
    plt.axvline(point1, linestyle='--', color='black')
    plt.axvline(point2, linestyle='--', color='gold')
    plt.grid()
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.show()

def plot_log(x, y):

    plt.plot(x, y, marker='o')
    plt.yscale('log')
    plt.grid()
    plt.show()

def plot_compare(x1, y1, x2, y2):

    plt.plot(x1, y1, linewidth=3)
    plt.plot(x2, y2, linestyle='--')
    # plt.axvline(point1, linestyle='--', color='black')
    # plt.axvline(point2, linestyle='--', color='gold')
    plt.grid()
    # plt.xlabel(xlabel)
    # plt.ylabel(ylabel)
    plt.show()
