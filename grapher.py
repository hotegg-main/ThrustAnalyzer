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


def plot_thrust_curve(path, time_raw, thrust_raw, time, thrust, info_thrust, tan1, tan2, nor):

    plt.figure(figsize=(8, 4))
    plt.plot(time_raw, thrust_raw, color=color[0])
    plt.plot(time, thrust, color=color[1])
    plt.plot(tan1[0], tan1[1], linestyle=':', color='black')
    plt.plot(tan2[0], tan2[1], linestyle=':', color='black')
    plt.plot(nor[0] , nor[1], linestyle='-.', color='black')
    plt.axvline(info_thrust['Burn. End, Time']  , linestyle='--', color='aqua')
    plt.axvline(info_thrust['Act. End Time']    , linestyle='--', color='aqua')
    plt.axhline(info_thrust['Max Thrust']       , linestyle='--', color='lime')
    plt.axhline(info_thrust['Ave. Thrust(Act)'] , linestyle='--', color='lime')
    plt.axhline(info_thrust['Ave. Thrust(Burn)'], linestyle='--', color='lime')
    plt.text(info_thrust['Burn. End, Time'], .0     , 'Burn. Time'          , color='aqua', rotation=90., verticalalignment='bottom', horizontalalignment='right')
    plt.text(info_thrust['Act. End Time'], .0       , 'Act. Time'           , color='aqua', rotation=90., verticalalignment='bottom', horizontalalignment='right')
    plt.text(.2, info_thrust['Max Thrust']          , 'Max. Thrust'         , color='lime', rotation=0., verticalalignment='top', horizontalalignment='left')
    plt.text(.2, info_thrust['Ave. Thrust(Act)']    , 'Ave. Thrust(Act)'    , color='lime', rotation=0., verticalalignment='top', horizontalalignment='left')
    plt.text(.2, info_thrust['Ave. Thrust(Burn)']   , 'Ave. Thrust(Burn)'   , color='lime', rotation=0., verticalalignment='top', horizontalalignment='left')
    plt.xlim(left=0., right=time[-1])
    plt.xlabel('Time [sec]')
    plt.ylabel('Force [N]')
    plt.minorticks_on()
    plt.grid(linestyle='--')
    plt.savefig(path + os.sep + 'thrust.png')
    plt.close()

def compare_thrust_curve(path, time1, force1, time2, force2, time3, force3, index_peak):

    plt.figure(figsize=(8, 4))
    plt.plot(time1, force1, color=color[0], label='Raw', linewidth=3)
    plt.plot(time2, force2, color=color[1], label='LPF', linewidth=3)
    plt.scatter(time2[index_peak], force2[index_peak], color='black', facecolors='none', s=100, zorder=2)
    plt.plot(time3, force3, color=color[2], label='LPF+Ave')
    plt.xlim(left=0., right=time2[-1])
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
    plt.xlim(left=0., right=time[-1])
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
    plt.xlabel('Time [sec]')
    plt.ylabel('Gradient [N/s]')
    plt.minorticks_on()
    plt.grid(linestyle='--')
    plt.savefig(path + os.sep + 'gradient2.png')
    plt.close()

def plot_gradient_deb(time, grad, time2, grad2, time3, grad3):

    plt.figure(figsize=(8, 4))
    plt.plot(time, grad, color=color[0])
    plt.scatter(time2, grad2, color=color[1])
    plt.scatter(time3, grad3, color=color[2])
    plt.xlim(left=0., right=time[-1])
    plt.xlabel('Time [sec]')
    plt.ylabel('Gradient [N/s]')
    plt.minorticks_on()
    plt.grid(linestyle='--')

    plt.show()

def plot_std(path, time, std):

    plt.figure('Standard')
    plt.plot(time, std, color=color[0], marker='o')
    plt.xlim(left=0., right=time[-1])
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
    plt.grid()
    plt.show()

def plot_compare_2ax(x1, y1, x2, y2):

    fig = plt.figure()
    ax1 = fig.add_subplot(111)
    ax2 = ax1.twinx()
    ax1.plot(x1, y1, color='royalblue')
    ax2.plot(x2, y2, color='orangered')
    ax1.grid()
    plt.show()

def plot_compare_2ax_3series(x1, y1, x2, y2, x3, y3):

    fig = plt.figure()
    ax1 = fig.add_subplot(111)
    ax2 = ax1.twinx()
    ax1.plot(x1, y1, color='royalblue')
    ax1.plot(x2, y2, color='orangered')
    ax2.plot(x3, y3, color='limegreen')
    ax1.grid()
    plt.show()
