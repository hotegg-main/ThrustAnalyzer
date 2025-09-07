import numpy as np

def load_thrust(path):

    raw_data = np.loadtxt(path, delimiter=',', skiprows=1, usecols=(0, 1), encoding='utf-8')

    time_raw = raw_data[:, 0]
    thrust_raw = raw_data[:, 1]
    dt = time_raw[1] - time_raw[0]

    return time_raw, thrust_raw, dt