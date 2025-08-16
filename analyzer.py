import numpy as np
from scipy import integrate, signal
from scipy.ndimage import gaussian_filter

def calc_total_impulse(time, thrust):

    return integrate.simpson(thrust, time)

def calc_gradient(time, thrust):
    
    return np.gradient(thrust, time)

def serch_static_point(time, thrust, grad):
    '''
    推力履歴の振動が収まる時刻を探す
    '''

    from grapher import plot_std

    time_fin = time[-1]
    # 何分割でサーチする？
    num_block = 50
    time_step = time_fin / 50

    time_arr = np.zeros(num_block)
    std_arr = np.zeros(num_block)

    time_sta = 0.
    time_end = time_sta + time_step
    point = 0
    time_static = 0.

    for i in range(num_block):

        idx_sta = np.argmax(time >= time_sta)
        idx_end = np.argmax(time >= time_end)

        time_arr[i] = time_end
        std_arr[i] = np.std(grad[idx_sta:idx_end])

        if std_arr[i] >= 200.:
            point = i
            time_static = time_end

        time_sta += time_step
        time_end += time_step

    plot_std('output', time_arr, std_arr)

    point = np.argmax(time >= time_static)

    return point


def __butter_lowpass(lowcut, fs, order=4):
    '''バターワースローパスフィルタを設計する関数
    '''
    nyq = 0.5 * fs
    low = lowcut / nyq
    b, a = signal.butter(order, low, btype='low')
    return b, a


def butter_lowpass_filter(x, lowcut, fs, order=4):
    '''
    データにローパスフィルタをかける関数
    Args:
        x: 元データ
        lowcut:カットする周波数
        fs: サンプリングレート
    '''
    b, a = __butter_lowpass(lowcut, fs, order=order)
    y = signal.filtfilt(b, a, x)
    return y

def detect_positive_slope(x, y):

    x_out = []
    y_out = []
    i_out = []

    x_temp = 0.
    y_temp = 0.

    for i in range(len(y)):

        if y[i] > 0.:

            if y[i] > y_temp:
                x_temp = x[i]
                y_temp = y[i]

        else:
            
            if y_temp != 0.:
                i_out.append(i)
                x_out.append(x_temp)
                y_out.append(y_temp)
            
            x_temp = 0.
            y_temp = 0.
        

    return np.array(i_out), np.array(x_out), np.array(y_out)

def detect_negative_slope(x, y):

    x_out = []
    y_out = []
    i_out = []

    x_temp = 0.
    y_temp = 0.

    for i in range(len(y) - 1):

        if y[i] < 0.:

            if y[i] < y_temp:
                x_temp = x[i]
                y_temp = y[i]

        else:
            
            if y_temp != 0.:
                i_out.append(i)
                x_out.append(x_temp)
                y_out.append(y_temp)
            x_temp = 0.
            y_temp = 0.
        

    return np.array(i_out), np.array(x_out), np.array(y_out)

def merge_index(index1, index2):
    '''
    配列を結合して、昇順ソート
    '''

    return np.sort(np.append(index1, index2))

def thin_out_data(index_p, x, y, index_max):
    '''
    
    '''

    x_out = []
    y_out = []

    ista = 0
    
    for idx_p in index_p:

        if idx_p <= index_max:

            iend = idx_p

            x_out.append(np.mean(x[ista:iend]))
            y_out.append(np.mean(y[ista:iend]))

            ista = iend

    x_out = np.array(x_out)
    y_out = gaussian_filter(np.array(y_out), 3)
    x_out = np.append(x_out, x[index_max:])
    y_out = np.append(y_out, y[index_max:])
    # y_out = np.append(y_out, gaussian_filter(y[index_max:], 3))

    return x_out, y_out