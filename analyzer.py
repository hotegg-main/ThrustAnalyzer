import numpy as np
from scipy import integrate, signal, interpolate
from scipy.ndimage import gaussian_filter

from grapher import plot_temp as plot
from grapher import plot_log

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

def serch_max_point(thrust):
    '''
    最大推力となるインデックスを探す
    '''

    return np.argmax(thrust)

def __butter_lowpass(lowcut, fs, order=4):
    '''
    バターワースローパスフィルタを設計する関数
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

def thin_out_data(index, x, y, index_sta, index_end):
    '''
    平均処理を行って、データを間引く
    Args:
        index       :間引くデータの区切りのインデックスの配列
        x           :間引く元データのx
        y           :間引く元データのy
        index_sta   :
        index_end   :
    '''

    x_out = []
    y_out = []

    ista = 0

    time_step_arr = []
    index_break = len(index) - 1
    
    # 平均処理
    for i in range(len(index) - 1):

        ista = index[i]
        iend = index[i+1]

        time_step_arr.append(x[iend] - x[ista])

        if ista > index_sta and iend <= index_end:
            # 間引く対象の範囲内なら、平均処理実施

            x_thin = np.mean(x[ista:iend])
            y_thin = np.mean(y[ista:iend])
            x_out.append(x_thin)
            y_out.append(y_thin)

    # 平均処理前後で誤差が大きい場合は切り詰める
    error = compare_data(x[index_sta:index_end], y[index_sta:index_end], x_out, y_out)
    # plot(x[index_sta:index_end], error)
    
    index_max_in  = np.min([index_sta + np.argmax(error >= 10.), index_end])
    if index_max_in < int(index_end / 2):
        index_max_in = index_end
    index_max_out = np.argmax(x_out >= np.min([x_out[-1], x[index_max_in]]))
    
    # print('Time Step Ave.', np.mean(time_step_arr) + np.std(time_step_arr))
    # plot(index[:len(time_step_arr)], time_step_arr, 'Time [sec]', 'Elapse [sec]')
    
    x_out = np.array(x_out[:index_max_out])
    y_out = gaussian_filter(np.array(y_out[:index_max_out]), 1.5)
    # y_out = np.array(y_out[:index_max_out])
    x_out = np.append(np.append(x[:index_sta], x_out), x[index_max_in:])
    y_out = np.append(np.append(y[:index_sta], y_out), y[index_max_in:])

    return x_out, y_out

def compare_data(time_src, thrust_src, time_dst, thrust_dst):
    '''
    Args:
    '''

    func_thrust1 = interpolate.interp1d(time_dst, thrust_dst, kind='linear', bounds_error=False, fill_value=(thrust_dst[0], thrust_dst[-1]))

    error = np.array([np.abs(thrust - func_thrust1(time)) for time, thrust in zip(time_src, thrust_src)])
    error /= thrust_src
    error *= 100.

    return error

def calc_fft(time, thrust):

    T = np.mean(time[1:] - time[:-1])
    N = len(thrust)
    window = np.hanning(N)

    # FFT実行
    fft_result = np.fft.fft(thrust * window)
    fft_result = fft_result / (N / 2)
    amp = np.abs(fft_result)

    # 周波数軸の計算
    freqs = np.fft.fftfreq(len(thrust), T)

    # plot(freqs[:N//2], amp[:N//2])
    plot_log(freqs[:N//2], amp[:N//2])

if __name__=='__main__':

    N = 1024            # サンプル数
    dt = 0.0005          # サンプリング周期 [s]
    f1, f2 = 100, 350    # 周波数 [Hz]

    t = np.arange(0, N * dt, dt) # 時間 [s]
    x = 1.5 * np.sin(2 * np.pi * f1 * t) + np.sin(2 * np.pi * f2 * t) # データ

    calc_fft(t, x)