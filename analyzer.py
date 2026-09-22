import numpy as np
from scipy import integrate, signal, interpolate
from scipy.ndimage import gaussian_filter
from scipy import optimize

from grapher import plot_temp2 as plot
from grapher import plot_temp, plot_compare, plot_gradient_deb, plot_compare_2ax, plot_compare_2ax_3series, plot_compare_2ax_4series
from grapher import plot_log

def calc_thrust_info(time, thrust):

    index_max       = np.argmax(thrust)                 # 最大推力時インデックス
    thrust_max      = thrust[index_max]                 # 最大推力
    thrust_min      = thrust_max * .05                  # 立ち上がり推力（審査書基準）
    index_sta       = np.argmax(thrust >= thrust_min)   # 立ち上がりインデックス
    time_sta        = time[index_sta]                   # 立ち上がり時刻
    if thrust[-1] <= thrust_min:
        index_end   = index_max + np.argmax(thrust[index_max:] <= thrust_min)   # 作動終了時インデックス
    else:
        index_end   = len(thrust) - 1
    time_end        = time[index_end]                   # 作動終了時刻
    time_act        = time_end - time_sta               # 作動時間
    total_impulse   = integrate.simpson(thrust[index_sta:index_end], time[index_sta:index_end]) # トータルインパルス
    thrust_ave      = total_impulse / time_act          # 平均推力（作動時間基準）

    info_thrust = \
    {
        'Max Thrust'        :thrust_max,
        'Act. Start Time'   :time_sta,
        'Act. End Time'     :time_end,
        'Act Time, Time'    :time_act,
        'Act Time, Force'   :thrust_min,
        'Total Impulse(Act)':total_impulse,
        'Ave. Thrust(Act)'  :thrust_ave,
    }

    return info_thrust

def update_info_buntout(info_thrust, time_burn, time, thrust):
    '''
    推力情報に燃焼時間の情報を更新する
    '''

    info_new        = info_thrust
    index_burn      = np.argmax(time >= time_burn)
    total_impulse   = integrate.simpson(thrust[:index_burn], time[:index_burn])
    thrust_ave      = total_impulse / (time_burn)
    info_new['Burn. End, Time']     = time_burn
    info_new['Burnout Time, Time']  = time_burn
    info_new['Total Impulse(Burn)'] = total_impulse
    info_new['Ave. Thrust(Burn)']   = thrust_ave

    return info_new

def calc_gradient(time, thrust):
    
    return np.gradient(thrust, time)

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
    '''
    dydxが正かつdydxが変曲点(d2y/dx2=0)となるポイントの検出
    '''

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
    '''
    dydxが負かつdydxが変曲点(d2y/dx2=0)となるポイントの検出
    '''

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

def func_linear_function(x0, y0, dydx0):
    '''
    ラムダ式で1次関数を返す
    '''

    return lambda x: dydx0 * (x - x0) + y0

def calc_tangent_aft_tangent_bisector(x, y, grad, index1, index2):
    '''
    後方折線角二等分線との交点を求める
    '''

    # 後方接線
    x1      = x[index1]
    y1      = y[index1]
    grad1   = grad[index1]
    x2      = x[index2]
    y2      = y[index2]
    grad2   = grad[index2]
    func_tan1 = func_linear_function(x1, y1, grad1)
    func_tan2 = func_linear_function(x2, y2, grad2)

    # 後方接線角二等分線
    x_nor = \
    optimize.bisect(lambda x:  func_tan1(x) - func_tan2(x),
                    x1, x2)
    y_nor = func_tan1(x_nor)
    theta_bisec   = np.arctan(grad2) \
                + .5 * np.arctan2(grad1 - grad2, 1. + grad1*grad2)
    grad_nor = - 1. / np.tan(theta_bisec)
    func_nor = func_linear_function(x_nor, y_nor, grad_nor)

    return func_tan1, func_tan2, func_nor

def calc_cross_point_of_backward_bisect(x, y, func_nor, x_sta, x_end):
    
    # 元の関数の線形近似
    func_org = interpolate.interp1d(x, y, kind='linear', bounds_error=False, fill_value=(y[0], y[-1]))
    
    # 交点
    x_cross = \
    optimize.bisect(lambda x:  func_org(x) - func_nor(x),
                    x_sta, x_end)
    
    return x_cross

def detect_inflection_point(y):
    '''
    変曲点を探す
    '''

    index_p = []
    index_n = []

    for i in range(len(y) - 2):

        if (y[i+1] > y[i]) and (y[i+1] > y[i+2]):
            # 正の変曲点
            index_p.append(i+1)

        elif (y[i+1] < y[i]) and (y[i+1] < y[i+2]):
            # 負の変曲点
            index_n.append(i+1)

    return np.array(index_p), np.array(index_n)

def serch_decrease_peak(x, y, grad):
    '''
    燃焼時間前後の傾きの変曲点を求める
    '''

    index_sta = np.argmax(y) # 最大推力から解析開始
    index_p, index_n = detect_inflection_point(grad[index_sta:]) # 変曲点の算出
    index_p += index_sta
    index_n += index_sta
    index_n = np.delete(index_n, 0)
    index = merge_index(index_p, index_n)   # 正負の変曲点を統合
    index_b = np.argmin(grad[index])        # 勾配が最小となる変曲点を算出
    index_a = index_b - 1                   # 勾配が最小となる変曲点の一つ前の変曲点を指定

    return index[index_a], index[index_b]

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

    # print(x[index_sta], x[index_end])
    
    # 平均処理
    for i in range(len(index) - 1):

        ista = index[i]
        iend = index[i+1]

        time_step_arr.append(x[iend] - x[ista])

        if ista > index_sta and iend <= index_end:
            # 間引く対象の範囲内なら、平均処理実施

            # x_thin = np.mean(x[ista:iend])
            x_thin = x[ista] if y[ista] >= y[iend] else x[iend]
            y_thin = np.mean(y[ista:iend])
            x_out.append(x_thin)
            y_out.append(y_thin)

    # 平均処理前後で誤差が大きい場合は切り詰める
    error = compare_data(x[index_sta:index_end], y[index_sta:index_end], x_out, y_out)
    # plot_compare_2ax_4series(x, y, x_out, y_out, x[index], y[index], x[index_sta:index_end], error)
    
    # index_max_in: 平滑前後の差が10%以上になる時刻　→定常燃焼終了区間を検出
    index_max_in  = np.argmax(error[int(len(error)/4):] >= 10.) # 前1/4以降の区間で差が10%以上になる時刻
    index_max_in  += index_sta + int(len(error)/4)
    index_max_in  = np.min([index_max_in, index_end])
    # if index_max_in < int(index_end / 2):
        # 平滑前後で差が10%以上になる時刻が、全時刻の前半に存在
        # print('True', x[index_max_in], x[index_end])
        # index_max_in = index_end
        # index_max_in = int(index_end / 2)

    index_max_out = np.argmin(x_out < np.min([x_out[-1], x[index_max_in]]))

    x_out = np.array(x_out[:index_max_out])
    y_out = np.array(y_out[:index_max_out])
    x_out = np.append(np.append(x[:index_sta], x_out), x[index_max_in:])
    y_out = np.append(np.append(y[:index_sta], y_out), y[index_max_in:])
    
    return x_out, y_out

    # plot(x, y, x[index_end], x[index_max_in])
    # plot_compare(x, y, x_out, y_out)
    # plot_compare_2ax(x, y, x[index_sta:index_end], error)
    # plot_compare_2ax_3series(x, y, x_out, y_out, x[index_sta:index_end], error)
    # plot_compare(x, y, x_filter, y_filter)

    # return x_filter, y_filter

def compare_data(time_src, thrust_src, time_dst, thrust_dst):
    '''
    Args:
    '''

    func_thrust1 = interpolate.interp1d(time_dst, thrust_dst, kind='linear', bounds_error=False, fill_value=(thrust_dst[0], thrust_dst[-1]))

    error = np.array([np.abs(thrust - func_thrust1(time)) for time, thrust in zip(time_src, thrust_src)])
    error /= thrust_src
    error *= 100.

    return error