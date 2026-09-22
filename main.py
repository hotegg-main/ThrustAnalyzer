import scipy.ndimage
from reader import *
from grapher import *
from analyzer import *
import scipy
import os

def prepare(path_file, path_result):
    '''
    前処理
    Args:
        path_file   :推力履歴のファイルパス
        path_result :結果出力先のフォルダパス
    Returns:
        time_raw    :時間
        thrust_raw  :推力
        dt          :時間刻み
        info_thrust :推力情報
    '''

    if not os.path.isdir(path_result):
        os.mkdir(path_result)

    time_raw, thrust_raw, dt = load_thrust(path_file)
    info_thrust = calc_thrust_info(time_raw, thrust_raw)

    index_sta = np.argmax(time_raw >= info_thrust['Act. Start Time'])
    index_end = np.argmax(time_raw >= info_thrust['Act Time, Time'] + info_thrust['Act. Start Time'])
    time_raw -= time_raw[index_sta]
    info_thrust['Act. Start Time'] = 0.
    
    return time_raw[index_sta:index_end], thrust_raw[index_sta:index_end], dt, info_thrust

def apply_LPF(time, thrust, dt):
    '''
    LPFの適用
    Args:
        time    :時間
        thrust  :推力
        dt      :時間刻み
    '''

    thrust_lpf = butter_lowpass_filter(thrust, 10, 1./ dt)
    time_lpf   = time

    return time_lpf, thrust_lpf

def apply_thin_out(time, thrust):

    index_max    = serch_max_point(thrust)

    grad = calc_gradient(time, thrust)

    # Positive Slope
    index_p, time_p, grad_p = detect_positive_slope(time, grad)
    # Negative Slope
    index_n, time_n, grad_n = detect_negative_slope(time, grad)
    # Merge Peak Index Info.
    index_peak = merge_index(index_p, index_n)
    
    if len(index_peak) > 10:
        index_thin_sta = index_p[np.argmax(index_p > index_max)]
        time_thin, thrust_thin = thin_out_data(index_peak, time, thrust, index_thin_sta, index_peak[-1])
        return time_thin, thrust_thin, index_peak
    
    else:
        return time, thrust, index_peak

def apply_gaussian(time, thrust):

    y_out_interp = interpolate.interp1d(time, thrust, kind='linear', bounds_error=False, fill_value=(thrust[0], thrust[-1]))
    time_gauss = np.arange(time[0], time[-1] + 0.05, 0.05)
    thrust_gauss = gaussian_filter(y_out_interp(time_gauss), 1.)

    return time_gauss, thrust_gauss

def calc_burnout(time, thrust):
    '''
    燃焼時間の計算
    '''

    grad = calc_gradient(time, thrust)                          # 勾配の算出
    index_1, index_2 = serch_decrease_peak(time, thrust, grad)  # 燃焼時間前後の変曲点の算出
    tan1, tan2, nor  = calc_tangent_aft_tangent_bisector(time, thrust, grad, index_1, index_2)
                                                                # 後方接線角二等分線の算出
    time_sta = time[index_1] - 2. * (time[index_2] - time[index_1])
    time_end = time[index_2] + 2. * (time[index_2] - time[index_1])
    x_cross = calc_cross_point_of_backward_bisect(time, thrust, nor, time_sta, time_end)
                                                                # 後方接戦角二等分線と推力カーブとの交点の算出
    time_array = np.array([time[index_1], time[index_2]])
    tan2_array = np.array([time_array, tan2(time_array)])       # 接戦2の座標
    time_array = np.array([time_sta, time_end])
    tan1_array = np.array([time_array, tan1(time_array)])       # 接戦1の座標
    nor_array  = np.array([time_array, nor(time_array)])        # 後方接戦角二等分線と推力カーブの座標

    return x_cross, tan1_array, tan2_array, nor_array

def make_summay(path, info_thrust):

    text = \
    [
        'Burn Start Time [sec]      :', str(round(info_thrust['Act. Start Time'], 3)),'\n',
        'Actuation Time [sec]       :', str(round(info_thrust['Act Time, Time'], 3)),'\n',
        'Burnout Time [sec]         :', str(round(info_thrust['Burnout Time, Time'], 3)),'\n',
        'Total Impulse(Act) [Ns]    :', str(round(info_thrust['Total Impulse(Act)'], 3)),'\n',
        'Ave. Thrust(Act) [N]       :', str(round(info_thrust['Ave. Thrust(Act)'], 3)),'\n',
        'Total Impulse(Burn) [Ns]   :', str(round(info_thrust['Total Impulse(Burn)'], 3)),'\n',
        'Ave. Thrust(Burn) [N]      :', str(round(info_thrust['Ave. Thrust(Burn)'], 3)),'\n',
    ]

    with open(path + os.sep +'summary.txt', mode='w') as f:

        f.writelines(text)

def output_thrust(path, time_raw, thrust_raw, time_filter, thrust_filter):

    data_raw    = np.c_[time_raw, thrust_raw]
    data_filter = np.c_[time_filter, thrust_filter]
    
    np.savetxt(path + os.sep +'thrust_raw.csv'      , data_raw      , delimiter=',', header='time[s],thrust[N]', comments='')
    np.savetxt(path + os.sep +'thrust_filter.csv'   , data_filter   , delimiter=',', header='time[s],thrust[N]', comments='')

def main(path_file, path_result):

    #############################################################################
    # Pre. Proc.                                                                #
    #############################################################################
    time_raw, thrust_raw, dt_raw, info_thrust = prepare(path_file, path_result)
    
    #############################################################################
    # Filter                                                                    #
    #############################################################################
    # 1st STEP: Low Pass Filter
    time_lpf, thrust_lpf = apply_LPF(time_raw, thrust_raw, dt_raw)

    # 2nd STEP: Data Thin Out (Average + Gaussian Filter)
    time_thin, thrust_thin, index_peak = apply_thin_out(time_lpf, thrust_lpf)
    time_gauss, thrust_gauss = apply_gaussian(time_thin, thrust_thin)

    time_filter   = time_gauss
    thrust_filter = thrust_gauss

    #############################################################################
    # Post Proc.                                                                #
    #############################################################################
    # Calc. Time Burnout
    time_burnout, tan1, tan2, nor = calc_burnout(time_filter, thrust_filter)
    info_thrust = update_info_buntout(info_thrust, time_burnout, time_raw, thrust_raw)

    make_summay(path_result, info_thrust)
    output_thrust(path_result, time_raw, thrust_raw, time_filter, thrust_filter)
    
    # Plot
    compare_thrust_curve(path_result, time_raw, thrust_raw, time_lpf, thrust_lpf, time_filter, thrust_filter, index_peak)
    plot_thrust_curve(path_result, time_raw, thrust_raw, time_filter, thrust_filter, info_thrust, tan1, tan2, nor)
    # plot_gradient(path_result, time_thin, grad_thin)

if __name__=='__main__':

    main('example/thrust.csv', 'output/sample')