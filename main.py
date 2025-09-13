import scipy.ndimage
from reader import *
from grapher import *
from analyzer import *
import scipy
import os

def prepare(path_file, path_result):

    if not os.path.isdir(path_result):
        os.mkdir(path_result)

    time_raw, thrust_raw, dt = load_thrust(path_file)
    info_thrust = calc_thrust_info(time_raw, thrust_raw)
    
    return time_raw, thrust_raw, dt, info_thrust

def apply_LPF(time, thrust, dt):
    '''
    Low Pass Filter(LPF)の適用
    '''

    thrust_lpf = butter_lowpass_filter(thrust, 10, 1./ dt)
    time_lpf = time

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
    index_thin_sta = index_p[np.argmax(index_p > index_max)]
    time_thin, thrust_thin = thin_out_data(index_peak, time, thrust, index_thin_sta, index_peak[-1])

    return time_thin, thrust_thin

def search_burnout(time, thrust, grad):
    '''
    燃焼時間の計算
    '''

    index_1, index_2 = serch_decrease_peak(time, thrust, grad)
    tan1, tan2, nor = calc_tangent_aft_tangent_bisector(time, thrust, grad, index_1, index_2)
    time_sta = time[index_1] - 2. * (time[index_2] - time[index_1])
    time_end = time[index_2] + 2. * (time[index_2] - time[index_1])
    x_cross = calc_cross_point_of_backward_bisect(time, thrust, nor, time_sta, time_end)
    time_array = np.array([time[index_1], time[index_2]])
    tan2_array = np.array([time_array, tan2(time_array)])
    time_array = np.array([time_sta, time_end])
    tan1_array = np.array([time_array, tan1(time_array)])
    nor_array  = np.array([time_array, nor(time_array)])

    return x_cross, tan1_array, tan2_array, nor_array

def main(path_file, path_result):

    time_raw, thrust_raw, dt_raw, info_thrust = prepare(path_file, path_result)
    
    # 1st STEP: Low Pass Filter
    time_lpf, thrust_lpf = apply_LPF(time_raw, thrust_raw, dt_raw)

    # 2nd STEP: Data Thin Out (Average + Gaussian Filter)
    time_thin, thrust_thin = apply_thin_out(time_lpf, thrust_lpf)

    grad_lpf = calc_gradient(time_lpf, thrust_lpf)
    grad_thin = calc_gradient(time_thin, thrust_thin)

    time_burnout, tan1, tan2, nor = search_burnout(time_thin, thrust_thin, grad_thin)
    info_thrust = updata_info_buntout(info_thrust, time_burnout, time_raw, thrust_raw)
    
    # Plot
    compare_thrust_curve(path_result, time_raw, thrust_raw, time_lpf, thrust_lpf, time_thin, thrust_thin)
    plot_thrust_curve(path_result, time_raw, thrust_raw, time_thin, thrust_thin, info_thrust, tan1, tan2, nor)
    plot_gradient(path_result, time_thin, grad_thin)

if __name__=='__main__':

    main('example/thrust.csv', 'output/sample')