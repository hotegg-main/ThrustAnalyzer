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
    
    print('Total Impulse [Ns]', calc_total_impulse(time_raw, thrust_raw))

    return time_raw, thrust_raw, dt

def apply_LPF(time, thrust, dt):
    '''
    Low Pass Filter(LPF)の適用
    '''

    thrust_lpf = butter_lowpass_filter(thrust, 10, 1./ dt)
    time_lpf = time

    return time_lpf, thrust_lpf

def apply_thin_out(time, thrust):

    # index_static = serch_static_point(time_raw, thrust_lpf, grad_lpf)
    index_max    = serch_max_point(thrust)

    grad = calc_gradient(time, thrust)

    # Positive Slope
    index_p, time_p, grad_p = detect_positive_slope(time, grad)
    # Negative Slope
    index_n, time_n, grad_n = detect_negative_slope(time, grad)
    # Merge Peak Index Info.
    index_peak = merge_index(index_p, index_n)
    index_thin_sta = index_n[np.argmax(index_n >= index_max)]
    time_thin, thrust_thin = thin_out_data(index_peak, time, thrust, index_thin_sta, index_peak[-1])

    return time_thin, thrust_thin

def main(path_file, path_result):

    time_raw, thrust_raw, dt_raw = prepare(path_file, path_result)
    
    # 1st STEP: Low Pass Filter
    time_lpf, thrust_lpf = apply_LPF(time_raw, thrust_raw, dt_raw)

    # 2nd STEP: Data Thin Out (Average + Gaussian Filter)
    time_thin, thrust_thin = apply_thin_out(time_lpf, thrust_lpf)

    grad_lpf = calc_gradient(time_lpf, thrust_lpf)
    grad_thin = calc_gradient(time_thin, thrust_thin)
    
    # Plot
    plot_thrust_curve(path_result, time_raw, thrust_raw)
    compare_thrust_curve(path_result, time_raw, thrust_raw, time_lpf, thrust_lpf, time_thin, thrust_thin)
    plot_thrust_curve(path_result, time_thin, thrust_thin)
    plot_gradient(path_result, time_thin, grad_thin)
    # plot_gradient2(path_result, time_raw, grad_lpf, time_lpf_p, grad_lpf_p, time_lpf_n, grad_lpf_n)


if __name__=='__main__':

    main('example/thrust.csv', 'output/sample')