import scipy.ndimage
from reader import *
from grapher import *
from analyzer import *
import scipy
import os

def main(file, path_result):

    if not os.path.isdir(path_result):
        os.mkdir(path_result)

    time_raw, thrust_raw, dt = load_thrust(file)
    
    grad_raw = calc_gradient(time_raw, thrust_raw)
    print('Total Impulse [Ns]', calc_total_impulse(time_raw, thrust_raw))
    
    # 1st STEP: Low Pass Filter
    thrust_lpf = butter_lowpass_filter(thrust_raw, 10, 1./ dt)
    time_lpf = time_raw
    grad_lpf = calc_gradient(time_lpf, thrust_lpf)
    
    index_static = serch_static_point(time_raw, thrust_lpf, grad_lpf)
    index_max    = serch_max_point(thrust_lpf)

    # 2nd STEP: Data Thin Out (Average)
    # Positive Slope
    index_p, time_lpf_p, grad_lpf_p = detect_positive_slope(time_lpf, grad_lpf)
    # Negative Slope
    index_n, time_lpf_n, grad_lpf_n = detect_negative_slope(time_lpf, grad_lpf)
    index_peak = merge_index(index_p, index_n)
    index_thin_sta = index_n[np.argmax(index_n >= index_max)]
    time_thin, thrust_thin = thin_out_data(index_peak, time_lpf, thrust_lpf, index_thin_sta, index_static)

    # Plot
    plot_thrust_curve(path_result, time_raw, thrust_raw)
    compare_thrust_curve(path_result, time_raw, thrust_raw, time_lpf, thrust_lpf, time_thin, thrust_thin)
    plot_thrust_curve(path_result, time_thin, thrust_thin)
    plot_gradient(path_result, time_raw, grad_raw)
    plot_gradient2(path_result, time_raw, grad_lpf, time_lpf_p, grad_lpf_p, time_lpf_n, grad_lpf_n)


if __name__=='__main__':

    main('example/thrust.csv', 'output/sample')