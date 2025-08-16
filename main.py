import scipy.ndimage
from reader import *
from grapher import *
from analyzer import *
import scipy

def main(file):

    time_raw, thrust_raw, dt = load_thrust(file)
    
    grad_raw = calc_gradient(time_raw, thrust_raw)
    print(calc_total_impulse(time_raw, thrust_raw))
    
    # Low Pass Filter
    thrust_lpf = butter_lowpass_filter(thrust_raw, 10, 1./ dt)
    time_lpf = time_raw
    # index_static = serch_static_point(time_raw, thrust_raw, grad_raw)
    grad_lpf = calc_gradient(time_lpf, thrust_lpf)
    
    index_static = serch_static_point(time_raw, thrust_lpf, grad_lpf)
    print('Time Static:', time_lpf[index_static])

    # Positive Slope
    index_p, time_lpf_p, grad_lpf_p = detect_positive_slope(time_lpf, grad_lpf)
    # Negative Slope
    index_n, time_lpf_n, grad_lpf_n = detect_negative_slope(time_lpf, grad_lpf)
    index_peak = merge_index(index_p, index_n)
    time_thin, thrust_thin = thin_out_data(index_peak, time_lpf, thrust_lpf, index_static)
    grad_thin = calc_gradient(time_thin, thrust_thin)

    serch_static_point(time_raw, thrust_raw, grad_raw)
    plot_thrust_curve('example', time_raw, thrust_raw)
    compare_thrust_curve('example', time_raw, thrust_raw, time_lpf, thrust_lpf)
    # compare_thrust_curve('example', time_raw, thrust_raw, time_thin, thrust_thin)
    # compare_thrust_curve('example', time_raw, thrust_raw, time_lpf, thrust_lpf, time_thin, thrust_thin)
    # plot_thrust_curve('example', time_lpf, thrust_lpf)
    plot_thrust_curve('example', time_thin, thrust_thin)
    plot_gradient('example', time_raw, grad_raw)
    plot_gradient('example', time_lpf, grad_lpf)
    # plot_gradient('example', time_thin, grad_thin)
    plot_gradient2('example', time_raw, grad_lpf, time_lpf_p, grad_lpf_p, time_lpf_n, grad_lpf_n)


if __name__=='__main__':

    main()