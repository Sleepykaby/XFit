# -*- coding: utf-8 -*-
"================================================================================================================================="
import numpy as np
from scipy.optimize import curve_fit
from typing import *


def model(q, A, q0, xigma):
    I0 = A / (xigma * np.pi * np.sqrt(2))
    b = 1 / (2 * xigma**2)
    return I0 * np.exp(-b * (q - q0)**2)

def peak(q, params):
    num_func = int(len(params) / 3)
    I_list = list()
    for i in range(num_func):
        I = np.zeros_like(q)
        amp, ctr, wid= params[i*3: (i+1)*3]
        I += model(q, amp, ctr, wid)
        I_list.append(I)
    return np.array(I_list)

def total_fit(q, *params):
    I_list = peak(q, params)
    I_sum = np.sum(I_list, axis=0) + params[-1]
    return I_sum

def single_fit(q, *params):
    I_list = peak(q, params)
    return [I + params[-1] for I in I_list]

def peak_fit(function: Any, xdata: List[float], ydata: List[float], amorp_info: List[float], crystal_info: List[float]) -> Iterable[List]:
    peaks_total = list()
    lower_bounds = list()
    upper_bounds = list()
    background = min(ydata)+0.001
    
    peaks_total.append(amorp_info)
    for i in range(0, len(amorp_info), 3):
        # lower_bounds.append([amorp_info[i] - 1, amorp_info[i+1] - 0.2, amorp_info[i+2] - 0.1])
        # upper_bounds.append([amorp_info[i] + 1, amorp_info[i+1] + 0.2, amorp_info[i+2] + 0.1])
        lower_bounds.append([amorp_info[i] - 1000, amorp_info[i+1] - 0.2, amorp_info[i+2] - 10])
        upper_bounds.append([amorp_info[i] + 1000, amorp_info[i+1] + 0.2, amorp_info[i+2] + 10])
    
    for j in crystal_info:
        peaks_total.append([ydata[j], xdata[j], 0.1])
        lower_bounds.append([0, xdata[j] - 0.2, 0])
        upper_bounds.append([np.inf, xdata[j] + 0.2, 0.2])
    
    peaks_total = np.hstack(peaks_total + [background])
    lower_bounds = np.hstack(lower_bounds + [min(ydata)])
    upper_bounds = np.hstack(upper_bounds + [min(ydata)+0.002])

    popt, _ = curve_fit(function, xdata, ydata, 
                            bounds=(lower_bounds, upper_bounds), 
                            p0=peaks_total, maxfev=99999999)

    return popt

def correlation_coefficient(origin_data: List[float], fit_value: List[float]) -> Iterable[List]:
    residuals = origin_data - fit_value
    rss = np.sum(residuals**2)
    tss = np.sum((origin_data - np.mean(origin_data))**2)
    r_squared = 1 - (rss / tss)

    return residuals, r_squared

"================================================================================================================================="
if __name__ == "__main__":
    from DirFileRead import file_reader

    PATH = "F://1//230929//1D_data//deducted_and_masked//PA6T_WAXS_T80//"
    my_dict = file_reader(dir_path=PATH, start_value=5, end_value=22)
    x_test = my_dict["q"][1]
    y_test = my_dict["I"][1]
    
    from DataPreProcess import BaselineCorrector, smoothed_and_normalized, baseline_calculation
    
    y_normalized = smoothed_and_normalized(y_test, 10, 3)
    baseline_corrector = BaselineCorrector(x_test, y_normalized)

    selected_points_index = baseline_corrector.selected_points
    x1, x2 = x_test[selected_points_index[0]], x_test[selected_points_index[1]]
    y1, y2 = y_normalized[selected_points_index[0]], y_normalized[selected_points_index[1]]

    baseline = baseline_calculation(x_test, x1, x2, y1, y2)
    baseline_corrected = [max(0, y_val - baseline_val) if selected_points_index[0] <= i <= selected_points_index[1] else y_val for i, (y_val, baseline_val) in enumerate(zip(y_normalized, baseline))]

    from CrystalInfo import PeakInfo

    crystal_info = PeakInfo(x_test, baseline_corrected)
    crystal_info.cry_peaks_find()
    
    from InputCheck import get_valid_amorphous_params
    
    amorphous_params, _ = get_valid_amorphous_params()
    print(amorphous_params)

    popt = peak_fit(total_fit, x_test, baseline_corrected, amorphous_params, crystal_info.crystal_indices)
    print(popt)
    
    y_fit = total_fit(x_test, *popt)
    fit_list = single_fit(x_test, *popt)
    residuals, r_squared = correlation_coefficient(baseline_corrected, y_fit)

    print(f"拟合参数: {popt}")
    print(f"相关系数: {r_squared}")

    import matplotlib.pyplot as plt
    import matplotlib.cm as cm

    plt.scatter(x_test, baseline_corrected, marker='o', c="None", edgecolors='gray', label='Origin Data')
    plt.plot(x_test, y_fit, ls='--', lw=1.2, c='r', label='Fit Data')
    
    for i, single_peak in enumerate(fit_list):
        plt.plot(x_test, single_peak, ls='--', c=cm.rainbow(i/len(fit_list)), lw=1.0)
    
    plt.xlabel("q (nm$^{-1}$)", fontdict={'size': 12})
    plt.ylabel("Intenisty (a.u.)", fontdict={'size': 12})
    plt.legend(loc='upper right')
    plt.show()
