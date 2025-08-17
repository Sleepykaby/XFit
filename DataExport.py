# -*- coding: utf-8 -*-
"================================================================================================================================="
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from scipy import integrate
import csv
from typing import *


# def post_processing(function: Any, lower_bound: float, upper_bound: float, params: Any, file_name: str) -> Iterable[List]:
#     peak_area = list()
#     fwhm = list()
#     num = int(len(params)/3)
#     for i in range(num):
#         amp, ctr, wid= params[i*3: (i+1)*3]
#         area, _ = integrate.quad(function, lower_bound, upper_bound, args=(amp, ctr, wid))
#         peak_area.append(area)
#         fwhm.append(wid)

#     with open("crystallinity.csv", mode='a', newline='') as file:
#         csv_writer = csv.writer(file)
#         csv_writer.writerow([file_name] + peak_area)

#     with open("grain_size.csv", mode='a', newline='') as file:
#         csv_writer = csv.writer(file)
#         csv_writer.writerow([file_name] + fwhm)

def post_processing(function: Any, lower_bound: float, upper_bound: float, 
                   params: Any, file_name: str, x_points: int = 1000) -> Iterable[List]:
    """
    后处理函数：计算峰面积、半高宽，并保存拟合数据
    
    Parameters:
    - function: 拟合函数
    - lower_bound: 积分下界
    - upper_bound: 积分上界  
    - params: 拟合参数
    - file_name: 文件名
    - x_points: x轴数据点数，默认1000
    """
    peak_area = []
    fwhm = []
    num_peaks = len(params) // 3
    
    # 生成x轴数据
    x_data = np.linspace(lower_bound, upper_bound, x_points)
    
    # 计算每个峰的参数和总拟合曲线
    individual_peaks = []
    total_fit = np.zeros_like(x_data)
    
    for i in range(num_peaks):
        amp, ctr, wid = params[i*3:(i+1)*3]
        
        # 计算峰面积
        area, _ = integrate.quad(function, lower_bound, upper_bound, args=(amp, ctr, wid))
        peak_area.append(area)
        
        # FWHM
        fwhm.append(wid)
        
        # 计算单个峰的y值
        peak_y = np.array([function(x, amp, ctr, wid) for x in x_data])
        individual_peaks.append(peak_y)
        total_fit += peak_y
    
    # 保存结晶度数据 (峰面积)
    with open("crystallinity.csv", mode='a', newline='', encoding='utf-8') as file:
        csv_writer = csv.writer(file)
        csv_writer.writerow([file_name] + peak_area)
    
    # 保存晶粒尺寸数据 (半高宽)
    with open("grain_size.csv", mode='a', newline='', encoding='utf-8') as file:
        csv_writer = csv.writer(file)
        csv_writer.writerow([file_name] + fwhm)
    
    # 保存拟合曲线数据 (x值和所有峰的y值)
    fit_file_name = f"fit_curves_{file_name.replace('.', '_')}.csv"
    with open(fit_file_name, mode='w', newline='', encoding='utf-8') as file:
        csv_writer = csv.writer(file)
        
        # 写入表头
        headers = ['x_values', 'total_fit'] + [f'peak_{i+1}' for i in range(num_peaks)]
        csv_writer.writerow(headers)
        
        # 写入数据
        for j in range(len(x_data)):
            row = [x_data[j], total_fit[j]] + [peak[j] for peak in individual_peaks]
            csv_writer.writerow(row)


def fitplot_window(xdata, ydata, cry_indices, fit, fit_list, filename, residuals, corr_coefficient, params):
    baseline = np.full_like(xdata, fill_value=params[-1])
    
    fig, (ax1, ax2) = plt.subplots(nrows=2, sharex=False, gridspec_kw={'height_ratios': [3, 1.8]})
    fig.subplots_adjust(hspace=0.4)

    ax1.plot(xdata, ydata, ls='-', c='gray', lw=1, label='origin')
    ax1.plot(xdata, fit, ls='--', c='red', lw=0.7, label='fit')
    ax1.set_xlabel("q (nm$^{-1}$)", fontdict={'size': 12})
    ax1.set_ylabel("Intensity (a.u.)", fontdict={'size': 12})
    ax1.legend(loc='upper left')
    ax1.grid(True, which='both')
    
    if fit_list is not None:
        for i, single_fit in enumerate(fit_list):
            ax1.plot(xdata, single_fit, ls='--', c=cm.rainbow(i/len(fit_list)), lw=0.8)
            ax1.fill_between(xdata, single_fit, baseline, facecolor=cm.rainbow(i/len(fit_list)), alpha=0.5)
        
    for index in cry_indices:
        ax1.plot(xdata[index], ydata[index], 'bo')

    ax2.scatter(xdata, residuals, s=7, color='None', marker='o', edgecolors='black')
    ax2.axhline(y=0, c='red', ls='--', lw=0.8)
    ax2.set_xlabel("q (nm$^{-1}$)", fontdict={'size': 12})
    ax2.set_ylabel("Residuals", fontdict={'size': 12})
    ax2.text(0.03, 0.9, f'R\u00B2 = {corr_coefficient}', transform=ax2.transAxes, fontsize=10,
            verticalalignment='top', bbox=dict(facecolor='white', edgecolor='gray', alpha=0.5))
    
    plt.savefig(f"{filename}.png")
    plt.close()


"================================================================================================================================="
if __name__ == "__main__":
    from DirFileRead import file_reader
    from pathlib import Path

    PATH = "F://1//230929//1D_data//deducted_and_masked//PA6T_WAXS_T80//"
    my_dict = file_reader(dir_path=PATH, start_value=5, end_value=22)
    x_test = my_dict["q"][1]
    y_test = my_dict["I"][1]

    dir_path = Path(PATH)
    files_list = [file.name for file in dir_path.iterdir() if file.is_file()]
    file_name = files_list[0]
    
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

    header = [f"q = {loc}" for loc in amorphous_params[1::3]]
    for index in crystal_info.crystal_indices:
        header.append(f"q = {x_test[index]}")

    with open("crystallinity.csv", mode='w', newline='') as file:
        csv_writer = csv.writer(file)        
        csv_writer.writerow(header)

    with open("grain_size.csv", mode='w', newline='') as file:
        csv_writer = csv.writer(file)
        csv_writer.writerow(header)
    
    from CurveFit import *

    popt = peak_fit(total_fit, x_test, baseline_corrected, amorphous_params, crystal_info.crystal_indices)
    
    y_fit = total_fit(x_test, *popt)
    fit_list = single_fit(x_test, *popt)
    residuals, r_squared = correlation_coefficient(baseline_corrected, y_fit)

    post_processing(model, x_test[0], x_test[-1], popt)
    fitplot_window(x_test, baseline_corrected, crystal_info.crystal_indices, y_fit, fit_list, file_name, residuals, r_squared, popt)