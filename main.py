# -*- coding: utf-8 -*-
"================================================================================================================================="
from InputCheck import get_valid_directory_input, get_valid_file_input, get_valid_range_input, get_valid_amorphous_params
from DirFileRead import files_reader, single_file_reader
from DataPreProcess import BaselineCorrector, smoothed_and_normalized, baseline_calculation
from CrystalInfo import PeakInfo
from CurveFit import peak_fit, correlation_coefficient
from DataExport import post_processing, fitplot_window
from pathlib import Path
from scipy import integrate
from scipy.signal import savgol_filter
from scipy.optimize import curve_fit
import numpy as np
import matplotlib.pyplot as plt
import csv, itertools


def model(q, A, q0, xigma):
    # I0 = A / (xigma * np.pi * np.sqrt(2))
    I0 = A
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

def main_1(file_name, q_data, I_data, amorp_info, crystal_info, selected_points_index):
    I_normalized = smoothed_and_normalized(I_data, 10, 3)
    q1, q2 = q_data[selected_points_index[0]], q_data[selected_points_index[1]]
    I1, I2 = I_normalized[selected_points_index[0]], I_normalized[selected_points_index[1]]
    baseline = baseline_calculation(q_data, q1, q2, I1, I2)
    I_corrected = [max(0, I_val - baseline_val) if selected_points_index[0] <= i <= selected_points_index[1] else I_val for i, (I_val, baseline_val) in enumerate(zip(I_normalized, baseline))]

    popt = peak_fit(total_fit, q_data, I_corrected, amorp_info, crystal_info)
    I_fit = total_fit(q_data, *popt)
    fit_list = single_fit(q_data, *popt)
    residuals, r_squared = correlation_coefficient(I_corrected, I_fit)

    post_processing(model, q_data[0], q_data[-1], popt, file_name)
    fitplot_window(q_data, I_corrected, crystal_info, I_fit, fit_list, file_name, residuals, r_squared, popt)


def Maier_Saupe(x, A, φ0, β, I0):
    '''
    Maier-Saupe 分布函数, A表示峰强因子, φ0表示峰位, β表示峰宽因子, I0表示自由基线
    '''
    return I0 + A*np.exp(β*(np.cos((x-φ0)*np.pi/180))**2)

def Pnc_calculatation(popt):
    '''
    Maier-Saupe 平均场理论计算取向因子
    '''
    def f1(x,β):
        return ((3*(x**2)-1)/2)*np.exp(β*(x**2))
    def f2(x,β):
        return np.exp(β*(x**2))
    _, _, wid, _ = popt
    value1, _ = integrate.quad(f1, -1, 1, args=(wid))
    value2, _ = integrate.quad(f2, -1, 1, args=(wid))
    Pnc = value1/value2
    return Pnc

"================================================================================================================================="
if __name__ == "__main__":
    string = '--*欢迎使用 WAXD Batch Processing version2.0 批处理软件,\n请按照提示进行输入,\n当前版本支持结晶度、晶粒尺寸、取向因子的批量分析功能*--'
    print(string)
    print('=' * len(string))
    
    key = input(f"请输入你想要执行的分析功能('结晶度和晶粒尺寸'，或'取向因子'): ")
    if key == "结晶度和晶粒尺寸":
        # Initialization
        PATH = get_valid_directory_input("请输入要处理文件夹的路径(例如: F://1//test//): ")
        calibrate_file = get_valid_file_input("请输入结晶峰标定文件的路径(例如: F://1//test//Calibration.cake): ")
        start_id, end_id = get_valid_range_input("请输入文件读取的范围(例如: 5, 23): ")

        # Caliberate File Read
        q_caliberate, I_caliberate = single_file_reader(file_path=calibrate_file, start_value=start_id, end_value=end_id)
        I_normalized = smoothed_and_normalized(I_caliberate, 10, 3)
        baseline_corrector = BaselineCorrector(q_caliberate, I_normalized)

        # Caliberate Crystal Peaks
        selected_points_index = baseline_corrector.selected_points
        q1, q2 = q_caliberate[selected_points_index[0]], q_caliberate[selected_points_index[1]]
        I1, I2 = I_normalized[selected_points_index[0]], I_normalized[selected_points_index[1]]
        baseline = baseline_calculation(q_caliberate, q1, q2, I1, I2)
        baseline_corrected = [max(0, I_val - baseline_val) if selected_points_index[0] <= i <= selected_points_index[1] else I_val for i, (I_val, baseline_val) in enumerate(zip(I_normalized, baseline))]
        crystal_info = PeakInfo(q_caliberate, baseline_corrected)
        crystal_info.cry_peaks_find()

        # Caliberate Amorphous Peaks
        amorphous_params, amorphous_num = get_valid_amorphous_params()

        # Create Output Files
        header = [f"q = {loc}" for loc in amorphous_params[1::3]]
        for index in crystal_info.crystal_indices:
            header.append(f"q = {q_caliberate[index]}")

        with open("crystallinity.csv", mode='w', newline='') as file:
            csv_writer = csv.writer(file)        
            csv_writer.writerow(["file name"] + header)

        with open("grain_size.csv", mode='w', newline='') as file:
            csv_writer = csv.writer(file)
            csv_writer.writerow(["file name"] + header)

        # Process All Data
        dir_path = Path(PATH)
        files_list = [file.name for file in dir_path.iterdir() if file.is_file()]
        my_dict = files_reader(dir_path=PATH, start_value=start_id, end_value=end_id)
        filename_iterator = iter(files_list)
        q_iterator = iter(my_dict["x"])
        I_iterator = iter(my_dict["y"])

        while True:
            try:
                q_name = next(q_iterator)
                I_name = next(I_iterator)
                file_name = next(filename_iterator)
                q = my_dict["x"][q_name]
                I = my_dict["y"][I_name]
                main_1(file_name, q, I, amorphous_params, crystal_info.crystal_indices, selected_points_index)
            except StopIteration:
                break
    
    elif key == "取向因子":
        PATH = get_valid_directory_input("请输入要处理文件夹的路径(例如: F://1//test//): ")
        p = Path(PATH)
        file_list = list(p.glob("*.cake"))

        Pnc_list = []

        for idx, data_file in enumerate(file_list, 1):
            fai_list, I_list = [], []

            with open(data_file, 'r') as file:
                for line in itertools.islice(file, 6, None):
                    columns = line.strip().split(',')
                    fai_column = float(columns[0])
                    I_column = float(columns[1])
                    fai_list.append(fai_column)
                    I_list.append(I_column)
            
            fai_array = np.array(fai_list)
            I_array = np.array(I_list)

            nan_indices = np.isnan(I_list)
            fai = fai_array[~nan_indices]
            I = I_array[~nan_indices]
            I_smoothed = savgol_filter(I, 200, 2)
            
            peak_index = np.argmax(I_smoothed)
            position = fai[peak_index]
            start_index = np.argmin(np.abs(fai - (position-20)))
            end_index = np.argmin(np.abs(fai - (position+20)))
            fai_range = fai[start_index : end_index+1]
            I_range = I_smoothed[start_index : end_index+1]

            popt, _ = curve_fit(Maier_Saupe, fai_range, I_range, p0=[0, position, 0, 0], maxfev=9999999)
            A, φ0, β, I0 = popt
            I_fit = Maier_Saupe(fai, A, φ0, β, I0)
            Pnc = Pnc_calculatation(popt)
            
            plt.plot(fai, I_smoothed, ls='-', lw=3, c='blue', label='origin')
            plt.plot(fai_range, I_range, ls='--', lw=3, c='green', label='fit_range')
            plt.plot(position, I_smoothed[peak_index], 'ro')
            plt.plot(fai, I_fit, ls='--', lw=1, c='red', label='fit') 
            plt.legend(loc='upper right')
            plt.show()
            print(f"{idx}-{data_file}: {Pnc}")

            Pnc_list.append((data_file, Pnc))

        # 将结果保存到txt文件
        np.savetxt('Pnc.txt', Pnc_list, fmt='%s, %.6f', delimiter=',', header='index,Pnc', comments='')
        print("Pnc values have been saved to Pnc.txt file.")

    else:
        print("！警告，输入内容仅限'结晶度和晶粒尺寸'或'取向因子'，请重新输入.")
        key = input(f"请输入你想要执行的分析功能('结晶度和晶粒尺寸'，或'取向因子'): ")