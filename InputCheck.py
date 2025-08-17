# -*- coding: utf-8 -*-
"================================================================================================================================="
import os


def get_valid_directory_input(prompt):
    while True:
        directory_path = input(prompt).strip()
        if os.path.isdir(fr'{directory_path}'):
            return directory_path
        else:
            print("无效文件夹路径，请重新输入。")

def get_valid_file_input(prompt):
    while True:
        file_path = input(prompt).strip()
        if os.path.isfile(fr'{file_path}'):
            return file_path
        else:
            print("无效文件路径，请重新输入。")

def get_valid_range_input(prompt):
    while True:
        try:
            start_id, end_id = input(prompt).strip().split(',')
            start_id, end_id = int(start_id), int(end_id)
            return start_id, end_id
        except ValueError:
            print("输入无效，请输入两个整数，用逗号分隔。")

def get_valid_amorphous_params():
    num_peaks = __get_valid_integer_input("请输入非晶峰的个数: ")
    amorphous_params = __get_amorphous_parameters(num_peaks)
    return amorphous_params, num_peaks

def __get_valid_integer_input(prompt):
    while True:
        try:
            value = int(input(prompt).strip())
            return value
        except ValueError:
            print("输入无效，请输入一个整数。")

def __get_amorphous_parameters(num_peaks):
    amorphous_parameters = []

    for i in range(num_peaks):
        height = float(input(f"请输入第{i+1}个非晶峰的峰高: "))
        center = float(input(f"请输入第{i+1}个非晶峰的峰位: "))
        width = float(input(f"请输入第{i+1}个非晶峰的半高全宽(FHWM): "))

        amorphous_parameters.extend([height, center, width])

    return amorphous_parameters

"================================================================================================================================="
if __name__ == "__main__":
    string = '--*欢迎使用 WAXD Batch Processing version1.0 批处理软件,\n请按照提示进行输入,\n当前版本仅支持结晶度的批处理功能*--'
    print(string)
    print('=' * len(string))

    dir_path = get_valid_directory_input("请输入要处理文件夹的路径(例如: F://1//test//): ")
    file_path = get_valid_file_input("请输入要处理的文件路径(例如: F://1//test//Calibration.cake): ")

    start_id, end_id = get_valid_range_input("请输入文件读取的范围(例如: 5, 23): ")
    print(start_id)
    print(end_id)

    amorphous_params, amorphous_num = get_valid_amorphous_params()
