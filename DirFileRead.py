# -*- coding: utf-8 -*-
"================================================================================================================================="
from pathlib import Path
from collections import defaultdict
from typing import Union, Iterable, Dict, List


def files_reader(dir_path: str, start_value: Union[int, float], end_value: Union[int, float]) -> Iterable[Dict]:
    p = Path(dir_path)
    data_dict = defaultdict(lambda: defaultdict(list))
    COUNT = 1

    for file_path in p.iterdir():
        list_generator = __line_reader(file_path, start_value, end_value)
        q_list, I_list = next(list_generator)
        data_dict["x"][COUNT] = q_list
        data_dict["y"][COUNT] = I_list
        data_dict["metadata"][COUNT] = file_path
        COUNT += 1

    print(f"{dir_path}内{COUNT+1}个文件已经读取完毕")
    return data_dict

def single_file_reader(file_path: str, start_value: Union[int, float], end_value: Union[int, float]) -> Iterable[Dict]:
    list_generator = __line_reader(file_path, start_value, end_value)
    q_list, I_list = next(list_generator)
    return q_list, I_list

def __line_reader(file_path: dir, start_value: Union[int, float], end_value: Union[int, float]) -> Iterable[List]:
    with open(file_path, 'r') as file:
        q_list = list()
        I_list = list()
        
        if str(file_path).endswith('.dat'):
            spiltter = '\t'
        elif str(file_path).endswith('.cake'):
            spiltter = ','

        for line in file:
            if line.__eq__("\n") or line.startswith('#'):
                continue
            num_alpha = sum(char.isalpha() for char in line)
            if num_alpha > 2:
                continue
            
            point_generator = __extractor(line, spiltter)
            q_point, I_point = next(point_generator) 
            
            if q_point is None and I_point is None:
                continue

            q_list.append(q_point)
            I_list.append(I_point)

        if start_value and end_value:
            start_id = min(range(len(q_list)), key=lambda i: abs(q_list[i] - start_value))
            end_id = min(range(len(q_list)), key=lambda i: abs(q_list[i] - end_value))
            q_sliced = q_list[start_id: end_id + 1]
            I_sliced = I_list[start_id: end_id + 1]

            yield q_sliced, I_sliced
        else:
            yield q_list, I_list

def __extractor(line, splitter=','):
    columns = line.strip().split(splitter)
    
    if len(columns) > 2:
        q = float(columns[0])
        I = float(columns[1])
        yield q, I
    
    else:
        yield None, None

"================================================================================================================================="
if __name__ == "__main__":
    PATH = "E://MyData//4_Project//2_intelligence molding//1D_data//deducted_and_masked//WAXD//" 
    my_dict = files_reader(dir_path=PATH, start_value=5, end_value=23)
    
    import matplotlib.pyplot as plt
    
    sequence = 1
    # print(my_dict["y"][sequence])
    print(my_dict["metadata"][9])

    plt.scatter(my_dict["x"][sequence], my_dict["y"][sequence], s=7, marker='o', c='None', edgecolor='red', label='I - q')
    plt.xlabel("q(nm$^{-1}$)")
    plt.ylabel("Intenisty(a.u.)")
    plt.legend(loc='upper right')
    plt.show()
