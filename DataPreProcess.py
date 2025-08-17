# -*- coding: utf-8 -*-
"================================================================================================================================="
from scipy.signal import savgol_filter
import matplotlib.pyplot as plt
from functools import partial
import numpy as np
from typing import Iterable, Union, List, Any


class BaselineCorrector:

    def __init__(self, x_data: List[float], y_data: List[float]):
        self.x_data = x_data
        self.y_data = y_data
        self.selected_points = []
        self.__initialize_plot()

    def __initialize_plot(self) -> None:
        _, self.__ax = plt.subplots()
        self.__ax.set_title('Baseline Correction')
        self.__lines = self.__ax.plot(self.x_data, self.y_data, ls='-', lw=0.7, label='Radius Integral')
        self.__ax.set_xlabel("q (nm$^{-1}$)", fontdict={'size': 12})
        self.__ax.set_ylabel("Intensity (a.u.)", fontdict={'size': 12})
        self.__vline = self.__ax.axvline(x=self.x_data[0], color='r', linestyle='-', lw=0.5)
        self.__hline = self.__ax.axhline(y=self.y_data[0], color='r', linestyle='-', lw=0.5)
        on_mouse_move_partial = partial(self.__on_mouse_move)
        plt.gcf().canvas.mpl_connect('motion_notify_event', on_mouse_move_partial)
        on_double_click_partial = partial(self.__on_double_click)
        plt.gcf().canvas.mpl_connect('button_press_event', on_double_click_partial)
        plt.show()

    def __on_mouse_move(self, event: Any) -> None:
        if event.inaxes:
            x = event.xdata
            y = event.ydata

            if x < min(self.x_data):
                x = min(self.x_data)
            elif x > max(self.x_data):
                x = max(self.x_data)

            index = min(range(len(self.x_data)), key=lambda i: abs(self.x_data[i] - x))
            x = self.x_data[index]
            y = self.y_data[index]

            self.__vline.set_xdata([x])
            self.__hline.set_ydata([y])
            plt.draw()

    def __on_double_click(self, event: Any) -> None:
        if event.dblclick:
            x = event.xdata
            index = min(range(len(self.x_data)), key=lambda i: abs(self.x_data[i] - x))
            self.selected_points.append(index)

            if len(self.selected_points) == 2:
                x1, x2 = self.x_data[self.selected_points[0]], self.x_data[self.selected_points[1]]
                y1, y2 = self.y_data[self.selected_points[0]], self.y_data[self.selected_points[1]]
                baseline = baseline_calculation(self.x_data, x1, x2, y1, y2)
                baseline_corrected = [max(0, y_val - baseline_val) if i >= self.selected_points[0] else y_val for i, (y_val, baseline_val) in enumerate(zip(self.y_data, baseline))]
                self.__lines[0].remove()
                self.__lines[0] = self.__ax.plot(self.x_data, baseline_corrected, ls='-', lw=1)[0]
                self.__ax.set_title('Baseline Correction')
                self.__reset_lines()

    def __reset_lines(self) -> None:
        if self.__vline:
            self.__vline.remove()
        if self.__hline:
            self.__hline.remove()

def smoothed_and_normalized(y_data: List[float], window_size: int, poly_order: int) -> Iterable[List]:
    y_smoothed = savgol_filter(y_data, window_size, poly_order)
    # y_normalized = [(y - min(y_smoothed)) / (max(y_smoothed)- min(y_smoothed)) for y in y_smoothed]
    return y_smoothed

def outlier_filter(data: List[float], threshold: Union[int, float]) -> List:
    filtered_data = [data[0]]
    for prev, current, nxt in zip(data, data[1:], data[2:]):
        if ((abs(current - prev) > threshold) or (abs(current - nxt) > threshold)) and ((current <= prev) or (current <= nxt)):
            filtered_data.append(np.nan)
        else:
            filtered_data.append(current)
    filtered_data.append(data[-1])
    return filtered_data

def baseline_calculation(x_data: List[float], x1: float, x2: float, y1: float, y2: float) -> Iterable[List]:
    line_slope = (y2 - y1) / (x2 - x1)
    line_intercept = y1 - line_slope * x1
    baseline = [line_slope * x_val + line_intercept for x_val in x_data]
    return baseline

"================================================================================================================================="
if __name__ == "__main__":
    from DirFileRead import file_reader

    PATH = "F://1//240222//1D_data//deducted_and_masked//PA6T_WAXS_T80//"
    my_dict = file_reader(dir_path=PATH, start_value=0, end_value=1)
    x_test = my_dict["q"][1]
    y_test = my_dict["I"][1]

    y_normalized = smoothed_and_normalized(y_test, 10, 3)

    baseline_corrector = BaselineCorrector(x_test, y_normalized)
    selected_points_index = baseline_corrector.selected_points
    x1, x2 = x_test[selected_points_index[0]], x_test[selected_points_index[1]]
    y1, y2 = y_normalized[selected_points_index[0]], y_normalized[selected_points_index[1]]
    baseline = baseline_calculation(x_test, x1, x2, y1, y2)
    baseline_corrected = [max(0, y_val - baseline_val) if selected_points_index[0] <= i <= selected_points_index[1] else y_val for i, (y_val, baseline_val) in enumerate(zip(y_normalized, baseline))]
    