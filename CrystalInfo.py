# -*- coding: utf-8 -*-
"================================================================================================================================="
import matplotlib.pyplot as plt
from typing import *


class PeakInfo:
    def __init__(self, x_data: List[float], y_data: List[float]):
        self.x_data = x_data
        self.y_data = y_data
        self.crystal_indices = list()
        self.vline = None
        self.hline = None

    def on_mouse_move(self, event: Any):
        if event.inaxes:
            x = min(max(event.xdata, min(self.x_data)), max(self.x_data))
            index = min(range(len(self.x_data)), key=lambda i: abs(self.x_data[i] - x))
            x, y = self.x_data[index], self.y_data[index]

            self.vline.set_xdata([x])
            self.hline.set_ydata([y])
            plt.draw()

    def on_double_click(self, event: Any):
        if event.dblclick:
            x, y = event.xdata, event.ydata
            index = min(range(len(self.x_data)), key=lambda i: abs(self.x_data[i] - x))
            print(f"coordinate: (x={x}, y={y}), index: {index}")
            self.crystal_indices.append(index)

    def cry_peaks_find(self):
        plt.plot(self.x_data, self.y_data, ls='-', c='black', lw=0.5, label='Radius Intergral')
        plt.title("Select Crystal Peak Position by Double Click")
        plt.xlabel("q (nm$^{-1}$)", fontdict={'size': 12})
        plt.ylabel("Intensity (a.u.)", fontdict={'size': 12})
        self.vline = plt.axvline(x=self.x_data[0], c='r', ls='-', lw=0.5)
        self.hline = plt.axhline(y=self.y_data[0], c='r', ls='-', lw=0.5)

        plt.gcf().canvas.mpl_connect('motion_notify_event', self.on_mouse_move)
        plt.gcf().canvas.mpl_connect('button_press_event', self.on_double_click)
        plt.show()

"================================================================================================================================="
if __name__ == "__main__":
    import numpy as np
    
    xdata = np.linspace(0, 10, 100)
    ydata = np.sin(xdata) + np.random.normal(0, 0.1, size=100)

    peak_info = PeakInfo(xdata, ydata)
    peak_info.cry_peaks_find()
    
    print(peak_info.crystal_indices)