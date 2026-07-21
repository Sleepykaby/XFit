from __future__ import annotations

import numpy as np
from scipy.signal import savgol_filter


class Preprocessor:
    def __init__(self, smooth_window: int = 11, smooth_poly_order: int = 3):
        self.smooth_window = smooth_window
        self.smooth_poly_order = smooth_poly_order

    def smooth_and_normalize(self, y_values: list[float]) -> np.ndarray:
        y_array = np.asarray(y_values, dtype=float)
        window = self._valid_window(len(y_array))
        if window > self.smooth_poly_order:
            y_array = savgol_filter(y_array, window, self.smooth_poly_order)
        data_range = float(np.max(y_array) - np.min(y_array))
        if data_range == 0:
            return np.zeros_like(y_array)
        return (y_array - np.min(y_array)) / data_range

    def baseline_correct(self, x_values: list[float], y_values: np.ndarray, indices: tuple[int, int] | None) -> np.ndarray:
        if indices is None:
            return y_values
        start_index, end_index = sorted(indices)
        self._validate_index(start_index, x_values)
        self._validate_index(end_index, x_values)
        x1, x2 = x_values[start_index], x_values[end_index]
        y1, y2 = y_values[start_index], y_values[end_index]
        if x1 == x2:
            raise ValueError("基线两个点的 x 坐标不能相同")
        slope = (y2 - y1) / (x2 - x1)
        intercept = y1 - slope * x1
        baseline = slope * np.asarray(x_values, dtype=float) + intercept
        corrected = y_values.copy()
        corrected[start_index : end_index + 1] = np.maximum(0, corrected[start_index : end_index + 1] - baseline[start_index : end_index + 1])
        return corrected

    def filter_outliers(self, y_values: list[float], threshold: float) -> np.ndarray:
        data = np.asarray(y_values, dtype=float).copy()
        for index in range(1, len(data) - 1):
            prev_value, current_value, next_value = data[index - 1], data[index], data[index + 1]
            if (abs(current_value - prev_value) > threshold or abs(current_value - next_value) > threshold) and (
                current_value <= prev_value or current_value <= next_value
            ):
                data[index] = np.nan
        if np.isnan(data).any():
            valid = ~np.isnan(data)
            data[~valid] = np.interp(np.flatnonzero(~valid), np.flatnonzero(valid), data[valid])
        return data

    def _valid_window(self, data_length: int) -> int:
        if data_length <= self.smooth_poly_order + 2:
            return 0
        window = min(self.smooth_window, data_length if data_length % 2 else data_length - 1)
        if window % 2 == 0:
            window -= 1
        return max(window, self.smooth_poly_order + 2 + ((self.smooth_poly_order + 2) % 2 == 0))

    def _validate_index(self, index: int, values: list[float]) -> None:
        if index < 0 or index >= len(values):
            raise IndexError(f"索引超出数据范围: {index}")
