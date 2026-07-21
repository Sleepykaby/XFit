from __future__ import annotations

import numpy as np
from scipy import integrate
from scipy.optimize import curve_fit

from .models import PeakParameter, flatten_peaks


def gaussian(q: np.ndarray | float, amplitude: float, center: float, width: float) -> np.ndarray | float:
    return amplitude * np.exp(-((q - center) ** 2) / (2 * width**2))


class CrystallinityAnalyzer:
    def fit(
        self,
        x_values: list[float],
        y_values: np.ndarray,
        amorphous_peaks: tuple[PeakParameter, ...],
        crystal_peak_indices: tuple[int, ...],
    ) -> dict[str, object]:
        if not amorphous_peaks and not crystal_peak_indices:
            raise ValueError("至少需要设置一个非晶峰或晶体峰")

        x_array = np.asarray(x_values, dtype=float)
        y_array = np.asarray(y_values, dtype=float)
        initial, lower, upper = self._build_bounds(x_array, y_array, amorphous_peaks, crystal_peak_indices)
        params, _ = curve_fit(self.total_fit, x_array, y_array, p0=initial, bounds=(lower, upper), maxfev=999999)
        fitted = self.total_fit(x_array, *params)
        component_fits = self.single_fits(x_array, *params)
        residuals, r_squared = correlation_coefficient(y_array, fitted)
        areas = self.peak_areas(float(x_array[0]), float(x_array[-1]), params)
        return {
            "params": params,
            "fit": fitted,
            "components": component_fits,
            "residuals": residuals,
            "r_squared": r_squared,
            "areas": areas,
            "widths": [float(width) for width in params[2:-1:3]],
        }

    def total_fit(self, q: np.ndarray, *params: float) -> np.ndarray:
        components = self._components(q, params[:-1])
        return np.sum(components, axis=0) + params[-1]

    def single_fits(self, q: np.ndarray, *params: float) -> list[np.ndarray]:
        return [component + params[-1] for component in self._components(q, params[:-1])]

    def peak_areas(self, lower_bound: float, upper_bound: float, params: np.ndarray) -> list[float]:
        areas: list[float] = []
        for amplitude, center, width in np.asarray(params[:-1]).reshape(-1, 3):
            area, _ = integrate.quad(gaussian, lower_bound, upper_bound, args=(amplitude, center, width))
            areas.append(float(area))
        return areas

    def _components(self, q: np.ndarray, peak_params: tuple[float, ...] | np.ndarray) -> np.ndarray:
        return np.asarray([gaussian(q, amplitude, center, max(width, 1e-9)) for amplitude, center, width in np.asarray(peak_params).reshape(-1, 3)])

    def _build_bounds(
        self,
        x_values: np.ndarray,
        y_values: np.ndarray,
        amorphous_peaks: tuple[PeakParameter, ...],
        crystal_peak_indices: tuple[int, ...],
    ) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        initial = flatten_peaks(amorphous_peaks)
        lower: list[float] = []
        upper: list[float] = []

        for peak in amorphous_peaks:
            lower.extend([max(0, peak.amplitude - 5), peak.center - 1, max(1e-6, peak.width - 0.2)])
            upper.extend([peak.amplitude + 5, peak.center + 1, peak.width + 0.2])

        for index in crystal_peak_indices:
            if index < 0 or index >= len(x_values):
                raise IndexError(f"晶体峰索引超出范围: {index}")
            initial.extend([float(y_values[index]), float(x_values[index]), 0.02])
            lower.extend([0, float(x_values[index] - 1), 1e-6])
            upper.extend([20, float(x_values[index] + 1), 0.2])

        background = float(np.min(y_values) + 0.001)
        initial.append(background)
        lower.append(float(np.min(y_values)))
        upper.append(float(np.min(y_values) + 0.01))
        return np.asarray(initial), np.asarray(lower), np.asarray(upper)


class OrientationAnalyzer:
    def fit(self, angle_values: list[float], intensity_values: np.ndarray) -> dict[str, object]:
        x_array = np.asarray(angle_values, dtype=float)
        y_array = np.asarray(intensity_values, dtype=float)
        peak_index = int(np.argmax(y_array))
        initial = [float(y_array[peak_index]), float(x_array[peak_index]), 1.0, float(np.min(y_array))]
        params, _ = curve_fit(self.maier_saupe, x_array, y_array, p0=initial, maxfev=999999)
        fitted = self.maier_saupe(x_array, *params)
        residuals, r_squared = correlation_coefficient(y_array, fitted)
        return {
            "params": params,
            "fit": fitted,
            "components": None,
            "residuals": residuals,
            "r_squared": r_squared,
            "pnc": self.orientation_factor(params),
            "peak_index": peak_index,
        }

    def maier_saupe(self, x: np.ndarray | float, amplitude: float, center: float, beta: float, baseline: float) -> np.ndarray | float:
        return baseline + amplitude * np.exp(beta * (np.cos((x - center) * np.pi / 180)) ** 2)

    def orientation_factor(self, params: np.ndarray) -> float:
        beta = float(params[2])

        def numerator(x: float) -> float:
            return ((3 * (x**2) - 1) / 2) * np.exp(beta * (x**2))

        def denominator(x: float) -> float:
            return np.exp(beta * (x**2))

        value1, _ = integrate.quad(numerator, 0, 1)
        value2, _ = integrate.quad(denominator, 0, 1)
        return float(value1 / value2)


def correlation_coefficient(origin_data: np.ndarray, fit_value: np.ndarray) -> tuple[np.ndarray, float]:
    residuals = origin_data - fit_value
    rss = float(np.sum(residuals**2))
    tss = float(np.sum((origin_data - np.mean(origin_data)) ** 2))
    r_squared = 1.0 if tss == 0 else 1 - (rss / tss)
    return residuals, float(r_squared)
