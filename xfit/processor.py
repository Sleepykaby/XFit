from __future__ import annotations

from pathlib import Path
from typing import Callable

from .analysis import CrystallinityAnalyzer, OrientationAnalyzer
from .data_io import DataReader
from .models import AnalysisMode, FileResult, ProcessingConfig
from .preprocess import Preprocessor
from .reports import ReportManager

ProgressCallback = Callable[[str], None]


class BatchProcessor:
    def __init__(self, config: ProcessingConfig, progress_callback: ProgressCallback | None = None):
        self.config = config
        self.progress_callback = progress_callback or (lambda message: None)
        self.reader = DataReader()
        self.preprocessor = Preprocessor(config.smooth_window, config.smooth_poly_order)
        self.reports = ReportManager(config.output_dir)

    def run(self) -> list[FileResult]:
        files = self.reader.discover(self.config.input_path)
        results: list[FileResult] = []
        self.progress_callback(f"发现 {len(files)} 个文件，开始处理...")
        for index, file_path in enumerate(files, start=1):
            self.progress_callback(f"[{index}/{len(files)}] 正在处理 {file_path.name}")
            try:
                result = self._process_file(file_path)
            except Exception as error:  # noqa: BLE001 - must keep batch running and report failures.
                self.reports.write_error(file_path, error)
                result = FileResult(file_path=file_path, success=False, message=str(error))
            self.reports.write_result(result)
            results.append(result)
        self.progress_callback("处理完成。")
        return results

    def _process_file(self, file_path: Path) -> FileResult:
        series = self.reader.read(file_path, self.config.start_value, self.config.end_value)
        if self.config.mode == AnalysisMode.CRYSTALLINITY:
            return self._process_crystallinity(series.file_path, series.x, series.y)
        return self._process_orientation(series.file_path, series.x, series.y)

    def _process_crystallinity(self, file_path: Path, x_values: list[float], y_values: list[float]) -> FileResult:
        normalized = self.preprocessor.smooth_and_normalize(y_values)
        corrected = self.preprocessor.baseline_correct(x_values, normalized, self.config.baseline_indices)
        analyzer = CrystallinityAnalyzer()
        fit_result = analyzer.fit(x_values, corrected, self.config.amorphous_peaks, self.config.crystal_peak_indices)
        params = fit_result["params"]
        centers = [float(center) for center in params[1:-1:3]]
        metric_files = self.reports.save_crystallinity_metrics(file_path.name, centers, fit_result["areas"], fit_result["widths"])
        plot_file = self.reports.save_fit_plot(
            x_values,
            corrected,
            fit_result["fit"],
            fit_result["residuals"],
            fit_result["r_squared"],
            file_path.stem,
            "q (nm⁻¹)",
            self.config.crystal_peak_indices,
            fit_result["components"],
            params,
        )
        return FileResult(
            file_path=file_path,
            success=True,
            message="结晶度与晶粒尺寸分析完成",
            outputs=[plot_file, *metric_files],
            metrics={"r_squared": round(float(fit_result["r_squared"]), 6)},
        )

    def _process_orientation(self, file_path: Path, angle_values: list[float], y_values: list[float]) -> FileResult:
        filtered = self.preprocessor.filter_outliers(y_values, self.config.orientation_outlier_threshold)
        normalized = self.preprocessor.smooth_and_normalize(filtered.tolist())
        analyzer = OrientationAnalyzer()
        fit_result = analyzer.fit(angle_values, normalized)
        metric_file = self.reports.save_orientation_metric(file_path.name, float(fit_result["pnc"]))
        plot_file = self.reports.save_fit_plot(
            angle_values,
            normalized,
            fit_result["fit"],
            fit_result["residuals"],
            fit_result["r_squared"],
            file_path.stem,
            "φ (degree)",
            (int(fit_result["peak_index"]),),
            None,
            fit_result["params"],
        )
        return FileResult(
            file_path=file_path,
            success=True,
            message="取向因子分析完成",
            outputs=[plot_file, metric_file],
            metrics={"Pnc": round(float(fit_result["pnc"]), 6), "r_squared": round(float(fit_result["r_squared"]), 6)},
        )
