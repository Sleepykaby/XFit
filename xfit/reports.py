from __future__ import annotations

import csv
import logging
import traceback
from datetime import datetime
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.cm as cm
import matplotlib.pyplot as plt
import numpy as np

from .models import FileResult


class ReportManager:
    def __init__(self, output_dir: Path):
        self.output_dir = Path(output_dir).expanduser().resolve()
        self.output_dir.mkdir(parents=True, exist_ok=True)
        self.plots_dir = self.output_dir / "plots"
        self.plots_dir.mkdir(exist_ok=True)
        self.logs_dir = self.output_dir / "logs"
        self.logs_dir.mkdir(exist_ok=True)
        self.error_csv = self.output_dir / "error_report.csv"
        self.summary_csv = self.output_dir / "summary.csv"
        self.logger = self._create_logger()
        self._init_csv_files()

    def write_result(self, result: FileResult) -> None:
        with self.summary_csv.open("a", newline="", encoding="utf-8-sig") as file:
            writer = csv.writer(file)
            writer.writerow([
                result.file_path.name,
                "成功" if result.success else "失败",
                result.message,
                "; ".join(f"{key}={value}" for key, value in result.metrics.items()),
                "; ".join(str(path) for path in result.outputs),
            ])

    def write_error(self, file_path: Path, error: BaseException) -> None:
        self.logger.exception("处理文件失败: %s", file_path)
        with self.error_csv.open("a", newline="", encoding="utf-8-sig") as file:
            writer = csv.writer(file)
            writer.writerow([datetime.now().isoformat(timespec="seconds"), file_path, type(error).__name__, str(error), traceback.format_exc()])

    def save_crystallinity_metrics(self, file_name: str, centers: list[float], areas: list[float], widths: list[float]) -> tuple[Path, Path]:
        crystallinity_file = self.output_dir / "crystallinity.csv"
        grain_file = self.output_dir / "grain_size.csv"
        self._ensure_metric_header(crystallinity_file, centers, "area")
        self._ensure_metric_header(grain_file, centers, "width")
        with crystallinity_file.open("a", newline="", encoding="utf-8-sig") as file:
            csv.writer(file).writerow([file_name, *areas])
        with grain_file.open("a", newline="", encoding="utf-8-sig") as file:
            csv.writer(file).writerow([file_name, *widths])
        return crystallinity_file, grain_file

    def save_orientation_metric(self, file_name: str, pnc: float) -> Path:
        output = self.output_dir / "orientation_factor.csv"
        if not output.exists():
            with output.open("w", newline="", encoding="utf-8-sig") as file:
                csv.writer(file).writerow(["file_name", "Pnc"])
        with output.open("a", newline="", encoding="utf-8-sig") as file:
            csv.writer(file).writerow([file_name, pnc])
        return output

    def save_fit_plot(
        self,
        x_values: list[float],
        y_values: np.ndarray,
        fitted: np.ndarray,
        residuals: np.ndarray,
        r_squared: float,
        file_stem: str,
        x_label: str,
        peak_indices: tuple[int, ...] = (),
        component_fits: list[np.ndarray] | None = None,
        params: np.ndarray | None = None,
    ) -> Path:
        plot_file = self.plots_dir / f"{file_stem}.png"
        x_array = np.asarray(x_values, dtype=float)
        fig, (fit_axis, residual_axis) = plt.subplots(nrows=2, figsize=(9, 7), gridspec_kw={"height_ratios": [3, 1.5]})
        fig.subplots_adjust(hspace=0.35)
        fit_axis.plot(x_array, y_values, color="#667085", linewidth=1, label="Data")
        fit_axis.plot(x_array, fitted, "--", color="#d92d20", linewidth=1.2, label="Fit")
        fit_axis.set_xlabel(x_label)
        fit_axis.set_ylabel("Intensity (a.u.)")
        fit_axis.grid(True, alpha=0.25)
        fit_axis.legend(loc="upper left")

        if component_fits:
            baseline = np.full_like(x_array, fill_value=float(params[-1]) if params is not None else 0)
            for index, component in enumerate(component_fits):
                color = cm.viridis(index / max(1, len(component_fits) - 1))
                fit_axis.plot(x_array, component, "--", color=color, linewidth=0.9)
                fit_axis.fill_between(x_array, component, baseline, color=color, alpha=0.22)
        for index in peak_indices:
            if 0 <= index < len(x_array):
                fit_axis.plot(x_array[index], y_values[index], "o", color="#175cd3")

        residual_axis.scatter(x_array, residuals, s=10, facecolor="white", edgecolor="#344054")
        residual_axis.axhline(0, color="#d92d20", linestyle="--", linewidth=0.8)
        residual_axis.set_xlabel(x_label)
        residual_axis.set_ylabel("Residuals")
        residual_axis.text(0.02, 0.88, f"R² = {r_squared:.5f}", transform=residual_axis.transAxes, fontsize=10)
        fig.savefig(plot_file, dpi=180, bbox_inches="tight")
        plt.close(fig)
        return plot_file

    def _create_logger(self) -> logging.Logger:
        logger = logging.getLogger(f"xfit.{id(self)}")
        logger.setLevel(logging.INFO)
        logger.handlers.clear()
        handler = logging.FileHandler(self.logs_dir / "xfit.log", encoding="utf-8")
        handler.setFormatter(logging.Formatter("%(asctime)s [%(levelname)s] %(message)s"))
        logger.addHandler(handler)
        return logger

    def _init_csv_files(self) -> None:
        with self.summary_csv.open("w", newline="", encoding="utf-8-sig") as file:
            csv.writer(file).writerow(["file_name", "status", "message", "metrics", "outputs"])
        with self.error_csv.open("w", newline="", encoding="utf-8-sig") as file:
            csv.writer(file).writerow(["time", "file", "error_type", "message", "traceback"])

    def _ensure_metric_header(self, output: Path, centers: list[float], metric_name: str) -> None:
        if output.exists():
            return
        with output.open("w", newline="", encoding="utf-8-sig") as file:
            csv.writer(file).writerow(["file_name", *[f"q={center:.6g}_{metric_name}" for center in centers]])
