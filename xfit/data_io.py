from __future__ import annotations

import math
from pathlib import Path

from .models import DataSeries


class DataReader:
    """Read two-column WAXD data from .dat, .cake, .csv or .txt files."""

    SUPPORTED_SUFFIXES = {".dat", ".cake", ".csv", ".txt"}

    def discover(self, input_path: Path) -> list[Path]:
        path = Path(input_path).expanduser().resolve()
        if path.is_file():
            return [path]
        if not path.is_dir():
            raise FileNotFoundError(f"输入路径不存在: {path}")
        files = sorted(p for p in path.iterdir() if p.is_file() and p.suffix.lower() in self.SUPPORTED_SUFFIXES)
        if not files:
            raise FileNotFoundError(f"未找到可处理文件: {path}")
        return files

    def read(self, file_path: Path, start_value: float | None = None, end_value: float | None = None) -> DataSeries:
        x_values: list[float] = []
        y_values: list[float] = []
        splitter = self._detect_splitter(file_path)

        with Path(file_path).open("r", encoding="utf-8-sig", errors="ignore") as file:
            for line in file:
                point = self._parse_line(line, splitter, 3 if Path(file_path).suffix.lower() == ".cake" else 2)
                if point is None:
                    continue
                x_value, y_value = point
                x_values.append(x_value)
                y_values.append(y_value)

        if len(x_values) < 5:
            raise ValueError(f"有效数据点过少: {file_path}")

        x_values, y_values = self._slice_range(x_values, y_values, start_value, end_value)
        if len(x_values) < 5:
            raise ValueError("当前横坐标范围内有效数据点过少，请扩大范围。")
        return DataSeries(Path(file_path), x_values, y_values)

    def _detect_splitter(self, file_path: Path) -> str | None:
        if file_path.suffix.lower() == ".dat":
            return "\t"
        if file_path.suffix.lower() in {".cake", ".csv"}:
            return ","
        return None

    def _parse_line(self, line: str, splitter: str | None, min_columns: int = 2) -> tuple[float, float] | None:
        stripped = line.strip()
        if not stripped or stripped.startswith("#"):
            return None
        columns = stripped.split(splitter) if splitter else stripped.replace(",", " ").split()
        if len(columns) < min_columns:
            return None
        try:
            x_value = float(columns[0])
            y_value = float(columns[1])
        except ValueError:
            return None
        if not (math.isfinite(x_value) and math.isfinite(y_value)):
            return None
        return x_value, y_value

    def _slice_range(
        self,
        x_values: list[float],
        y_values: list[float],
        start_value: float | None,
        end_value: float | None,
    ) -> tuple[list[float], list[float]]:
        if start_value is None or end_value is None:
            return x_values, y_values
        start_index = min(range(len(x_values)), key=lambda index: abs(x_values[index] - start_value))
        end_index = min(range(len(x_values)), key=lambda index: abs(x_values[index] - end_value))
        if start_index > end_index:
            start_index, end_index = end_index, start_index
        return x_values[start_index : end_index + 1], y_values[start_index : end_index + 1]
