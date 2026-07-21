from __future__ import annotations

from dataclasses import dataclass, field
from enum import Enum
from pathlib import Path
from typing import Iterable


class AnalysisMode(str, Enum):
    CRYSTALLINITY = "crystallinity"
    ORIENTATION = "orientation"


@dataclass(frozen=True)
class PeakParameter:
    amplitude: float
    center: float
    width: float

    def as_tuple(self) -> tuple[float, float, float]:
        return self.amplitude, self.center, self.width


@dataclass(frozen=True)
class ProcessingConfig:
    input_path: Path
    output_dir: Path
    mode: AnalysisMode
    start_value: float | None = None
    end_value: float | None = None
    calibration_file: Path | None = None
    baseline_indices: tuple[int, int] | None = None
    crystal_peak_indices: tuple[int, ...] = ()
    amorphous_peaks: tuple[PeakParameter, ...] = ()
    smooth_window: int = 11
    smooth_poly_order: int = 3
    orientation_outlier_threshold: float = 6.0


@dataclass(frozen=True)
class DataSeries:
    file_path: Path
    x: list[float]
    y: list[float]

    @property
    def name(self) -> str:
        return self.file_path.stem


@dataclass
class FileResult:
    file_path: Path
    success: bool
    message: str = ""
    outputs: list[Path] = field(default_factory=list)
    metrics: dict[str, float | str] = field(default_factory=dict)


def flatten_peaks(peaks: Iterable[PeakParameter]) -> list[float]:
    values: list[float] = []
    for peak in peaks:
        values.extend(peak.as_tuple())
    return values
