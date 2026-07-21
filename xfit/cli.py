from __future__ import annotations

import argparse
from pathlib import Path

from .models import AnalysisMode, PeakParameter, ProcessingConfig
from .processor import BatchProcessor


def parse_peaks(value: str) -> tuple[PeakParameter, ...]:
    if not value:
        return ()
    peaks: list[PeakParameter] = []
    for raw_peak in value.split(";"):
        amplitude, center, width = [float(part.strip()) for part in raw_peak.split(",")]
        peaks.append(PeakParameter(amplitude, center, width))
    return tuple(peaks)


def parse_indices(value: str | None) -> tuple[int, ...]:
    if not value:
        return ()
    return tuple(int(item.strip()) for item in value.split(",") if item.strip())


def main() -> None:
    parser = argparse.ArgumentParser(description="XFit WAXD batch analyzer")
    parser.add_argument("input_path", type=Path, help="输入文件或文件夹")
    parser.add_argument("output_dir", type=Path, help="输出目录")
    parser.add_argument("--mode", choices=[mode.value for mode in AnalysisMode], default=AnalysisMode.CRYSTALLINITY.value)
    parser.add_argument("--start", type=float, default=None, help="数据读取起点")
    parser.add_argument("--end", type=float, default=None, help="数据读取终点")
    parser.add_argument("--baseline", default=None, help="基线索引，例如 10,180")
    parser.add_argument("--crystal-peaks", default=None, help="晶体峰索引，例如 80,120")
    parser.add_argument("--amorphous-peaks", default="", help="非晶峰参数，例如 1.2,12.5,0.4;0.8,18,0.5")
    args = parser.parse_args()

    baseline_indices = parse_indices(args.baseline)
    config = ProcessingConfig(
        input_path=args.input_path,
        output_dir=args.output_dir,
        mode=AnalysisMode(args.mode),
        start_value=args.start,
        end_value=args.end,
        baseline_indices=(baseline_indices[0], baseline_indices[1]) if baseline_indices else None,
        crystal_peak_indices=parse_indices(args.crystal_peaks),
        amorphous_peaks=parse_peaks(args.amorphous_peaks),
    )
    results = BatchProcessor(config, print).run()
    failed = [result for result in results if not result.success]
    print(f"完成: 成功 {len(results) - len(failed)}/{len(results)}，输出目录 {config.output_dir}")


if __name__ == "__main__":
    main()
