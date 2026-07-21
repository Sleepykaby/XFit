from __future__ import annotations

import sys
from dataclasses import dataclass, field
from pathlib import Path
from typing import Callable

import matplotlib as mpl
import numpy as np
from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure
from PyQt6.QtCore import QThread, Qt, pyqtSignal
from PyQt6.QtGui import QFont
from PyQt6.QtWidgets import (
    QApplication,
    QButtonGroup,
    QCheckBox,
    QComboBox,
    QDoubleSpinBox,
    QFileDialog,
    QFrame,
    QGridLayout,
    QGroupBox,
    QHBoxLayout,
    QLabel,
    QLineEdit,
    QListWidget,
    QMessageBox,
    QPushButton,
    QProgressBar,
    QRadioButton,
    QSlider,
    QSpinBox,
    QStackedWidget,
    QTextEdit,
    QVBoxLayout,
    QWidget,
)

from .analysis import gaussian
from .data_io import DataReader
from .models import AnalysisMode, DataSeries, PeakParameter, ProcessingConfig
from .preprocess import Preprocessor
from .processor import BatchProcessor

mpl.rcParams["font.sans-serif"] = ["Microsoft YaHei", "SimHei", "Arial Unicode MS", "DejaVu Sans"]
mpl.rcParams["axes.unicode_minus"] = False


@dataclass
class WizardState:
    input_path: Path | None = None
    output_dir: Path | None = None
    mode: AnalysisMode = AnalysisMode.CRYSTALLINITY
    smooth_window: int = 11
    smooth_poly_order: int = 3
    outlier_threshold: float = 6.0
    calibration_file: Path | None = None
    calibration_raw: DataSeries | None = None
    calibration_sliced: DataSeries | None = None
    processed_x: list[float] = field(default_factory=list)
    processed_y: np.ndarray | None = None
    start_value: float | None = None
    end_value: float | None = None
    baseline_indices: tuple[int, int] | None = None
    crystal_peak_indices: tuple[int, ...] = ()
    amorphous_peaks: tuple[PeakParameter, ...] = ()


class ProcessingWorker(QThread):
    progress = pyqtSignal(str)
    finished = pyqtSignal(int, int, str)
    failed = pyqtSignal(str)

    def __init__(self, config: ProcessingConfig):
        super().__init__()
        self.config = config

    def run(self) -> None:
        try:
            results = BatchProcessor(self.config, self.progress.emit).run()
            success_count = sum(result.success for result in results)
            self.finished.emit(success_count, len(results), str(self.config.output_dir))
        except Exception as error:  # noqa: BLE001
            self.failed.emit(str(error))


class PlotPanel(QWidget):
    point_double_clicked = pyqtSignal(float)

    def __init__(self):
        super().__init__()
        self.figure = Figure(figsize=(7, 4), facecolor="white")
        self.canvas = FigureCanvas(self.figure)
        self.axis = self.figure.add_subplot(111)
        layout = QVBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.addWidget(self.canvas)
        self.canvas.mpl_connect("button_press_event", self._on_click)

    def plot_curve(
        self,
        x_values,
        y_values,
        title="Curve",
        x_label="q (nm^-1)",
        markers=(),
        amorphous_peaks=(),
        show_legend=False,
    ):
        self.axis.clear()
        x_array = np.asarray(x_values, dtype=float)
        y_array = np.asarray(y_values, dtype=float)
        if x_array.size and y_array.size:
            self.axis.plot(x_array, y_array, color="#175CD3", linewidth=1.5, label="data")
        for index in markers or ():
            if 0 <= index < len(x_array):
                self.axis.scatter(x_array[index], y_array[index], color="#D92D20", s=48, zorder=4)
        for peak_index, peak in enumerate(amorphous_peaks or (), start=1):
            component = gaussian(x_array, peak.amplitude, peak.center, peak.width)
            self.axis.plot(x_array, component, "--", linewidth=1.1, label=f"amorphous {peak_index}")
        if show_legend or amorphous_peaks:
            self.axis.legend(loc="upper right")
        self.axis.set_title(title, fontsize=12, pad=12)
        self.axis.set_xlabel(x_label)
        self.axis.set_ylabel("Intensity (a.u.)")
        self.axis.grid(True, alpha=0.25)
        self.figure.tight_layout()
        self.canvas.draw_idle()

    def clear(self, title="Waiting for data"):
        self.axis.clear()
        self.axis.set_title(title, fontsize=12, pad=12)
        self.axis.set_xlabel("q (nm^-1)")
        self.axis.set_ylabel("Intensity (a.u.)")
        self.axis.grid(True, alpha=0.2)
        self.canvas.draw_idle()

    def _on_click(self, event):
        if event.dblclick and event.inaxes and event.xdata is not None:
            self.point_double_clicked.emit(float(event.xdata))


def show_error(parent, title, message):
    QMessageBox.critical(parent, title, message)


def show_info(parent, title, message):
    QMessageBox.information(parent, title, message)


def make_card() -> QFrame:
    frame = QFrame()
    frame.setObjectName("Card")
    frame.setFrameShape(QFrame.Shape.StyledPanel)
    return frame


def section_title(title: str, subtitle: str) -> QWidget:
    widget = QWidget()
    layout = QVBoxLayout(widget)
    layout.setContentsMargins(0, 0, 0, 0)
    title_label = QLabel(title)
    title_label.setStyleSheet("font-size: 18px; font-weight: 800; color: #101828;")
    subtitle_label = QLabel(subtitle)
    subtitle_label.setStyleSheet("color: #667085;")
    subtitle_label.setWordWrap(True)
    layout.addWidget(title_label)
    layout.addWidget(subtitle_label)
    return widget


def primary_button(text: str) -> QPushButton:
    button = QPushButton(text)
    button.setObjectName("PrimaryButton")
    return button


def danger_button(text: str) -> QPushButton:
    button = QPushButton(text)
    button.setObjectName("DangerButton")
    return button


def path_row(label: str, edit: QLineEdit, first_callback: Callable[[], None], second_callback: Callable[[], None] | None) -> QVBoxLayout:
    layout = QVBoxLayout()
    layout.addWidget(QLabel(label))
    row = QHBoxLayout()
    row.addWidget(edit, stretch=1)
    first_button = QPushButton("选择文件" if second_callback else "选择")
    first_button.clicked.connect(first_callback)
    row.addWidget(first_button)
    if second_callback:
        second_button = QPushButton("选择文件夹")
        second_button.clicked.connect(second_callback)
        row.addWidget(second_button)
    layout.addLayout(row)
    return layout


def integer_spin(minimum: int, maximum: int, value: int, step: int) -> QSpinBox:
    spin = QSpinBox()
    spin.setRange(minimum, maximum)
    spin.setSingleStep(step)
    spin.setValue(value)
    return spin


def double_spin(minimum: float, maximum: float, value: float, decimals: int) -> QDoubleSpinBox:
    spin = QDoubleSpinBox()
    spin.setDecimals(decimals)
    spin.setRange(minimum, maximum)
    spin.setValue(value)
    spin.setSingleStep(0.1)
    return spin


def nearest_index(values: list[float], target: float) -> int:
    return min(range(len(values)), key=lambda index: abs(values[index] - target))


def parse_amorphous_text(text: str) -> tuple[PeakParameter, ...]:
    peaks: list[PeakParameter] = []
    for line_number, raw_line in enumerate(text.splitlines(), start=1):
        line = raw_line.strip()
        if not line:
            continue
        parts = [part.strip() for part in line.replace(";", ",").split(",") if part.strip()]
        if len(parts) != 3:
            raise ValueError(f"第 {line_number} 行必须包含三个数值：峰高,峰位置,峰宽。")
        try:
            amplitude, center, width = (float(parts[0]), float(parts[1]), float(parts[2]))
        except ValueError as error:
            raise ValueError(f"第 {line_number} 行存在非数字内容。") from error
        if amplitude <= 0 or center < 0 or width <= 0:
            raise ValueError(f"第 {line_number} 行参数不合理。峰高和峰宽必须大于 0，峰位置不能为负。")
        peaks.append(PeakParameter(amplitude, center, width))
    return tuple(peaks)

class XFitWindow(QWidget):
    def __init__(self):
        super().__init__()
        self.state = WizardState()
        self.worker: ProcessingWorker | None = None
        self.setWindowTitle("XFit WAXD Batch Analyzer")
        self.resize(1220, 800)
        self.setMinimumSize(1080, 720)
        self._build_ui()
        self._apply_style()

    def _build_ui(self):
        root = QVBoxLayout(self)
        root.setContentsMargins(28, 24, 28, 24)
        header = QHBoxLayout()
        title_box = QVBoxLayout()
        title = QLabel("XFit WAXD 批量分析")
        title.setObjectName("Title")
        subtitle = QLabel("向导式流程 | 标定曲线预览 | 交互选峰 | 自动错误报告")
        subtitle.setObjectName("Subtitle")
        title_box.addWidget(title)
        title_box.addWidget(subtitle)
        header.addLayout(title_box)
        header.addStretch()
        self.step_label = QLabel("步骤 1 / 4")
        self.step_label.setObjectName("StepBadge")
        header.addWidget(self.step_label)
        root.addLayout(header)
        self.stack = QStackedWidget()
        root.addWidget(self.stack, stretch=1)
        self.start_page = StartPage(self.state, self.go_next_from_start)
        self.calibration_page = CalibrationPage(self.state, self.go_start, self.go_peak_page)
        self.peak_page = PeakPage(self.state, self.go_calibration_page, self.run_crystallinity)
        self.orientation_page = OrientationPage(self.state, self.go_start, self.run_orientation)
        self.result_page = ResultPage(self.go_start)
        for page in [self.start_page, self.calibration_page, self.peak_page, self.orientation_page, self.result_page]:
            self.stack.addWidget(page)
        self.go_start()

    def _apply_style(self):
        self.setFont(QFont("Microsoft YaHei UI", 10))
        self.setStyleSheet("""
        QWidget { background: #F6F8FC; color: #101828; }
        QLabel#Title { font-size: 25px; font-weight: 800; color: #101828; }
        QLabel#Subtitle { color: #667085; font-size: 13px; }
        QLabel#StepBadge { background: #E0EAFF; color: #175CD3; padding: 8px 14px; border-radius: 14px; font-weight: 700; }
        QFrame#Card, QGroupBox { background: #FFFFFF; border: 1px solid #EAECF0; border-radius: 18px; }
        QGroupBox { margin-top: 12px; padding: 18px 14px 14px 14px; font-weight: 700; }
        QGroupBox::title { subcontrol-origin: margin; left: 16px; padding: 0 8px; color: #344054; }
        QLineEdit, QDoubleSpinBox, QSpinBox, QComboBox, QTextEdit, QListWidget { background: #FFFFFF; border: 1px solid #D0D5DD; border-radius: 10px; padding: 8px; }
        QLineEdit:focus, QDoubleSpinBox:focus, QSpinBox:focus, QTextEdit:focus { border: 1px solid #2E90FA; }
        QPushButton { background: #FFFFFF; border: 1px solid #D0D5DD; border-radius: 10px; padding: 9px 14px; font-weight: 600; }
        QPushButton:hover { background: #F2F4F7; }
        QPushButton#PrimaryButton { background: #175CD3; border-color: #175CD3; color: #FFFFFF; }
        QPushButton#PrimaryButton:hover { background: #1849A9; }
        QPushButton#DangerButton { color: #B42318; border-color: #FDA29B; }
        QProgressBar { border: 1px solid #D0D5DD; border-radius: 8px; height: 12px; background: #FFFFFF; }
        QProgressBar::chunk { border-radius: 8px; background: #2E90FA; }
        """)

    def go_start(self):
        self._set_page(self.start_page, "步骤 1 / 4")

    def go_next_from_start(self):
        if self.state.mode == AnalysisMode.CRYSTALLINITY:
            self.calibration_page.refresh_from_state()
            self._set_page(self.calibration_page, "结晶分析 2 / 4")
        else:
            self.orientation_page.refresh_from_state()
            self._set_page(self.orientation_page, "取向分析 2 / 3")

    def go_calibration_page(self):
        self.calibration_page.refresh_from_state()
        self._set_page(self.calibration_page, "结晶分析 2 / 4")

    def go_peak_page(self):
        self.peak_page.refresh_from_state()
        self._set_page(self.peak_page, "结晶分析 3 / 4")

    def _set_page(self, page, step_text):
        self.step_label.setText(step_text)
        self.stack.setCurrentWidget(page)

    def run_crystallinity(self):
        self._run_worker(ProcessingConfig(
            input_path=self.state.input_path,
            output_dir=self.state.output_dir,
            mode=AnalysisMode.CRYSTALLINITY,
            start_value=self.state.start_value,
            end_value=self.state.end_value,
            calibration_file=self.state.calibration_file,
            baseline_indices=self.state.baseline_indices,
            crystal_peak_indices=self.state.crystal_peak_indices,
            amorphous_peaks=self.state.amorphous_peaks,
            smooth_window=self.state.smooth_window,
            smooth_poly_order=self.state.smooth_poly_order,
            orientation_outlier_threshold=self.state.outlier_threshold,
        ), "开始结晶度 / 晶粒尺寸批处理...")

    def run_orientation(self):
        self._run_worker(ProcessingConfig(
            input_path=self.state.input_path,
            output_dir=self.state.output_dir,
            mode=AnalysisMode.ORIENTATION,
            start_value=self.state.start_value,
            end_value=self.state.end_value,
            smooth_window=self.state.smooth_window,
            smooth_poly_order=self.state.smooth_poly_order,
            orientation_outlier_threshold=self.state.outlier_threshold,
        ), "开始取向因子批处理...")

    def _run_worker(self, config, message):
        self.result_page.reset(message)
        self._set_page(self.result_page, "结果")
        self.worker = ProcessingWorker(config)
        self.worker.progress.connect(self.result_page.append_log)
        self.worker.finished.connect(self._on_worker_finished)
        self.worker.failed.connect(self._on_worker_failed)
        self.worker.start()

    def _on_worker_finished(self, success_count, total_count, output_dir):
        self.result_page.finish(f"处理完成：成功 {success_count}/{total_count}\n输出目录：{output_dir}")

    def _on_worker_failed(self, message):
        self.result_page.fail(f"处理失败：{message}")


class StartPage(QWidget):
    def __init__(self, state, next_callback):
        super().__init__()
        self.state = state
        self.next_callback = next_callback
        self._build()

    def _build(self):
        layout = QVBoxLayout(self)
        card = make_card()
        card_layout = QVBoxLayout(card)
        card_layout.addWidget(section_title("基础设置", "选择输入数据、输出目录和分析功能。"))
        self.input_edit = QLineEdit()
        self.output_edit = QLineEdit(str(Path.cwd() / "xfit_output"))
        card_layout.addLayout(path_row("输入文件或文件夹", self.input_edit, self.choose_input_file, self.choose_input_folder))
        card_layout.addLayout(path_row("输出文件夹", self.output_edit, self.choose_output_dir, None))
        mode_group = QGroupBox("分析功能")
        mode_layout = QHBoxLayout(mode_group)
        self.crystallinity_radio = QRadioButton("结晶度 / 晶粒尺寸")
        self.orientation_radio = QRadioButton("取向分析")
        self.crystallinity_radio.setChecked(True)
        self.mode_buttons = QButtonGroup(self)
        self.mode_buttons.addButton(self.crystallinity_radio)
        self.mode_buttons.addButton(self.orientation_radio)
        mode_layout.addWidget(self.crystallinity_radio)
        mode_layout.addWidget(self.orientation_radio)
        mode_layout.addStretch()
        card_layout.addWidget(mode_group)
        options = QGroupBox("通用参数")
        options_layout = QGridLayout(options)
        self.window_spin = integer_spin(5, 101, 11, 2)
        self.poly_spin = integer_spin(1, 9, 3, 1)
        self.outlier_spin = double_spin(0.1, 1000, 6.0, 2)
        options_layout.addWidget(QLabel("平滑窗口"), 0, 0)
        options_layout.addWidget(self.window_spin, 0, 1)
        options_layout.addWidget(QLabel("多项式阶数"), 0, 2)
        options_layout.addWidget(self.poly_spin, 0, 3)
        options_layout.addWidget(QLabel("取向异常值阈值"), 1, 0)
        options_layout.addWidget(self.outlier_spin, 1, 1)
        card_layout.addWidget(options)
        actions = QHBoxLayout()
        actions.addStretch()
        next_button = primary_button("下一步")
        next_button.clicked.connect(self.validate_and_next)
        actions.addWidget(next_button)
        card_layout.addLayout(actions)
        layout.addWidget(card)
        layout.addStretch()

    def choose_input_file(self):
        file_path, _ = QFileDialog.getOpenFileName(self, "选择输入文件", "", "Data Files (*.dat *.cake *.csv *.txt);;All Files (*)")
        if file_path:
            self.input_edit.setText(file_path)

    def choose_input_folder(self):
        folder = QFileDialog.getExistingDirectory(self, "选择输入文件夹")
        if folder:
            self.input_edit.setText(folder)

    def choose_output_dir(self):
        folder = QFileDialog.getExistingDirectory(self, "选择输出文件夹")
        if folder:
            self.output_edit.setText(folder)

    def validate_and_next(self):
        input_path = Path(self.input_edit.text().strip())
        output_dir_text = self.output_edit.text().strip()
        if not input_path.exists():
            show_error(self, "输入路径错误", "请选择存在的输入文件或文件夹。")
            return
        if not output_dir_text:
            show_error(self, "输出路径错误", "请指定输出文件夹。")
            return
        if self.window_spin.value() <= self.poly_spin.value():
            show_error(self, "平滑参数错误", "平滑窗口必须大于多项式阶数。")
            return
        if int(self.window_spin.value()) % 2 == 0:
            show_error(self, "平滑参数错误", "平滑窗口必须是奇数。")
            return
        self.state.input_path = input_path
        self.state.output_dir = Path(output_dir_text)
        self.state.mode = AnalysisMode.CRYSTALLINITY if self.crystallinity_radio.isChecked() else AnalysisMode.ORIENTATION
        self.state.smooth_window = int(self.window_spin.value())
        self.state.smooth_poly_order = int(self.poly_spin.value())
        self.state.outlier_threshold = float(self.outlier_spin.value())
        self.next_callback()

class CalibrationPage(QWidget):
    def __init__(self, state, back_callback, next_callback):
        super().__init__()
        self.state = state
        self.back_callback = back_callback
        self.next_callback = next_callback
        self.reader = DataReader()
        self.selected_baseline: list[int] = []
        self._syncing_range = False
        self._build()

    def _build(self):
        layout = QVBoxLayout(self)
        main = QHBoxLayout()
        control_card = make_card()
        controls = QVBoxLayout(control_card)
        controls.addWidget(section_title("标定曲线处理", "加载标定文件，设置横坐标范围，可选基线校准。"))
        self.calibration_edit = QLineEdit()
        controls.addLayout(path_row("标定文件", self.calibration_edit, self.choose_calibration_file, None))
        range_box = QGroupBox("横坐标范围")
        range_layout = QGridLayout(range_box)
        self.start_spin = double_spin(0.0, 1_000_000, 0.0, 4)
        self.end_spin = double_spin(0.0, 1_000_000, 30.0, 4)
        self.start_slider = QSlider(Qt.Orientation.Horizontal)
        self.end_slider = QSlider(Qt.Orientation.Horizontal)
        for slider in [self.start_slider, self.end_slider]:
            slider.setRange(0, 10000)
        range_layout.addWidget(QLabel("起点"), 0, 0)
        range_layout.addWidget(self.start_spin, 0, 1)
        range_layout.addWidget(self.start_slider, 1, 0, 1, 2)
        range_layout.addWidget(QLabel("终点"), 2, 0)
        range_layout.addWidget(self.end_spin, 2, 1)
        range_layout.addWidget(self.end_slider, 3, 0, 1, 2)
        controls.addWidget(range_box)
        self.baseline_check = QCheckBox("开启基线校准模式：在图中双击选择两个基线点")
        self.baseline_check.stateChanged.connect(self.toggle_baseline_mode)
        controls.addWidget(self.baseline_check)
        redraw_button = QPushButton("重新绘制 / 检查处理效果")
        redraw_button.clicked.connect(self.apply_range_and_baseline)
        reset_button = danger_button("重置当前设置")
        reset_button.clicked.connect(self.reset_calibration)
        controls.addWidget(redraw_button)
        controls.addWidget(reset_button)
        controls.addStretch()
        main.addWidget(control_card, 1)
        plot_card = make_card()
        plot_layout = QVBoxLayout(plot_card)
        self.plot = PlotPanel()
        self.plot.point_double_clicked.connect(self.select_baseline_point)
        self.baseline_list = QListWidget()
        self.baseline_list.setMaximumHeight(76)
        plot_layout.addWidget(self.plot, stretch=1)
        plot_layout.addWidget(QLabel("已选基线点"))
        plot_layout.addWidget(self.baseline_list)
        main.addWidget(plot_card, 2)
        layout.addLayout(main, stretch=1)
        actions = QHBoxLayout()
        back_button = QPushButton("上一步")
        back_button.clicked.connect(self.back_callback)
        next_button = primary_button("下一步：选择结晶峰")
        next_button.clicked.connect(self.validate_and_next)
        actions.addWidget(back_button)
        actions.addStretch()
        actions.addWidget(next_button)
        layout.addLayout(actions)
        self.start_spin.valueChanged.connect(self._spin_range_changed)
        self.end_spin.valueChanged.connect(self._spin_range_changed)
        self.start_slider.valueChanged.connect(lambda value: self._slider_changed(value, self.start_spin))
        self.end_slider.valueChanged.connect(lambda value: self._slider_changed(value, self.end_spin))
        self.plot.clear("请选择标定文件")

    def refresh_from_state(self):
        if self.state.calibration_file:
            self.calibration_edit.setText(str(self.state.calibration_file))
        if self.state.calibration_raw:
            self._configure_range_widgets(self.state.calibration_raw.x)
            if self.state.processed_y is not None and self.state.processed_x:
                self.plot.plot_curve(self.state.processed_x, self.state.processed_y, "处理后的标定曲线", markers=self.state.baseline_indices or ())
            else:
                self.plot.plot_curve(self.state.calibration_raw.x, self.state.calibration_raw.y, "标定文件原始曲线")

    def choose_calibration_file(self):
        file_path, _ = QFileDialog.getOpenFileName(self, "选择标定文件", "", "Data Files (*.dat *.cake *.csv *.txt);;All Files (*)")
        if not file_path:
            return
        try:
            series = self.reader.read(Path(file_path), None, None)
        except Exception as error:
            show_error(self, "标定文件读取失败", str(error))
            return
        if min(series.x) < 0:
            show_error(self, "标定文件不合理", "横坐标必须大于等于 0。")
            return
        self.state.calibration_file = Path(file_path)
        self.state.calibration_raw = series
        self.state.calibration_sliced = None
        self.state.baseline_indices = None
        self.state.processed_x = []
        self.state.processed_y = None
        self.selected_baseline.clear()
        self.baseline_list.clear()
        self.calibration_edit.setText(file_path)
        self._configure_range_widgets(series.x)
        self.plot.plot_curve(series.x, series.y, "标定文件原始曲线")

    def toggle_baseline_mode(self):
        if self.baseline_check.isChecked():
            if self.state.calibration_sliced is None:
                self.apply_range_and_baseline(show_popup=False)
            show_info(self, "基线校准", "请在绘图区域双击选择两个基线点；选择后点击重新绘制查看效果。")

    def select_baseline_point(self, x_value):
        if not self.baseline_check.isChecked():
            return
        if self.state.calibration_sliced is None:
            self.apply_range_and_baseline(show_popup=False)
        if self.state.calibration_sliced is None:
            return
        x_values = self.state.calibration_sliced.x
        index = nearest_index(x_values, x_value)
        if len(self.selected_baseline) >= 2:
            self.selected_baseline.clear()
            self.baseline_list.clear()
        self.selected_baseline.append(index)
        self.baseline_list.addItem(f"index={index}, x={x_values[index]:.6g}")
        display_y = self._current_sliced_normalized()
        self.plot.plot_curve(x_values, display_y, "基线点选择", markers=tuple(self.selected_baseline))

    def apply_range_and_baseline(self, show_popup=True):
        if self.state.calibration_file is None:
            if show_popup:
                show_error(self, "缺少标定文件", "请先选择标定文件。")
            return
        start_value, end_value = self._validated_range()
        if start_value is None:
            return
        try:
            series = self.reader.read(self.state.calibration_file, start_value, end_value)
            preprocessor = Preprocessor(self.state.smooth_window, self.state.smooth_poly_order)
            normalized = preprocessor.smooth_and_normalize(series.y)
            baseline_indices = tuple(self.selected_baseline) if len(self.selected_baseline) == 2 else None
            corrected = preprocessor.baseline_correct(series.x, normalized, baseline_indices)
        except Exception as error:
            show_error(self, "曲线处理失败", str(error))
            return
        self.state.calibration_sliced = series
        self.state.processed_x = series.x
        self.state.processed_y = corrected
        self.state.start_value = start_value
        self.state.end_value = end_value
        self.state.baseline_indices = baseline_indices
        title = "已基线校准的标定曲线" if baseline_indices else "已应用横坐标范围的标定曲线"
        self.plot.plot_curve(series.x, corrected, title, markers=baseline_indices or ())

    def reset_calibration(self):
        self.selected_baseline.clear()
        self.baseline_list.clear()
        self.baseline_check.setChecked(False)
        self.state.baseline_indices = None
        self.state.calibration_sliced = None
        self.state.processed_x = []
        self.state.processed_y = None
        if self.state.calibration_raw:
            self._configure_range_widgets(self.state.calibration_raw.x)
            self.plot.plot_curve(self.state.calibration_raw.x, self.state.calibration_raw.y, "标定文件原始曲线")

    def validate_and_next(self):
        if self.state.processed_y is None:
            self.apply_range_and_baseline()
        if self.state.processed_y is None:
            return
        self.next_callback()

    def _configure_range_widgets(self, x_values):
        self._syncing_range = True
        minimum = max(0.0, float(min(x_values)))
        maximum = float(max(x_values))
        self.start_spin.setRange(minimum, maximum)
        self.end_spin.setRange(minimum, maximum)
        self.start_spin.setValue(self.state.start_value if self.state.start_value is not None else minimum)
        self.end_spin.setValue(self.state.end_value if self.state.end_value is not None else maximum)
        self.start_slider.setValue(self._value_to_slider(self.start_spin.value()))
        self.end_slider.setValue(self._value_to_slider(self.end_spin.value()))
        self._syncing_range = False

    def _spin_range_changed(self):
        if self._syncing_range or self.state.calibration_raw is None:
            return
        self._syncing_range = True
        self.start_slider.setValue(self._value_to_slider(self.start_spin.value()))
        self.end_slider.setValue(self._value_to_slider(self.end_spin.value()))
        self._syncing_range = False

    def _slider_changed(self, slider_value, spin):
        if self._syncing_range or self.state.calibration_raw is None:
            return
        self._syncing_range = True
        spin.setValue(self._slider_to_value(slider_value))
        self._syncing_range = False

    def _slider_to_value(self, slider_value):
        return self.start_spin.minimum() + (self.start_spin.maximum() - self.start_spin.minimum()) * slider_value / 10000

    def _value_to_slider(self, value):
        width = self.start_spin.maximum() - self.start_spin.minimum()
        return int((value - self.start_spin.minimum()) / width * 10000) if width else 0

    def _validated_range(self):
        start_value = self.start_spin.value()
        end_value = self.end_spin.value()
        if start_value < 0 or end_value < 0:
            show_error(self, "横坐标范围错误", "横坐标必须大于等于 0。")
            return None, None
        if start_value >= end_value:
            show_error(self, "横坐标范围错误", "起点必须小于终点。")
            return None, None
        return start_value, end_value

    def _current_sliced_normalized(self):
        if self.state.calibration_sliced is None:
            return np.asarray([])
        return Preprocessor(self.state.smooth_window, self.state.smooth_poly_order).smooth_and_normalize(self.state.calibration_sliced.y)

class PeakPage(QWidget):
    def __init__(self, state, back_callback, run_callback):
        super().__init__()
        self.state = state
        self.back_callback = back_callback
        self.run_callback = run_callback
        self.selected_peaks: list[int] = []
        self.amorphous_buffer: list[PeakParameter] = []
        self.current_amorphous_index = 0
        self._build()

    def _build(self):
        layout = QVBoxLayout(self)
        main = QHBoxLayout()
        plot_card = make_card()
        plot_layout = QVBoxLayout(plot_card)
        self.plot = PlotPanel()
        self.plot.point_double_clicked.connect(self.select_crystal_peak)
        plot_layout.addWidget(self.plot)
        main.addWidget(plot_card, 2)
        controls = make_card()
        control_layout = QVBoxLayout(controls)
        control_layout.addWidget(section_title("峰信息设置", "双击曲线选择结晶峰；按组输入非晶峰参数并可绘图检查。"))
        self.peak_list = QListWidget()
        reset_peaks_button = danger_button("重置结晶峰选择")
        reset_peaks_button.clicked.connect(self.reset_peaks)
        control_layout.addWidget(QLabel("已选结晶峰"))
        control_layout.addWidget(self.peak_list)
        control_layout.addWidget(reset_peaks_button)

        amorphous_group = QGroupBox("非晶峰参数")
        amorphous_layout = QGridLayout(amorphous_group)
        self.amorphous_count_spin = integer_spin(1, 20, 1, 1)
        self.amorphous_count_spin.valueChanged.connect(self.resize_amorphous_buffer)
        self.amorphous_index_label = QLabel("第 1 / 1 组")
        self.amorphous_height_spin = double_spin(0.000001, 1_000_000, 1.0, 4)
        self.amorphous_center_spin = double_spin(0.0, 1_000_000, 1.0, 4)
        self.amorphous_width_spin = double_spin(0.000001, 1_000_000, 0.2, 4)
        amorphous_layout.addWidget(QLabel("非晶峰个数"), 0, 0)
        amorphous_layout.addWidget(self.amorphous_count_spin, 0, 1)
        amorphous_layout.addWidget(self.amorphous_index_label, 0, 2)
        amorphous_layout.addWidget(QLabel("峰高"), 1, 0)
        amorphous_layout.addWidget(self.amorphous_height_spin, 1, 1, 1, 2)
        amorphous_layout.addWidget(QLabel("峰位置"), 2, 0)
        amorphous_layout.addWidget(self.amorphous_center_spin, 2, 1, 1, 2)
        amorphous_layout.addWidget(QLabel("峰宽"), 3, 0)
        amorphous_layout.addWidget(self.amorphous_width_spin, 3, 1, 1, 2)
        prev_button = QPushButton("上一组")
        next_button = QPushButton("下一组")
        prev_button.clicked.connect(self.previous_amorphous)
        next_button.clicked.connect(self.next_amorphous)
        plot_amorphous_button = QPushButton("绘制非晶峰")
        clear_amorphous_button = danger_button("清空非晶峰绘图")
        plot_amorphous_button.clicked.connect(self.plot_amorphous_peaks)
        clear_amorphous_button.clicked.connect(self.clear_amorphous_plot)
        amorphous_layout.addWidget(prev_button, 4, 0)
        amorphous_layout.addWidget(next_button, 4, 1)
        amorphous_layout.addWidget(plot_amorphous_button, 5, 0, 1, 2)
        amorphous_layout.addWidget(clear_amorphous_button, 5, 2)
        control_layout.addWidget(amorphous_group)
        control_layout.addStretch()
        main.addWidget(controls, 1)
        layout.addLayout(main, stretch=1)
        actions = QHBoxLayout()
        back_button = QPushButton("上一步")
        back_button.clicked.connect(self.back_callback)
        run_button = primary_button("确定并开始处理")
        run_button.clicked.connect(self.validate_and_run)
        actions.addWidget(back_button)
        actions.addStretch()
        actions.addWidget(run_button)
        layout.addLayout(actions)
        self.resize_amorphous_buffer(1)

    def refresh_from_state(self):
        self.selected_peaks = list(self.state.crystal_peak_indices)
        self._refresh_peak_list()
        if self.state.processed_y is not None and self.state.processed_x:
            self.plot.plot_curve(self.state.processed_x, self.state.processed_y, "处理后的标定曲线：双击选择结晶峰", markers=tuple(self.selected_peaks))
        else:
            self.plot.clear("请先完成标定曲线处理")
        if self.state.amorphous_peaks:
            self.amorphous_buffer = list(self.state.amorphous_peaks)
            self.amorphous_count_spin.setValue(len(self.amorphous_buffer))
            self.current_amorphous_index = 0
            self.load_current_amorphous()

    def select_crystal_peak(self, x_value):
        if self.state.processed_y is None or not self.state.processed_x:
            show_error(self, "缺少曲线", "请返回上一步完成标定曲线处理。")
            return
        index = nearest_index(self.state.processed_x, x_value)
        if index not in self.selected_peaks:
            self.selected_peaks.append(index)
        self._refresh_peak_list()
        self.redraw_base_with_current_peaks()

    def reset_peaks(self):
        self.selected_peaks.clear()
        self.state.crystal_peak_indices = ()
        self._refresh_peak_list()
        self.redraw_base_with_current_peaks()

    def resize_amorphous_buffer(self, count):
        self.save_current_amorphous()
        while len(self.amorphous_buffer) < count:
            self.amorphous_buffer.append(PeakParameter(1.0, 1.0, 0.2))
        self.amorphous_buffer = self.amorphous_buffer[:count]
        self.current_amorphous_index = min(self.current_amorphous_index, count - 1)
        self.load_current_amorphous()

    def save_current_amorphous(self):
        if not self.amorphous_buffer:
            return
        self.amorphous_buffer[self.current_amorphous_index] = PeakParameter(
            float(self.amorphous_height_spin.value()),
            float(self.amorphous_center_spin.value()),
            float(self.amorphous_width_spin.value()),
        )

    def load_current_amorphous(self):
        if not self.amorphous_buffer:
            return
        peak = self.amorphous_buffer[self.current_amorphous_index]
        self.amorphous_index_label.setText(f"第 {self.current_amorphous_index + 1} / {len(self.amorphous_buffer)} 组")
        self.amorphous_height_spin.setValue(peak.amplitude)
        self.amorphous_center_spin.setValue(peak.center)
        self.amorphous_width_spin.setValue(peak.width)

    def previous_amorphous(self):
        self.save_current_amorphous()
        self.current_amorphous_index = max(0, self.current_amorphous_index - 1)
        self.load_current_amorphous()

    def next_amorphous(self):
        self.save_current_amorphous()
        self.current_amorphous_index = min(len(self.amorphous_buffer) - 1, self.current_amorphous_index + 1)
        self.load_current_amorphous()

    def plot_amorphous_peaks(self):
        self.save_current_amorphous()
        if self.state.processed_y is None or not self.state.processed_x:
            show_error(self, "缺少曲线", "请返回上一步完成标定曲线处理。")
            return
        self.plot.plot_curve(
            self.state.processed_x,
            self.state.processed_y,
            "结晶峰与非晶峰预览",
            markers=tuple(self.selected_peaks),
            amorphous_peaks=tuple(self.amorphous_buffer),
            show_legend=True,
        )

    def clear_amorphous_plot(self):
        self.redraw_base_with_current_peaks()

    def redraw_base_with_current_peaks(self):
        if self.state.processed_y is not None and self.state.processed_x:
            self.plot.plot_curve(self.state.processed_x, self.state.processed_y, "处理后的标定曲线：双击选择结晶峰", markers=tuple(self.selected_peaks))

    def validate_and_run(self):
        if not self.selected_peaks:
            show_error(self, "缺少结晶峰", "请在绘图区域双击选择至少一个结晶峰。")
            return
        self.save_current_amorphous()
        try:
            self._validate_amorphous_buffer()
        except ValueError as error:
            show_error(self, "非晶峰参数错误", str(error))
            return
        self.state.crystal_peak_indices = tuple(sorted(self.selected_peaks))
        self.state.amorphous_peaks = tuple(self.amorphous_buffer)
        self.run_callback()

    def _validate_amorphous_buffer(self):
        for index, peak in enumerate(self.amorphous_buffer, start=1):
            if peak.amplitude <= 0:
                raise ValueError(f"第 {index} 组峰高必须大于 0。")
            if peak.center < 0:
                raise ValueError(f"第 {index} 组峰位置不能为负。")
            if peak.width <= 0:
                raise ValueError(f"第 {index} 组峰宽必须大于 0。")

    def _refresh_peak_list(self):
        self.peak_list.clear()
        for index in sorted(self.selected_peaks):
            if 0 <= index < len(self.state.processed_x):
                self.peak_list.addItem(f"index={index}, q={self.state.processed_x[index]:.6g}")

class OrientationPage(QWidget):
    def __init__(self, state, back_callback, run_callback):
        super().__init__()
        self.state = state
        self.back_callback = back_callback
        self.run_callback = run_callback
        self.reader = DataReader()
        self.preview_series: DataSeries | None = None
        self._syncing = False
        self._build()

    def _build(self):
        layout = QVBoxLayout(self)
        main = QHBoxLayout()
        controls = make_card()
        control_layout = QVBoxLayout(controls)
        control_layout.addWidget(section_title("取向分析", "设置角度范围，预览过滤和平滑后的曲线，然后批量计算 Pnc。"))
        self.preview_combo = QComboBox()
        load_button = QPushButton("加载预览曲线")
        load_button.clicked.connect(self.load_preview)
        control_layout.addWidget(QLabel("预览文件"))
        control_layout.addWidget(self.preview_combo)
        control_layout.addWidget(load_button)
        range_box = QGroupBox("角度范围")
        range_layout = QGridLayout(range_box)
        self.start_spin = double_spin(0.0, 360.0, 0.0, 3)
        self.end_spin = double_spin(0.0, 360.0, 180.0, 3)
        self.start_slider = QSlider(Qt.Orientation.Horizontal)
        self.end_slider = QSlider(Qt.Orientation.Horizontal)
        for slider in [self.start_slider, self.end_slider]:
            slider.setRange(0, 10000)
        range_layout.addWidget(QLabel("起始角度"), 0, 0)
        range_layout.addWidget(self.start_spin, 0, 1)
        range_layout.addWidget(self.start_slider, 1, 0, 1, 2)
        range_layout.addWidget(QLabel("结束角度"), 2, 0)
        range_layout.addWidget(self.end_spin, 2, 1)
        range_layout.addWidget(self.end_slider, 3, 0, 1, 2)
        control_layout.addWidget(range_box)
        redraw_button = QPushButton("重新绘制预览")
        redraw_button.clicked.connect(self.redraw_preview)
        reset_button = danger_button("重置取向设置")
        reset_button.clicked.connect(self.reset_orientation)
        control_layout.addWidget(redraw_button)
        control_layout.addWidget(reset_button)
        control_layout.addStretch()
        main.addWidget(controls, 1)
        plot_card = make_card()
        plot_layout = QVBoxLayout(plot_card)
        self.plot = PlotPanel()
        plot_layout.addWidget(self.plot)
        main.addWidget(plot_card, 2)
        layout.addLayout(main, stretch=1)
        actions = QHBoxLayout()
        back_button = QPushButton("上一步")
        back_button.clicked.connect(self.back_callback)
        run_button = primary_button("开始取向分析")
        run_button.clicked.connect(self.validate_and_run)
        actions.addWidget(back_button)
        actions.addStretch()
        actions.addWidget(run_button)
        layout.addLayout(actions)
        self.start_spin.valueChanged.connect(self._spin_changed)
        self.end_spin.valueChanged.connect(self._spin_changed)
        self.start_slider.valueChanged.connect(lambda value: self._slider_changed(value, self.start_spin))
        self.end_slider.valueChanged.connect(lambda value: self._slider_changed(value, self.end_spin))
        self.plot.clear("请选择或加载预览文件")

    def refresh_from_state(self):
        self.preview_combo.clear()
        if self.state.input_path is None:
            return
        try:
            files = self.reader.discover(self.state.input_path)
        except Exception as error:
            show_error(self, "输入路径错误", str(error))
            return
        for file_path in files:
            self.preview_combo.addItem(file_path.name, str(file_path))

    def load_preview(self):
        file_path = self._selected_preview_file()
        if file_path is None:
            return
        try:
            self.preview_series = self.reader.read(file_path, None, None)
        except Exception as error:
            show_error(self, "预览文件读取失败", str(error))
            return
        self._configure_range(self.preview_series.x)
        self.redraw_preview()

    def redraw_preview(self):
        if self.preview_series is None:
            self.load_preview()
            return
        start_value, end_value = self._validated_range()
        if start_value is None:
            return
        try:
            series = self.reader.read(self.preview_series.file_path, start_value, end_value)
            preprocessor = Preprocessor(self.state.smooth_window, self.state.smooth_poly_order)
            filtered = preprocessor.filter_outliers(series.y, self.state.outlier_threshold)
            normalized = preprocessor.smooth_and_normalize(filtered.tolist())
        except Exception as error:
            show_error(self, "预览绘制失败", str(error))
            return
        self.plot.plot_curve(series.x, normalized, "取向分析预览曲线", x_label="phi (degree)")

    def reset_orientation(self):
        self.start_spin.setValue(0.0)
        self.end_spin.setValue(180.0)
        self.state.start_value = None
        self.state.end_value = None
        if self.preview_series:
            self.plot.plot_curve(self.preview_series.x, self.preview_series.y, "取向原始预览曲线", x_label="phi (degree)")

    def validate_and_run(self):
        start_value, end_value = self._validated_range()
        if start_value is None:
            return
        self.state.start_value = start_value
        self.state.end_value = end_value
        self.run_callback()

    def _selected_preview_file(self):
        if self.preview_combo.count() == 0:
            show_error(self, "缺少文件", "输入路径中没有可预览的数据文件。")
            return None
        return Path(self.preview_combo.currentData())

    def _configure_range(self, x_values):
        minimum = max(0.0, float(min(x_values)))
        maximum = min(360.0, float(max(x_values)))
        self._syncing = True
        self.start_spin.setRange(minimum, maximum)
        self.end_spin.setRange(minimum, maximum)
        self.start_spin.setValue(self.state.start_value if self.state.start_value is not None else minimum)
        self.end_spin.setValue(self.state.end_value if self.state.end_value is not None else maximum)
        self.start_slider.setValue(self._value_to_slider(self.start_spin.value()))
        self.end_slider.setValue(self._value_to_slider(self.end_spin.value()))
        self._syncing = False

    def _spin_changed(self):
        if self._syncing:
            return
        self._syncing = True
        self.start_slider.setValue(self._value_to_slider(self.start_spin.value()))
        self.end_slider.setValue(self._value_to_slider(self.end_spin.value()))
        self._syncing = False

    def _slider_changed(self, slider_value, spin):
        if self._syncing:
            return
        self._syncing = True
        spin.setValue(self._slider_to_value(slider_value))
        self._syncing = False

    def _slider_to_value(self, slider_value):
        return self.start_spin.minimum() + (self.start_spin.maximum() - self.start_spin.minimum()) * slider_value / 10000

    def _value_to_slider(self, value):
        width = self.start_spin.maximum() - self.start_spin.minimum()
        return int((value - self.start_spin.minimum()) / width * 10000) if width else 0

    def _validated_range(self):
        start_value = self.start_spin.value()
        end_value = self.end_spin.value()
        if start_value < 0 or end_value < 0:
            show_error(self, "角度范围错误", "角度不能为负。")
            return None, None
        if start_value >= end_value:
            show_error(self, "角度范围错误", "起始角度必须小于结束角度。")
            return None, None
        return start_value, end_value


class ResultPage(QWidget):
    def __init__(self, back_callback):
        super().__init__()
        self.back_callback = back_callback
        self._build()

    def _build(self):
        layout = QVBoxLayout(self)
        card = make_card()
        card_layout = QVBoxLayout(card)
        card_layout.addWidget(section_title("运行结果", "处理过程会持续输出状态，失败文件会写入错误报告。"))
        self.progress = QProgressBar()
        self.progress.setRange(0, 0)
        self.log = QTextEdit()
        self.log.setReadOnly(True)
        self.log.setMinimumHeight(430)
        self.back_button = QPushButton("回到首页")
        self.back_button.clicked.connect(self.back_callback)
        card_layout.addWidget(self.progress)
        card_layout.addWidget(self.log, stretch=1)
        card_layout.addWidget(self.back_button, alignment=Qt.AlignmentFlag.AlignRight)
        layout.addWidget(card)

    def reset(self, message):
        self.progress.setRange(0, 0)
        self.log.clear()
        self.append_log(message)

    def append_log(self, message):
        self.log.append(message)

    def finish(self, message):
        self.progress.setRange(0, 1)
        self.progress.setValue(1)
        self.append_log(message)
        show_info(self, "处理完成", message)

    def fail(self, message):
        self.progress.setRange(0, 1)
        self.progress.setValue(0)
        self.append_log(message)
        show_error(self, "处理失败", message)


def main():
    app = QApplication(sys.argv)
    app.setApplicationName("XFit WAXD Batch Analyzer")
    window = XFitWindow()
    window.show()
    sys.exit(app.exec())


XFitApp = XFitWindow


if __name__ == "__main__":
    main()