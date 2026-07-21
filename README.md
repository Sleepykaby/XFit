# XFit WAXD Batch Analyzer

XFit 是一个面向 Windows 用户的 WAXD 数据批处理分析软件，支持结晶度、晶粒尺寸和取向因子分析。

## 目录结构

- `xfit/`：正式软件包，包含 UI、算法、批处理、报告导出。
- `scripts/`：Windows 安装辅助脚本。
- `file/`: 文件格式示例，支持'.dat'和'.cake'文件后缀。
- `main.py`：图形界面启动入口。
- `xfit.iss`：xfit 安装脚本（Inno Setup Compile 后使用）。

## 1. 安装前准备

1. 安装 Python 3.10、3.11、3.12 或 3.13。
2. 安装时勾选 `Add python.exe to PATH`。
3. 打开 PowerShell，验证：

```powershell
python --version
```

> 如果科学计算依赖安装失败，建议使用 Python 3.11 或 3.12。

## 2. 安装 XFit

进入项目目录：

```powershell
cd E:\MS_Data\4_Project\XFit
```

创建虚拟环境：

```powershell
python -m venv .venv
```

激活虚拟环境：

```powershell
.\.venv\Scripts\Activate.ps1
```

或使用 conda 创建并激活环境：

```powershell
conda create -n venv python=3.10
conda activate venv
```

如果提示执行策略限制，先运行：

```powershell
Set-ExecutionPolicy -Scope CurrentUser RemoteSigned
```

然后重新激活虚拟环境。

升级 pip：

```powershell
python -m pip install --upgrade pip
```

安装软件和依赖：

```powershell
python -m pip install .
```

安装内容包括：

- `numpy`
- `scipy`
- `matplotlib`
- `PyQt6`
- `xfit-waxd`

验证命令行入口：

```powershell
xfit --help
```

## 3. 启动图形界面

在项目目录中运行：

```powershell
python main.py
```

或安装后运行：

```powershell
xfit-gui
```

## 4. 图形界面使用流程

### 4.1 首页：选择路径和功能

1. 选择输入文件或输入文件夹。
2. 选择输出文件夹，所有结果、图片、日志和错误报告都会保存到这里。
3. 选择分析功能：
   - `结晶度 / 晶粒尺寸`
   - `取向分析`
4. 设置通用参数：
   - 平滑窗口：必须为奇数，且大于多项式阶数。
   - 多项式阶数：用于 Savitzky-Golay 平滑。
   - 取向异常值阈值：仅取向分析使用。
5. 点击 `下一步`。

如果输入路径不存在、输出路径为空、平滑参数不合理，软件会弹窗提示并要求重新设置。

### 4.2 结晶度 / 晶粒尺寸：标定曲线界面

选择 `结晶度 / 晶粒尺寸` 后进入该界面。

1. 点击 `选择`，选择标定文件。
2. 软件读取标定文件并绘制原始曲线。
3. 使用滑动条或数值输入框设置横坐标范围，支持小数。
4. 横坐标范围必须大于 0，且起点必须小于终点；否则会弹窗提示。
5. 点击 `重新绘制 / 检查处理效果`，软件会按当前范围截取曲线并进行平滑归一化。
6. 如需基线校准：
   - 勾选 `开启基线校准模式`。
   - 在绘图区域双击选择两个基线点。
   - 点击 `重新绘制 / 检查处理效果` 查看基线扣除效果。
7. 点击 `重置当前设置` 可清除基线点并恢复原始标定曲线。
8. 点击 `下一步：选择结晶峰`。

该步骤处理后的横坐标范围、基线点和曲线结果会继承到下一步。

### 4.3 结晶度 / 晶粒尺寸：结晶峰和非晶峰界面

1. 界面显示上一步处理后的标定曲线。
2. 在绘图区域双击选择结晶峰位置。
3. 右侧列表会显示已选结晶峰的索引和 q 值。
4. 在非晶峰参数框中输入非晶峰信息，每行一个峰：

```text
1.0, 12.5, 0.4
0.8, 18.0, 0.5
```

每行格式为：

```text
峰高（即高斯峰中心位置的纵坐标）, 峰位（高斯峰中心）, 半峰宽一半（高斯峰标准差）
```

5. 点击 `重置结晶峰` 可清空已选结晶峰。
6. 点击 `上一步` 可返回标定曲线界面重新设置范围或基线。
7. 点击 `确定并开始处理`，软件会批量处理首页选择的输入文件或文件夹，并保存结果。

如果没有选择结晶峰、没有输入非晶峰、非晶峰参数不是数字或小于等于 0，软件会弹窗提示。

### 4.4 取向分析界面

选择 `取向分析` 后进入独立取向分析界面，不会进入结晶度标定流程。

1. 软件会列出输入文件夹或输入文件中的可处理文件。
2. 选择一个文件作为预览文件。
3. 点击 `加载预览曲线`。
4. 使用滑动条或输入框设置角度范围，支持小数。
5. 角度范围不能为负，结束角度必须大于起始角度。
6. 点击 `重新绘制预览`，软件会按原始逻辑进行：
   - 角度范围截取
   - 异常值过滤
   - 平滑归一化
   - 绘制预览曲线
7. 点击 `重置取向设置` 可恢复默认角度设置。
8. 点击 `开始取向分析`，软件会批量处理输入文件并保存取向因子结果。

取向分析的拟合逻辑参考旧版代码：自动寻找最大峰位，使用 Maier-Saupe 函数拟合，并计算 `Pnc`。

## 5. 输出文件说明

输出目录结构示例：

```text
输出目录/
├─ summary.csv
├─ error_report.csv
├─ crystallinity.csv
├─ grain_size.csv
├─ orientation_factor.csv
├─ logs/
│  └─ xfit.log
└─ plots/
   ├─ 文件1.png
   └─ 文件2.png
```

说明：

- `summary.csv`：每个文件的处理状态、指标和输出文件列表。
- `error_report.csv`：失败文件、错误类型、错误消息和完整 traceback。
- `logs/xfit.log`：运行日志。
- `plots/*.png`：每个数据文件的拟合图和残差图。
- `crystallinity.csv`：非晶峰和结晶峰的峰面积结果。
- `grain_size.csv`：峰宽度结果，后续可使用谢乐公式进一步处理得到晶粒尺寸。
- `orientation_factor.csv`：取向因子 `Pnc`。

## 6. 命令行批处理

结晶度 / 晶粒尺寸：

```powershell
xfit "D:\data\waxd" "D:\data\xfit_output" --mode crystallinity --start 5 --end 23 --baseline 10,180 --crystal-peaks 80,120 --amorphous-peaks "1.0,12.5,0.4;0.8,18,0.5"
```

取向分析：

```powershell
xfit "D:\data\orientation" "D:\data\xfit_output" --mode orientation --start 0 --end 180
```

## 7. Windows 安装包

```powershell
Releases --> XFit --> XFitSetup-v2.0.0.exe
```

## 8. 常见问题

### 8.1 无法启动 UI

请确认已经安装 PyQt6：

```powershell
python -m pip install PyQt6
```

### 8.2 依赖安装失败

建议使用 Python 3.11 或 3.12，并重新创建虚拟环境。

### 8.3 处理失败但程序继续运行

这是正常设计。请查看输出目录中的：

- `error_report.csv`
- `logs/xfit.log`

把这两个文件发给开发者即可定位问题。
