@echo off
REM 快速测试脚本 - 验证转座子分析环境 (合并版本)

echo ========================================
echo 转座子分析环境测试 (统一版本)
echo ========================================
echo.

REM 测试1: 检查Python
echo [1/4] 测试Python安装...
python --version >nul 2>nul
if %errorlevel% equ 0 (
    echo ✓ Python已安装
    python --version
) else (
    echo ✗ Python未安装或不在PATH中
    echo   请安装Python 3.7或更高版本
)
echo.

REM 测试2: 检查必需文件
echo [2/4] 检查必需文件...
set FILES_OK=1

if exist "transposon_analyzer.py" (
    echo ✓ transposon_analyzer.py (合并分析脚本)
) else (
    echo ✗ transposon_analyzer.py (合并分析脚本)
    set FILES_OK=0
)

if exist "transposon_analysis_simple.nf" (
    echo ✓ transposon_analysis_simple.nf (Nextflow流程)
) else (
    echo ✗ transposon_analysis_simple.nf (Nextflow流程)
    echo   注意: Nextflow流程可选，可以直接使用Python脚本
)

if exist "newplot.txt" (
    echo ✓ newplot.txt (输入数据文件)
) else (
    echo ✗ newplot.txt (输入数据文件)
    echo   警告: 请确保有输入数据文件
)

if exist "nextflow.config" (
    echo ✓ nextflow.config (Nextflow配置文件)
) else (
    echo ✗ nextflow.config (配置文件，可选)
)
echo.

REM 测试3: 验证Python脚本语法
echo [3/4] 验证Python脚本语法...
python -m py_compile transposon_analyzer.py >nul 2>nul
if %errorlevel% equ 0 (
    echo ✓ transposon_analyzer.py 语法正确
) else (
    echo ✗ transposon_analyzer.py 语法错误
    set FILES_OK=0
)
echo.

REM 测试4: 测试Python脚本帮助信息
echo [4/4] 测试脚本功能...
python transposon_analyzer.py --help >nul 2>nul
if %errorlevel% equ 0 (
    echo ✓ Python脚本可以正常执行
    echo.
    echo 脚本帮助信息:
    python transposon_analyzer.py --help
) else (
    echo ✗ Python脚本执行出错
    set FILES_OK=0
)
echo.

REM 检查输入数据格式
if exist "newplot.txt" (
    echo 输入数据检查:
    set /p first_line=<newplot.txt
    echo ✓ 输入文件首行: %first_line%
    
    for /f %%i in ('type "newplot.txt" ^| find /c /v ""') do set LINE_COUNT=%%i
    echo ✓ 数据行数: %LINE_COUNT%
)
echo.

REM 检查可选的Nextflow
echo 可选组件检查:
where nextflow >nul 2>nul
if %errorlevel% equ 0 (
    echo ✓ Nextflow已安装 (可以使用完整流程)
    nextflow -version | findstr "nextflow"
) else (
    echo ○ Nextflow未安装 (可以直接使用Python脚本)
    echo   如需安装: https://www.nextflow.io/docs/latest/getstarted.html
)
echo.

REM 总结
echo ========================================
echo 测试总结
echo ========================================
if %FILES_OK% equ 1 (
    echo ✓ 环境检查通过，可以运行分析流程
    echo.
    echo 推荐运行方式:
    echo 1. 直接使用Python脚本:
    echo    python transposon_analyzer.py --mode all --input newplot.txt --output results
    echo.
    echo 2. 使用批处理脚本:
    echo    run_analysis.bat
    echo.
    if exist "transposon_analysis_simple.nf" (
        echo 3. 使用Nextflow流程 (如果已安装Nextflow):
        echo    nextflow run transposon_analysis_simple.nf
    )
) else (
    echo ✗ 环境检查失败，请检查缺失的文件
    echo.
    echo 问题解决建议:
    echo 1. 确保transposon_analyzer.py文件存在且语法正确
    echo 2. 确保Python 3.7+已正确安装
    echo 3. 准备输入数据文件 (newplot.txt)
)
echo.

pause
