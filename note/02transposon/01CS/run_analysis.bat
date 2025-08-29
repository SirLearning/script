@echo off
REM Windows批处理脚本用于运行转座子分析流程 (合并版本)

echo ========================================
echo 转座子分析流程 - Windows运行脚本 (统一版本)
echo ========================================
echo.

REM 检查Python是否已安装
python --version >nul 2>nul
if %errorlevel% neq 0 (
    echo 错误: 未找到Python。请先安装Python 3.7+。
    pause
    exit /b 1
)

REM 显示当前目录和文件
echo 当前工作目录: %cd%
echo.
echo 检查必需文件...

if not exist "transposon_analyzer.py" (
    echo 错误: 未找到合并脚本文件 transposon_analyzer.py
    pause
    exit /b 1
)

if not exist "newplot.txt" (
    echo 警告: 未找到输入文件 newplot.txt
    echo 请确保数据文件在当前目录中
)

echo 所有必需文件检查完成。
echo.

REM 设置默认参数
set INPUT_DATA=newplot.txt
set OUTPUT_DIR=results
set MODE=all
set GFF_FILE=
set LIB_FILE=

REM 解析命令行参数
:parse_args
if "%1"=="--input" (
    set INPUT_DATA=%2
    shift
    shift
    goto parse_args
)
if "%1"=="--output" (
    set OUTPUT_DIR=%2
    shift
    shift
    goto parse_args
)
if "%1"=="--mode" (
    set MODE=%2
    shift
    shift
    goto parse_args
)
if "%1"=="--gff" (
    set GFF_FILE=%2
    shift
    shift
    goto parse_args
)
if "%1"=="--lib" (
    set LIB_FILE=%2
    shift
    shift
    goto parse_args
)
if "%1"=="--help" (
    goto show_help
)
if "%1"=="/?" (
    goto show_help
)
if not "%1"=="" (
    shift
    goto parse_args
)

echo 运行参数:
echo   输入文件: %INPUT_DATA%
echo   输出目录: %OUTPUT_DIR%
echo   分析模式: %MODE%
if not "%GFF_FILE%"=="" echo   GFF文件: %GFF_FILE%
if not "%LIB_FILE%"=="" echo   库文件: %LIB_FILE%
echo.

REM 提供运行选项
echo 选择运行方式:
echo 1) 直接运行Python脚本 (推荐)
echo 2) 运行Nextflow流程 (需要安装Nextflow)
echo 3) 退出
set /p choice="请选择 (1-3): "

if "%choice%"=="1" goto run_python
if "%choice%"=="2" goto run_nextflow
if "%choice%"=="3" goto exit_script

:run_python
echo.
echo 开始运行Python分析脚本...
echo ========================================

REM 构建Python命令
set PYTHON_CMD=python transposon_analyzer.py --mode %MODE% --input "%INPUT_DATA%" --output "%OUTPUT_DIR%"

if not "%GFF_FILE%"=="" (
    set PYTHON_CMD=%PYTHON_CMD% --gff "%GFF_FILE%"
)

if not "%LIB_FILE%"=="" (
    set PYTHON_CMD=%PYTHON_CMD% --lib "%LIB_FILE%"
)

echo 执行命令: %PYTHON_CMD%
echo.

%PYTHON_CMD%

if %errorlevel% equ 0 (
    echo.
    echo ========================================
    echo Python分析完成！
    echo 结果保存在: %OUTPUT_DIR%
    echo ========================================
    echo.
    echo 生成的文件:
    if exist "%OUTPUT_DIR%" (
        dir /b "%OUTPUT_DIR%"
    )
) else (
    echo.
    echo ========================================
    echo Python分析失败！错误代码: %errorlevel%
    echo ========================================
)

goto end_script

:run_nextflow
echo.
echo 检查Nextflow安装...
where nextflow >nul 2>nul
if %errorlevel% neq 0 (
    echo 错误: 未找到Nextflow。
    echo 请先安装Nextflow或选择直接运行Python脚本。
    pause
    exit /b 1
)

echo 开始运行Nextflow流程...
echo ========================================

nextflow run transposon_analysis_simple.nf ^
    --input_data "%INPUT_DATA%" ^
    --output_dir "%OUTPUT_DIR%"

if %errorlevel% equ 0 (
    echo.
    echo ========================================
    echo Nextflow分析完成！
    echo 结果保存在: %OUTPUT_DIR%
    echo ========================================
) else (
    echo.
    echo ========================================
    echo Nextflow分析失败！错误代码: %errorlevel%
    echo ========================================
)

goto end_script

:show_help
echo.
echo 用法: run_analysis.bat [选项]
echo.
echo 选项:
echo   --input FILE     指定输入数据文件 (默认: newplot.txt)
echo   --output DIR     指定输出目录 (默认: results)
echo   --mode MODE      指定分析模式 (默认: all)
echo                    可选: composition, gff, library, qc, all
echo   --gff FILE       指定GFF3文件 (可选)
echo   --lib FILE       指定库文件 (可选)
echo   --help, /?       显示此帮助信息
echo.
echo 分析模式说明:
echo   composition      仅进行组成分析
echo   gff              仅进行GFF3分析
echo   library          仅进行库文件分析
echo   qc               仅进行质量控制
echo   all              进行所有分析 (推荐)
echo.
echo 示例:
echo   run_analysis.bat
echo   run_analysis.bat --input mydata.txt --output my_results
echo   run_analysis.bat --mode composition --input newplot.txt
echo   run_analysis.bat --gff chr1A.gff3 --lib library.lib
echo.
pause
exit /b 0

:exit_script
echo 操作已取消。
goto end_script

:end_script
echo.
pause
exit /b %errorlevel%
