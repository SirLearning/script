@echo off
REM 一键运行转座子分析脚本

echo ========================================
echo 转座子分析工具 - 一键运行
echo ========================================
echo.

REM 检查Python是否可用
python --version >nul 2>nul
if %errorlevel% neq 0 (
    echo 错误: Python未安装或不可用
    echo 请先安装Python 3.7或更高版本
    pause
    exit /b 1
)

REM 检查脚本文件
if not exist "transposon_analyzer.py" (
    echo 错误: 找不到分析脚本 transposon_analyzer.py
    pause
    exit /b 1
)

REM 检查输入文件
if not exist "newplot.txt" (
    echo 错误: 找不到输入文件 newplot.txt
    echo 请确保数据文件在当前目录中
    pause
    exit /b 1
)

echo 开始分析转座子数据...
echo 输入文件: newplot.txt
echo 输出目录: results
echo.

REM 运行分析
python transposon_analyzer.py --mode all --input newplot.txt --output results

if %errorlevel% equ 0 (
    echo.
    echo ========================================
    echo 分析成功完成！
    echo ========================================
    echo.
    echo 结果文件位置:
    echo - 详细报告: results\transposon_analysis_report.html
    echo - 文本报告: results\transposon_analysis_report.txt
    echo - 组成汇总: results\composition_summary.txt
    echo.
    echo 要查看HTML报告，请在浏览器中打开:
    echo results\transposon_analysis_report.html
    echo.
    
    REM 询问是否打开报告
    set /p open_report="要现在打开HTML报告吗？(y/N): "
    if /i "%open_report%"=="y" (
        start "" "results\transposon_analysis_report.html"
    )
) else (
    echo.
    echo ========================================
    echo 分析失败！
    echo ========================================
    echo 请检查错误信息并重试
)

echo.
pause
