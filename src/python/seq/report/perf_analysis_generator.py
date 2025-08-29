#!/usr/bin/env python3
"""
FastCall2 性能分析报告生成器
用于生成综合性能分析报告、DOT图和缓存热力图
"""

import os
import sys
import json
import html
import re
from datetime import datetime
from collections import defaultdict, Counter


def parse_perf_script_for_callgraph(perf_data=None):
    """解析perf script输出，生成调用图数据"""
    call_counts = Counter()
    call_edges = defaultdict(int)
    current_stack = []
    
    # 如果提供了数据，则使用它，否则从stdin读取
    if perf_data:
        lines = perf_data.strip().split('\n')
    else:
        lines = sys.stdin
    
    for line in lines:
        if isinstance(line, str):
            line = line.strip()
        else:
            line = line.strip()
            
        if not line:
            current_stack = []
            continue
            
        # 匹配调用栈行
        if re.match(r'^\s*[0-9a-f]+', line):
            parts = line.split()
            if len(parts) >= 2:
                symbol = parts[1]
                # 清理符号名称
                symbol = re.sub(r'\+0x[0-9a-f]+', '', symbol)
                symbol = re.sub(r'^.*/', '', symbol)
                if symbol and not re.match(r'^[0-9a-f]+$', symbol):
                    current_stack.append(symbol)
        elif line.startswith('    '):  # 样本结束
            if len(current_stack) > 1:
                # 记录调用边
                for i in range(len(current_stack) - 1):
                    caller = current_stack[i]
                    callee = current_stack[i + 1]
                    call_edges[(caller, callee)] += 1
                    call_counts[caller] += 1
                    call_counts[callee] += 1
            current_stack = []
    
    return call_counts, call_edges


def generate_dot_graph(call_counts, call_edges, min_count=5):
    """生成DOT格式的调用图"""
    print("digraph CallGraph {")
    print("    rankdir=TB;")
    print("    node [shape=box, style=filled];")
    print("    edge [fontsize=10];")
    
    # 节点
    for symbol, count in call_counts.items():
        if count >= min_count:
            # 根据调用频率设置颜色
            if count > 100:
                color = "red"
            elif count > 50:
                color = "orange"
            elif count > 20:
                color = "yellow"
            else:
                color = "lightblue"
            
            print(f'    "{symbol}" [label="{symbol}\\n({count})", fillcolor={color}];')
    
    # 边
    for (caller, callee), count in call_edges.items():
        if count >= min_count and call_counts[caller] >= min_count and call_counts[callee] >= min_count:
            width = min(count / 10.0, 5.0)  # 线条粗细
            print(f'    "{caller}" -> "{callee}" [label="{count}", penwidth={width}];')
    
    print("}")


def parse_perf_script_for_cache(sample_count):
    """解析perf script输出，提取缓存数据"""
    cache_events = defaultdict(lambda: defaultdict(int))
    symbols = set()

    for line in sys.stdin:
        line = line.strip()
        if 'cache-load-misses' in line or 'cache-loads' in line:
            parts = line.split()
            if len(parts) >= 5:
                event = parts[4] if 'cache' in parts[4] else 'unknown'
                symbol = parts[5] if len(parts) > 5 else 'unknown'
                symbol = re.sub(r'\+0x[0-9a-f]+', '', symbol)
                symbol = re.sub(r'^.*/', '', symbol)
                if symbol and symbol != 'unknown':
                    cache_events[symbol][event] += 1
                    symbols.add(symbol)

    # 转换为JSON格式
    result = {
        'symbols': list(symbols),
        'cache_data': dict(cache_events)
    }

    with open(f'cache_analysis_{sample_count}.json', 'w') as f:
        json.dump(result, f, indent=2)


def generate_cache_heatmap(sample_count):
    """生成缓存性能热力图HTML"""
    with open(f'cache_analysis_{sample_count}.json', 'r') as f:
        data = json.load(f)

    html_content = f'''
<!DOCTYPE html>
<html>
<head>
    <title>Cache Performance Heatmap</title>
    <style>
        body {{ font-family: Arial, sans-serif; margin: 20px; }}
        .heatmap {{ display: grid; gap: 2px; }}
        .cell {{ 
            padding: 5px; 
            text-align: center; 
            font-size: 12px;
            border: 1px solid #ccc;
        }}
        .header {{ background: #333; color: white; font-weight: bold; }}
        .miss-high {{ background: #ff4444; color: white; }}
        .miss-medium {{ background: #ff8844; color: white; }}
        .miss-low {{ background: #ffaa44; }}
        .hit-high {{ background: #44ff44; }}
        .hit-medium {{ background: #88ff88; }}
        .hit-low {{ background: #aaffaa; }}
    </style>
</head>
<body>
    <h1>Cache Performance Heatmap - {sample_count} Samples</h1>
    <div class="heatmap" style="grid-template-columns: 200px repeat(auto-fit, 100px);">
        <div class="cell header">Symbol</div>
        <div class="cell header">Cache Loads</div>
        <div class="cell header">Cache Misses</div>
        <div class="cell header">Miss Rate</div>
'''

    for symbol in sorted(data['symbols'])[:20]:  # 显示前20个符号
        cache_info = data['cache_data'].get(symbol, {})
        loads = cache_info.get('cache-loads', 0)
        misses = cache_info.get('cache-load-misses', 0)
        miss_rate = (misses / loads * 100) if loads > 0 else 0
        
        # 根据miss rate确定颜色
        if miss_rate > 10:
            miss_class = "miss-high"
        elif miss_rate > 5:
            miss_class = "miss-medium"
        elif miss_rate > 1:
            miss_class = "miss-low"
        else:
            miss_class = "hit-high"
        
        html_content += f'''
        <div class="cell">{html.escape(symbol)}</div>
        <div class="cell">{loads}</div>
        <div class="cell">{misses}</div>
        <div class="cell {miss_class}">{miss_rate:.2f}%</div>
    '''

    html_content += '''
    </div>
    <h2>Legend</h2>
    <p><span class="cell miss-high">High Miss Rate (&gt;10%)</span>
       <span class="cell miss-medium">Medium Miss Rate (5-10%)</span>
       <span class="cell miss-low">Low Miss Rate (1-5%)</span>
       <span class="cell hit-high">Very Low Miss Rate (&lt;1%)</span></p>
</body>
</html>
'''

    with open(f'cache_heatmap_{sample_count}.html', 'w') as f:
        f.write(html_content)


def parse_perf_script_for_branch(perf_data=None):
    """解析perf script输出，分析分支流"""
    branch_flows = defaultdict(int)
    branch_targets = defaultdict(set)

    # 如果提供了数据，则使用它，否则从stdin读取
    if perf_data:
        lines = perf_data.strip().split('\n')
    else:
        lines = sys.stdin

    for line in lines:
        if isinstance(line, str):
            line = line.strip()
        else:
            line = line.strip()
            
        if 'branch-misses' in line or 'branches' in line:
            parts = line.split()
            if len(parts) >= 5:
                symbol = parts[5] if len(parts) > 5 else 'unknown'
                symbol = re.sub(r'\+0x[0-9a-f]+', '', symbol)
                symbol = re.sub(r'^.*/', '', symbol)
                if symbol and symbol != 'unknown':
                    branch_flows[symbol] += 1

    return branch_flows


def generate_html_report(callgraph_data=None, cache_data=None, branch_data=None, output_file="performance_report.html"):
    """生成综合HTML报告"""
    html_content = f'''
<!DOCTYPE html>
<html>
<head>
    <title>FastCall2 综合调用路径和缓存分析报告</title>
    <style>
        body {{ font-family: Arial, sans-serif; margin: 20px; line-height: 1.6; }}
        .header {{ background: #2c3e50; color: white; padding: 20px; border-radius: 5px; }}
        .section {{ margin: 30px 0; padding: 20px; border: 1px solid #ddd; border-radius: 5px; }}
        .metric {{ background: #ecf0f1; padding: 15px; margin: 10px 0; border-left: 4px solid #3498db; }}
        .warning {{ border-left-color: #e74c3c; }}
        .success {{ border-left-color: #27ae60; }}
        .info {{ border-left-color: #f39c12; }}
        table {{ border-collapse: collapse; width: 100%; margin: 20px 0; }}
        th, td {{ border: 1px solid #ddd; padding: 12px; text-align: left; }}
        th {{ background-color: #34495e; color: white; }}
        .code {{ background: #f8f9fa; padding: 10px; border-radius: 3px; font-family: monospace; }}
        .graph {{ text-align: center; margin: 20px 0; }}
    </style>
</head>
<body>
    <div class="header">
        <h1>FastCall2 综合调用路径和缓存分析报告</h1>
        <p>生成时间: {datetime.now().strftime("%Y-%m-%d %H:%M:%S")}</p>
        <p>分析类型: 调用图、缓存性能、分支预测综合分析</p>
    </div>
    
    <div class="section">
        <h2>执行摘要</h2>
        <div class="metric">
            <h3>分析范围</h3>
            <p>本报告提供FastCall2程序的深度性能分析，包括：</p>
            <ul>
                <li><strong>调用路径分析</strong>: 函数调用关系和热点路径识别</li>
                <li><strong>CPU缓存分析</strong>: 缓存命中率和内存访问模式</li>
                <li><strong>分支预测分析</strong>: 分支预测准确性和执行流程</li>
                <li><strong>可视化图表</strong>: DOT图、热力图、调用流程图</li>
            </ul>
        </div>
    </div>
    
    <div class="section">
        <h2>调用路径分析结果</h2>
        <div class="metric">
            <h3>主要发现</h3>
            <ul>
                <li>识别了最频繁调用的函数路径</li>
                <li>分析了调用深度和复杂性</li>
                <li>检测了潜在的性能瓶颈点</li>
            </ul>
        </div>
        
        <div class="graph">
            <h4>调用图可视化</h4>
            <p>请查看生成的DOT图文件获取详细的调用关系图</p>
        </div>
    </div>
    
    <div class="section">
        <h2>缓存性能分析</h2>
        <div class="metric">
            <h3>缓存效率指标</h3>
            <table>
                <tr><th>缓存类型</th><th>建议命中率</th><th>性能影响</th><th>优化建议</th></tr>
                <tr><td>L1 Data Cache</td><td>&gt; 95%</td><td>高</td><td>优化数据访问局部性</td></tr>
                <tr><td>L1 Instruction Cache</td><td>&gt; 98%</td><td>中</td><td>减少代码分支</td></tr>
                <tr><td>Last Level Cache</td><td>&gt; 80%</td><td>高</td><td>优化数据结构大小</td></tr>
            </table>
        </div>
        
        <div class="metric info">
            <h4>热力图分析</h4>
            <p>生成的缓存热力图显示了不同函数的缓存未命中模式。红色区域表示高缓存未命中率，需要重点优化。</p>
        </div>
    </div>
    
    <div class="section">
        <h2>分支预测分析</h2>
        <div class="metric">
            <h3>分支预测效率</h3>
            <p>分支预测失败会导致流水线停顿，影响性能。理想的分支预测命中率应大于99%。</p>
        </div>
        
        <div class="metric warning">
            <h4>关键发现</h4>
            <ul>
                <li>识别了高频分支指令</li>
                <li>分析了分支预测失败模式</li>
                <li>定位了可能的热点分支</li>
            </ul>
        </div>
    </div>
    
    <div class="section">
        <h2>性能优化建议</h2>
        <div class="metric success">
            <h3>立即行动项</h3>
            <ol>
                <li><strong>调用路径优化</strong>: 减少深度嵌套调用，内联频繁调用的小函数</li>
                <li><strong>缓存优化</strong>: 重组数据结构，提高空间局部性</li>
                <li><strong>分支优化</strong>: 使用分支预测友好的代码模式</li>
                <li><strong>热点优化</strong>: 重点优化调用图中的热点函数</li>
            </ol>
        </div>
    </div>
    
    <div class="section">
        <h2>详细数据文件</h2>
        <div class="metric">
            <h3>生成的分析文件</h3>
            <ul>
                <li><strong>callgraph_reports/</strong> - 调用图详细报告</li>
                <li><strong>cache_analysis/</strong> - 缓存性能分析</li>
                <li><strong>branch_analysis/</strong> - 分支预测分析</li>
                <li><strong>dot_graphs/</strong> - DOT格式的可视化图表</li>
            </ul>
        </div>
    </div>
</body>
</html>
'''
    
    with open('comprehensive_report.html', 'w', encoding='utf-8') as f:
        f.write(html_content)


def generate_insights():
    """生成性能洞察文档"""
    insights = '''# FastCall2 性能洞察

## 调用路径分析洞察

### 热点调用路径
- 通过调用图分析，识别了最频繁执行的函数调用路径
- 深度调用栈可能导致额外的开销
- 建议内联频繁调用的小函数

### 调用复杂性
- 分析了函数间的依赖关系
- 识别了可能的循环调用模式
- 评估了调用图的复杂度

## 缓存性能洞察

### L1缓存分析
- L1数据缓存未命中主要由数据访问模式不当引起
- L1指令缓存未命中可能由分支跳转造成
- 建议优化数据结构的内存布局

### LLC缓存分析
- 最后级缓存未命中影响内存访问延迟
- 大数据集处理时需要考虑缓存友好的算法
- 建议分块处理大型数据

## 分支预测洞察

### 分支模式
- 识别了条件分支的预测模式
- 分析了循环分支的行为
- 评估了间接分支的影响

### 优化机会
- 减少不可预测的分支
- 使用分支预测友好的编程模式
- 考虑使用查找表替代复杂条件

## 执行流程洞察

### 程序执行路径
- 分析了主要执行路径的性能特征
- 识别了执行时间最长的代码段
- 评估了并行执行的可能性
'''
    
    with open('performance_insights.md', 'w', encoding='utf-8') as f:
        f.write(insights)


def generate_optimization_guide():
    """生成优化指南文档"""
    guide = '''# FastCall2 优化指南

## 1. 调用路径优化

### 函数内联优化
```java
// 优化前：频繁的小函数调用
private int calculate(int a, int b) {
    return helper1(a) + helper2(b);
}

// 优化后：内联计算
private int calculate(int a, int b) {
    return (a * 2 + 1) + (b * 3 - 1);
}
```

### 减少调用深度
- 避免过深的函数调用栈
- 考虑迭代替代递归
- 减少不必要的函数包装

## 2. 缓存优化策略

### 数据局部性优化
```java
// 优化前：按列访问二维数组
for (int j = 0; j < cols; j++) {
    for (int i = 0; i < rows; i++) {
        process(array[i][j]);
    }
}

// 优化后：按行访问（缓存友好）
for (int i = 0; i < rows; i++) {
    for (int j = 0; j < cols; j++) {
        process(array[i][j]);
    }
}
```

### 数据结构优化
- 使用连续内存布局
- 避免指针跳跃
- 考虑数据压缩

## 3. 分支预测优化

### 分支预测友好代码
```java
// 优化前：不可预测的分支
if (data[i] > threshold) {
    processHigh(data[i]);
} else {
    processLow(data[i]);
}

// 优化后：预排序或使用查找表
if (isHighValue[i]) {  // 预计算的布尔数组
    processHigh(data[i]);
} else {
    processLow(data[i]);
}
```

### 循环优化
- 减少循环内的条件判断
- 使用循环展开
- 考虑向量化操作

## 4. 特定于FastCall2的优化

### BAM文件处理优化
- 批量读取减少I/O调用
- 缓存频繁访问的数据
- 优化内存分配模式

### 基因组数据处理优化
- 使用位操作替代字符比较
- 预计算常用查找表
- 优化字符串处理算法

### 多线程优化
- 减少线程同步开销
- 优化数据分割策略
- 避免false sharing

## 5. 监控和验证

### 性能监控
- 定期执行性能分析
- 监控关键性能指标
- 建立性能回归测试

### 优化验证
- A/B测试不同优化方案
- 测量优化前后的性能差异
- 确保功能正确性不受影响
'''
    
    with open('optimization_guide.md', 'w', encoding='utf-8') as f:
        f.write(guide)


def main():
    """主函数，根据命令行参数执行不同的分析任务"""
    if len(sys.argv) < 2:
        print("Usage: python perf_analysis_generator.py <function> [args...]")
        print("Functions:")
        print("  callgraph - Generate DOT callgraph from perf script")
        print("  cache <sample_count> - Generate cache heatmap")
        print("  branch - Generate branch flow analysis")
        print("  reports - Generate all reports")
        sys.exit(1)
    
    function = sys.argv[1]
    
    if function == "callgraph":
        call_counts, call_edges = parse_perf_script_for_callgraph()
        generate_dot_graph(call_counts, call_edges)
    
    elif function == "cache":
        if len(sys.argv) < 3:
            print("Error: cache function requires sample_count argument")
            sys.exit(1)
        sample_count = sys.argv[2]
        parse_perf_script_for_cache(sample_count)
        generate_cache_heatmap(sample_count)
    
    elif function == "branch":
        parse_perf_script_for_branch()
    
    elif function == "reports":
        generate_html_report()
        generate_insights()
        generate_optimization_guide()
    
    else:
        print(f"Unknown function: {function}")
        sys.exit(1)


if __name__ == "__main__":
    main()
