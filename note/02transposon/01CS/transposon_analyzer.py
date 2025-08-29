#!/usr/bin/env python3
"""
Transposon Analysis Tool
合并的转座子分析工具，包含多种分析功能

功能:
1. 转座子组成分析 (newplot.txt格式)
2. GFF3文件类别统计
3. 库文件统计分析
4. 数据质量控制

使用方法:
python transposon_analyzer.py --mode composition --input newplot.txt --output results/
python transposon_analyzer.py --mode gff --input chr1A.gff3 --output results/
python transposon_analyzer.py --mode library --input library.lib --output results/
python transposon_analyzer.py --mode all --input newplot.txt --gff chr1A.gff3 --lib library.lib --output results/
"""

import argparse
import os
import sys
import re
from datetime import datetime
from pathlib import Path
import json


class TransposonAnalyzer:
    def __init__(self, output_dir="results"):
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)
        self.results = {}
        
    def log_info(self, message):
        """记录信息到日志"""
        timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        print(f"[{timestamp}] {message}")
        
    def analyze_composition(self, input_file):
        """分析转座子组成 - 基于newplot.txt格式"""
        self.log_info(f"开始分析转座子组成: {input_file}")
        
        if not os.path.exists(input_file):
            raise FileNotFoundError(f"输入文件不存在: {input_file}")
            
        sums = {}
        details = []
        
        try:
            with open(input_file, 'r') as file:
                for line_num, line in enumerate(file, 1):
                    line = line.strip()
                    if not line:
                        continue
                        
                    # 分割行内容
                    parts = line.split()
                    if len(parts) < 2:
                        continue
                        
                    last_part = ''
                    max_percentage = 0
                    max_word = ''
                    count_number = 0
                    
                    # 解析每个部分
                    for part in parts:
                        match = re.match(r'(\w+)_\w+', last_part)
                        if match:
                            word = match.group(1)  # 修正: 使用group(1)而不是groups()
                            if ':' in part:
                                partc = part.split(":")
                                try:
                                    percentage = float(partc[0])
                                    if percentage > max_percentage:
                                        max_percentage = percentage
                                        max_word = word
                                except ValueError:
                                    continue
                        last_part = part
                    
                    # 提取计数数字
                    try:
                        count_number = int(parts[-1])
                    except ValueError:
                        continue
                    
                    # 记录详细信息
                    detail = {
                        'line': line_num,
                        'parts': parts,
                        'max_percentage': max_percentage,
                        'max_word': max_word,
                        'count_number': count_number
                    }
                    details.append(detail)
                    
                    # 累计统计 - 使用前三个字母作为键
                    if max_word:
                        key = max_word[:3] if len(max_word) >= 3 else max_word
                        sums[key] = sums.get(key, 0) + count_number
                        
        except Exception as e:
            raise RuntimeError(f"分析组成文件时出错: {e}")
        
        # 保存结果
        composition_results = {
            'summary': sums,
            'details': details,
            'total_lines': len(details),
            'analysis_time': datetime.now().isoformat()
        }
        
        # 写入详细统计文件
        stats_file = self.output_dir / "composition_stats.txt"
        with open(stats_file, 'w', encoding='utf-8') as f:
            f.write("转座子组成分析详细结果\n")
            f.write("=" * 50 + "\n")
            f.write(f"分析时间: {composition_results['analysis_time']}\n")
            f.write(f"总行数: {composition_results['total_lines']}\n\n")
            
            for detail in details:
                f.write(f"行 {detail['line']}: {detail['parts']}\n")
                f.write(f"  最高百分比: {detail['max_percentage']}\n")
                f.write(f"  最高词: {detail['max_word']}\n")
                f.write(f"  计数: {detail['count_number']}\n\n")
        
        # 写入汇总文件
        summary_file = self.output_dir / "composition_summary.txt"
        with open(summary_file, 'w', encoding='utf-8') as f:
            f.write("转座子组成分析汇总\n")
            f.write("=" * 30 + "\n")
            f.write(f"分析时间: {composition_results['analysis_time']}\n")
            f.write(f"输入文件: {input_file}\n")
            f.write(f"总行数: {composition_results['total_lines']}\n\n")
            f.write("类别统计:\n")
            for key, value in sorted(sums.items()):
                f.write(f"{key}: {value}\n")
        
        self.results['composition'] = composition_results
        self.log_info(f"组成分析完成，结果保存到: {stats_file}, {summary_file}")
        return composition_results
    
    def analyze_gff(self, gff_file):
        """分析GFF3文件中的转座子类别"""
        self.log_info(f"开始分析GFF3文件: {gff_file}")
        
        if not os.path.exists(gff_file):
            self.log_info(f"GFF文件不存在，跳过分析: {gff_file}")
            return None
            
        category_counts = {}
        total_lines = 0
        processed_lines = 0
        
        try:
            with open(gff_file, 'r') as file:
                for line in file:
                    total_lines += 1
                    
                    # 跳过注释行
                    if line.startswith('#'):
                        continue
                    
                    # 分割字段
                    fields = line.strip().split('\t')
                    if len(fields) < 9:
                        continue
                    
                    # 提取属性字段
                    attributes = fields[8]
                    attribute_pairs = []
                    
                    for pair in attributes.split(';'):
                        if '=' in pair:
                            key, value = pair.split('=', 1)
                            attribute_pairs.append((key.strip(), value.strip()))
                    
                    # 查找compo属性
                    category = None
                    for key, value in attribute_pairs:
                        if key == 'compo':
                            category = value
                            break
                    
                    # 如果没找到compo，使用最后一个键作为类别
                    if not category and attribute_pairs:
                        category = attribute_pairs[-1][0]
                    
                    if category:
                        category_counts[category] = category_counts.get(category, 0) + 1
                        processed_lines += 1
                        
        except Exception as e:
            raise RuntimeError(f"分析GFF文件时出错: {e}")
        
        # 保存结果
        gff_results = {
            'category_counts': category_counts,
            'total_lines': total_lines,
            'processed_lines': processed_lines,
            'analysis_time': datetime.now().isoformat()
        }
        
        # 写入类别计数文件
        counts_file = self.output_dir / "gff_category_counts.txt"
        with open(counts_file, 'w', encoding='utf-8') as f:
            for category, count in category_counts.items():
                f.write(f'{category}: {count}\n')
        
        # 写入分析报告
        report_file = self.output_dir / "gff_analysis_report.txt"
        with open(report_file, 'w', encoding='utf-8') as f:
            f.write("GFF3转座子类别分析报告\n")
            f.write("=" * 40 + "\n")
            f.write(f"分析时间: {gff_results['analysis_time']}\n")
            f.write(f"输入文件: {gff_file}\n")
            f.write(f"总行数: {gff_results['total_lines']}\n")
            f.write(f"处理行数: {gff_results['processed_lines']}\n\n")
            f.write("类别统计:\n")
            for category, count in sorted(category_counts.items()):
                f.write(f"{category}: {count}\n")
        
        self.results['gff'] = gff_results
        self.log_info(f"GFF分析完成，结果保存到: {counts_file}, {report_file}")
        return gff_results
    
    def analyze_library(self, lib_file):
        """分析库文件统计"""
        self.log_info(f"开始分析库文件: {lib_file}")
        
        if not os.path.exists(lib_file):
            self.log_info(f"库文件不存在，跳过分析: {lib_file}")
            return None
            
        lib_data = {}
        total_entries = 0
        
        try:
            with open(lib_file, 'r') as file:
                for line in file:
                    line = line.strip()
                    if not line:
                        continue
                        
                    elements = line.split('\t')
                    if len(elements) < 3:
                        continue
                    
                    chrom = elements[0]
                    try:
                        start = int(elements[1])
                        end = int(elements[2])
                        length = end - start + 1
                        
                        if chrom not in lib_data:
                            lib_data[chrom] = []
                        
                        lib_data[chrom].append({
                            'start': start,
                            'end': end,
                            'length': length
                        })
                        total_entries += 1
                        
                    except ValueError:
                        continue
                        
        except Exception as e:
            raise RuntimeError(f"分析库文件时出错: {e}")
        
        # 计算统计信息
        chrom_stats = {}
        for chrom, entries in lib_data.items():
            lengths = [entry['length'] for entry in entries]
            chrom_stats[chrom] = {
                'count': len(entries),
                'total_length': sum(lengths),
                'avg_length': sum(lengths) / len(lengths) if lengths else 0,
                'min_length': min(lengths) if lengths else 0,
                'max_length': max(lengths) if lengths else 0
            }
        
        lib_results = {
            'library_data': lib_data,
            'chromosome_stats': chrom_stats,
            'total_entries': total_entries,
            'analysis_time': datetime.now().isoformat()
        }
        
        # 写入库统计文件
        stats_file = self.output_dir / f"lib_stats_{Path(lib_file).stem}.txt"
        with open(stats_file, 'w', encoding='utf-8') as f:
            f.write("库文件统计分析\n")
            f.write("=" * 30 + "\n")
            f.write(f"分析时间: {lib_results['analysis_time']}\n")
            f.write(f"输入文件: {lib_file}\n")
            f.write(f"总条目数: {lib_results['total_entries']}\n\n")
            f.write("染色体统计:\n")
            for chrom, stats in chrom_stats.items():
                f.write(f"{chrom}:\n")
                f.write(f"  条目数: {stats['count']}\n")
                f.write(f"  总长度: {stats['total_length']}\n")
                f.write(f"  平均长度: {stats['avg_length']:.2f}\n")
                f.write(f"  最小长度: {stats['min_length']}\n")
                f.write(f"  最大长度: {stats['max_length']}\n\n")
        
        self.results['library'] = lib_results
        self.log_info(f"库分析完成，结果保存到: {stats_file}")
        return lib_results
    
    def quality_control(self, input_file):
        """数据质量控制分析"""
        self.log_info(f"开始质量控制分析: {input_file}")
        
        if not os.path.exists(input_file):
            self.log_info(f"文件不存在，跳过QC: {input_file}")
            return None
        
        # 基本统计
        file_size = os.path.getsize(input_file)
        total_lines = 0
        empty_lines = 0
        malformed_lines = 0
        sample_lines = []
        
        try:
            with open(input_file, 'r') as file:
                for line_num, line in enumerate(file, 1):
                    total_lines += 1
                    
                    # 收集前5行作为样本
                    if line_num <= 5:
                        sample_lines.append(line.strip())
                    
                    # 检查空行
                    if not line.strip():
                        empty_lines += 1
                        continue
                    
                    # 检查格式错误的行
                    if not re.match(r'^[A-Za-z]', line.strip()):
                        malformed_lines += 1
                        
        except Exception as e:
            raise RuntimeError(f"质量控制分析时出错: {e}")
        
        qc_results = {
            'file_size': file_size,
            'total_lines': total_lines,
            'empty_lines': empty_lines,
            'malformed_lines': malformed_lines,
            'sample_lines': sample_lines,
            'analysis_time': datetime.now().isoformat()
        }
        
        # 写入QC报告
        qc_report_file = self.output_dir / "qc_report.txt"
        with open(qc_report_file, 'w', encoding='utf-8') as f:
            f.write("质量控制报告\n")
            f.write("=" * 20 + "\n")
            f.write(f"分析时间: {qc_results['analysis_time']}\n")
            f.write(f"文件: {input_file}\n")
            f.write(f"文件大小: {qc_results['file_size']} 字节\n")
            f.write(f"总行数: {qc_results['total_lines']}\n")
            f.write(f"空行数: {qc_results['empty_lines']}\n")
            f.write(f"格式错误行数: {qc_results['malformed_lines']}\n")
        
        # 写入数据验证文件
        validation_file = self.output_dir / "data_validation.txt"
        with open(validation_file, 'w', encoding='utf-8') as f:
            f.write("数据验证结果\n")
            f.write("=" * 20 + "\n")
            f.write(f"空行数: {qc_results['empty_lines']}\n")
            f.write(f"潜在格式错误行数: {qc_results['malformed_lines']}\n\n")
            f.write("样本数据 (前5行):\n")
            for i, line in enumerate(qc_results['sample_lines'], 1):
                f.write(f"{i}: {line}\n")
        
        self.results['qc'] = qc_results
        self.log_info(f"质量控制完成，结果保存到: {qc_report_file}, {validation_file}")
        return qc_results
    
    def generate_comprehensive_report(self):
        """生成综合分析报告"""
        self.log_info("生成综合分析报告...")
        
        # 生成文本报告
        report_file = self.output_dir / "transposon_analysis_report.txt"
        with open(report_file, 'w', encoding='utf-8') as f:
            f.write("转座子分析综合报告\n")
            f.write("=" * 50 + "\n")
            f.write(f"生成时间: {datetime.now().isoformat()}\n\n")
            
            # 组成分析结果
            if 'composition' in self.results:
                f.write("1. 组成分析结果\n")
                f.write("-" * 30 + "\n")
                comp_results = self.results['composition']
                f.write(f"总行数: {comp_results['total_lines']}\n")
                f.write("类别统计:\n")
                for key, value in sorted(comp_results['summary'].items()):
                    f.write(f"  {key}: {value}\n")
                f.write("\n")
            
            # GFF分析结果
            if 'gff' in self.results:
                f.write("2. GFF3分析结果\n")
                f.write("-" * 30 + "\n")
                gff_results = self.results['gff']
                f.write(f"总行数: {gff_results['total_lines']}\n")
                f.write(f"处理行数: {gff_results['processed_lines']}\n")
                f.write("类别计数:\n")
                for category, count in sorted(gff_results['category_counts'].items()):
                    f.write(f"  {category}: {count}\n")
                f.write("\n")
            
            # 库分析结果
            if 'library' in self.results:
                f.write("3. 库文件分析结果\n")
                f.write("-" * 30 + "\n")
                lib_results = self.results['library']
                f.write(f"总条目数: {lib_results['total_entries']}\n")
                f.write("染色体统计:\n")
                for chrom, stats in lib_results['chromosome_stats'].items():
                    f.write(f"  {chrom}: {stats['count']} 条目\n")
                f.write("\n")
            
            # QC结果
            if 'qc' in self.results:
                f.write("4. 质量控制结果\n")
                f.write("-" * 30 + "\n")
                qc_results = self.results['qc']
                f.write(f"文件大小: {qc_results['file_size']} 字节\n")
                f.write(f"总行数: {qc_results['total_lines']}\n")
                f.write(f"空行数: {qc_results['empty_lines']}\n")
                f.write(f"格式错误行数: {qc_results['malformed_lines']}\n")
                f.write("\n")
        
        # 生成HTML报告
        html_file = self.output_dir / "transposon_analysis_report.html"
        with open(html_file, 'w', encoding='utf-8') as f:
            f.write("""<!DOCTYPE html>
<html>
<head>
    <title>转座子分析报告</title>
    <meta charset="utf-8">
    <style>
        body { font-family: Arial, sans-serif; margin: 20px; }
        h1 { color: #2c3e50; }
        h2 { color: #34495e; border-bottom: 2px solid #ecf0f1; }
        pre { background-color: #f8f9fa; padding: 10px; border-radius: 5px; overflow-x: auto; }
        .header { background-color: #3498db; color: white; padding: 20px; border-radius: 5px; }
        .section { margin: 20px 0; }
        table { border-collapse: collapse; width: 100%; }
        th, td { border: 1px solid #ddd; padding: 8px; text-align: left; }
        th { background-color: #f2f2f2; }
    </style>
</head>
<body>
    <div class="header">
        <h1>转座子分析报告</h1>
        <p>生成时间: """ + datetime.now().strftime("%Y-%m-%d %H:%M:%S") + """</p>
    </div>
    
    <div class="section">
        <h2>分析概述</h2>
        <p>本报告包含转座子组成分析、GFF3类别统计、库文件分析和质量控制结果。</p>
    </div>
    
    <div class="section">
        <h2>详细结果</h2>
        <pre>""")
            
            # 读取文本报告内容
            with open(report_file, 'r', encoding='utf-8') as txt_file:
                f.write(txt_file.read())
            
            f.write("""</pre>
    </div>
</body>
</html>""")
        
        # 保存JSON格式的结果
        json_file = self.output_dir / "analysis_results.json"
        with open(json_file, 'w', encoding='utf-8') as f:
            json.dump(self.results, f, indent=2, ensure_ascii=False)
        
        self.log_info(f"综合报告生成完成: {report_file}, {html_file}, {json_file}")


def main():
    parser = argparse.ArgumentParser(description='转座子分析工具')
    parser.add_argument('--mode', choices=['composition', 'gff', 'library', 'qc', 'all'], 
                       default='all', help='分析模式')
    parser.add_argument('--input', required=True, help='主输入文件')
    parser.add_argument('--gff', help='GFF3文件 (可选)')
    parser.add_argument('--lib', help='库文件 (可选)')
    parser.add_argument('--output', default='results', help='输出目录')
    
    args = parser.parse_args()
    
    # 创建分析器
    analyzer = TransposonAnalyzer(args.output)
    
    try:
        # 根据模式执行相应分析
        if args.mode in ['composition', 'all']:
            analyzer.analyze_composition(args.input)
            analyzer.quality_control(args.input)
        
        if args.mode in ['gff', 'all'] and args.gff:
            analyzer.analyze_gff(args.gff)
        
        if args.mode in ['library', 'all'] and args.lib:
            analyzer.analyze_library(args.lib)
        
        if args.mode == 'qc':
            analyzer.quality_control(args.input)
        
        # 生成综合报告
        analyzer.generate_comprehensive_report()
        
        print(f"\n分析完成！结果保存在: {analyzer.output_dir}")
        
    except Exception as e:
        print(f"错误: {e}", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
