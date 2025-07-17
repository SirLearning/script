#!/usr/bin/env python3
"""
小麦GWAS分析脚本
适用于分析小麦基因组关联研究数据

使用方法:
python wheat_gwas_analysis_script.py

依赖:
- pandas
- numpy
- matplotlib
- seaborn
- scipy
- boto3 (用于访问数据湖)
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
import boto3
from pathlib import Path

class WheatGWASAnalyzer:
    """小麦GWAS分析器"""
    
    def __init__(self, output_dir='./gwas_results'):
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(exist_ok=True)
        self.gwas_data = None
        self.significant_snps = None
        
    def load_data(self, data_source='s3://biomni-datalake/gwas_catalog.pkl'):
        """加载GWAS数据"""
        print(f"正在加载数据: {data_source}")
        try:
            self.gwas_data = pd.read_pickle(data_source)
            print(f"数据加载成功，形状: {self.gwas_data.shape}")
            return True
        except Exception as e:
            print(f"数据加载失败: {e}")
            return False
    
    def filter_wheat_data(self):
        """筛选小麦相关数据"""
        if self.gwas_data is None:
            print("请先加载数据")
            return False
            
        # 筛选小麦相关的GWAS数据
        wheat_mask = self.gwas_data['DISEASE/TRAIT'].str.contains('wheat', case=False, na=False)
        self.wheat_gwas = self.gwas_data[wheat_mask]
        
        print(f"找到 {len(self.wheat_gwas)} 个小麦相关的GWAS记录")
        return True
    
    def quality_control(self, pvalue_threshold=5e-8):
        """数据质量控制"""
        if not hasattr(self, 'wheat_gwas'):
            print("请先筛选小麦数据")
            return False
            
        # 筛选显著性SNP
        self.significant_snps = self.wheat_gwas[self.wheat_gwas['P-VALUE'] < pvalue_threshold]
        
        print(f"质量控制结果:")
        print(f"- 显著性SNP数量 (P < {pvalue_threshold}): {len(self.significant_snps)}")
        print(f"- 涉及染色体数量: {self.wheat_gwas['CHR_ID'].nunique()}")
        print(f"- 涉及研究数量: {self.wheat_gwas['PUBMEDID'].nunique()}")
        
        return True
    
    def create_manhattan_plot(self):
        """创建曼哈顿图"""
        if not hasattr(self, 'wheat_gwas'):
            print("请先筛选小麦数据")
            return False
            
        plt.figure(figsize=(12, 6))
        
        # 按染色体分组绘制
        chr_colors = plt.cm.tab10(np.linspace(0, 1, self.wheat_gwas['CHR_ID'].nunique()))
        chr_ids = self.wheat_gwas['CHR_ID'].unique()
        
        for i, chr_id in enumerate(chr_ids):
            chr_data = self.wheat_gwas[self.wheat_gwas['CHR_ID'] == chr_id]
            plt.scatter(chr_data['CHR_POS'], -np.log10(chr_data['P-VALUE']), 
                       c=[chr_colors[i]], label=f'Chr {chr_id}', alpha=0.7)
        
        # 添加显著性线
        plt.axhline(y=-np.log10(5e-8), color='red', linestyle='--', alpha=0.7, 
                   label='Genome-wide significance')
        
        plt.xlabel('Position')
        plt.ylabel('-log10(P-value)')
        plt.title('Manhattan Plot - Wheat GWAS Results')
        plt.legend()
        plt.grid(True, alpha=0.3)
        
        plot_path = self.output_dir / 'manhattan_plot.png'
        plt.tight_layout()
        plt.savefig(plot_path, dpi=300, bbox_inches='tight')
        plt.close()
        
        print(f"曼哈顿图已保存至: {plot_path}")
        return True
    
    def create_qq_plot(self):
        """创建QQ图"""
        if not hasattr(self, 'wheat_gwas'):
            print("请先筛选小麦数据")
            return False
            
        # 计算观测和期望P值
        observed_p = self.wheat_gwas['P-VALUE'].dropna().sort_values()
        n = len(observed_p)
        expected_p = np.arange(1, n+1) / (n+1)
        
        # 计算-log10值
        observed_log = -np.log10(observed_p)
        expected_log = -np.log10(expected_p)
        
        # 创建QQ图
        plt.figure(figsize=(8, 8))
        plt.scatter(expected_log, observed_log, alpha=0.6, color='blue')
        plt.plot([0, max(expected_log)], [0, max(expected_log)], 'r--', alpha=0.8, label='Expected')
        plt.xlabel('Expected -log10(P)')
        plt.ylabel('Observed -log10(P)')
        plt.title('QQ Plot - Wheat GWAS Results')
        plt.legend()
        plt.grid(True, alpha=0.3)
        
        # 计算lambda值
        lambda_gc = np.median(stats.chi2.ppf(1 - observed_p, df=1)) / stats.chi2.ppf(0.5, df=1)
        plt.text(0.05, 0.95, f'λ = {lambda_gc:.3f}', transform=plt.gca().transAxes, 
                bbox=dict(boxstyle="round,pad=0.3", facecolor="white", alpha=0.8))
        
        plot_path = self.output_dir / 'qq_plot.png'
        plt.tight_layout()
        plt.savefig(plot_path, dpi=300, bbox_inches='tight')
        plt.close()
        
        print(f"QQ图已保存至: {plot_path}")
        print(f"基因组膨胀因子(λ): {lambda_gc:.3f}")
        return True
    
    def save_results(self):
        """保存分析结果"""
        if not hasattr(self, 'wheat_gwas'):
            print("请先筛选小麦数据")
            return False
            
        # 保存完整数据
        wheat_data_file = self.output_dir / 'wheat_gwas_data.csv'
        self.wheat_gwas.to_csv(wheat_data_file, index=False)
        
        # 保存显著性SNP
        if hasattr(self, 'significant_snps') and len(self.significant_snps) > 0:
            sig_file = self.output_dir / 'significant_snps.csv'
            self.significant_snps.to_csv(sig_file, index=False)
            print(f"显著性SNP已保存至: {sig_file}")
        
        print(f"小麦GWAS数据已保存至: {wheat_data_file}")
        return True
    
    def run_analysis(self):
        """运行完整分析"""
        print("=== 开始小麦GWAS分析 ===")
        
        # 加载数据
        if not self.load_data():
            return False
            
        # 筛选小麦数据
        if not self.filter_wheat_data():
            return False
            
        # 质量控制
        if not self.quality_control():
            return False
            
        # 创建图表
        self.create_manhattan_plot()
        self.create_qq_plot()
        
        # 保存结果
        self.save_results()
        
        print("=== 小麦GWAS分析完成 ===")
        return True

if __name__ == "__main__":
    # 创建分析器实例
    analyzer = WheatGWASAnalyzer(output_dir='./wheat_gwas_results')
    
    # 运行分析
    analyzer.run_analysis()
