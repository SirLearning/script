#!/usr/bin/env python3
"""
小麦种质Taxonomy信息提取工具

该工具连接NCBI Taxonomy数据库，根据输入的物种名称提取taxonomy信息，
包括taxon ID、同义词、俗名、分类等级、谱系等信息。

作者: 生物信息学研究助手
版本: 1.0
"""

import pandas as pd
import requests
import json
import time
import os
import xml.etree.ElementTree as ET
from typing import Dict, List, Optional, Union
from urllib.parse import quote
import argparse


class WheatGermplasmTaxonomyExtractor:
    """
    小麦种质taxonomy信息提取器
    连接NCBI Taxonomy数据库，提取完整的taxonomy信息
    """
    
    def __init__(self, output_dir: str = './output'):
        """
        初始化提取器
        
        Args:
            output_dir: 输出目录路径
        """
        self.output_dir = output_dir
        self.base_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
        self.session = requests.Session()
        self.session.headers.update({
            'User-Agent': 'WheatGermplasmTaxonomyExtractor/1.0'
        })
        
        # 创建输出目录
        os.makedirs(self.output_dir, exist_ok=True)
        
        # 存储处理结果
        self.results = []
        
        # 小麦相关的常见属名
        self.wheat_genera = [
            'Triticum', 'Aegilops', 'Secale', 'Hordeum', 'Agropyron',
            'Elymus', 'Dasypyrum', 'Pseudoroegneria', 'Thinopyrum'
        ]
    
    def search_taxon_by_name(self, name: str, delay: float = 0.5) -> Dict:
        """
        根据名称搜索taxon ID
        
        Args:
            name: 物种名称（学名或俗名）
            delay: 请求间隔时间（秒）
            
        Returns:
            包含搜索结果的字典
        """
        try:
            # 添加延迟以避免请求过于频繁
            time.sleep(delay)
            
            # 构建搜索URL
            search_url = f"{self.base_url}esearch.fcgi"
            params = {
                'db': 'taxonomy',
                'term': name,
                'retmode': 'json',
                'retmax': 20
            }
            
            # 发送请求
            response = self.session.get(search_url, params=params, timeout=30)
            response.raise_for_status()
            
            data = response.json()
            
            if 'esearchresult' in data and 'idlist' in data['esearchresult']:
                id_list = data['esearchresult']['idlist']
                count = int(data['esearchresult']['count'])
                
                return {
                    'query': name,
                    'count': count,
                    'taxon_ids': id_list,
                    'status': 'success'
                }
            else:
                return {
                    'query': name,
                    'count': 0,
                    'taxon_ids': [],
                    'status': 'no_results'
                }
                
        except Exception as e:
            return {
                'query': name,
                'count': 0,
                'taxon_ids': [],
                'status': 'error',
                'error': str(e)
            }
    
    def fetch_taxon_details(self, taxon_id: str, delay: float = 0.5) -> Dict:
        """
        获取指定taxon ID的详细信息
        
        Args:
            taxon_id: NCBI taxon ID
            delay: 请求间隔时间（秒）
            
        Returns:
            包含详细信息的字典
        """
        try:
            # 添加延迟
            time.sleep(delay)
            
            # 构建获取URL
            fetch_url = f"{self.base_url}efetch.fcgi"
            params = {
                'db': 'taxonomy',
                'id': taxon_id,
                'retmode': 'xml'
            }
            
            # 发送请求
            response = self.session.get(fetch_url, params=params, timeout=30)
            response.raise_for_status()
            
            # 解析XML
            root = ET.fromstring(response.content)
            
            # 提取信息
            taxon_info = self._parse_taxon_xml(root, taxon_id)
            
            return taxon_info
            
        except Exception as e:
            return {
                'taxon_id': taxon_id,
                'status': 'error',
                'error': str(e)
            }
    
    def _parse_taxon_xml(self, root: ET.Element, taxon_id: str) -> Dict:
        """
        解析taxonomy XML数据
        
        Args:
            root: XML根元素
            taxon_id: taxon ID
            
        Returns:
            解析后的信息字典
        """
        result = {
            'taxon_id': taxon_id,
            'scientific_name': '',
            'common_names': [],
            'synonyms': [],
            'rank': '',
            'lineage': [],
            'parent_taxon_id': '',
            'genetic_code': '',
            'mitochondrial_genetic_code': '',
            'status': 'success'
        }
        
        # 查找Taxon元素
        taxon = root.find('.//Taxon')
        if taxon is None:
            result['status'] = 'no_data'
            return result
        
        try:
            # 提取基本信息
            elements = {
                'taxon_id': taxon.find('TaxId'),
                'scientific_name': taxon.find('ScientificName'),
                'rank': taxon.find('Rank'),
                'parent_taxon_id': taxon.find('ParentTaxId')
            }
            
            for key, elem in elements.items():
                if elem is not None and elem.text:
                    result[key] = elem.text
            
            # 提取遗传密码信息
            genetic_code_elem = taxon.find('.//GeneticCode/GCId')
            if genetic_code_elem is not None and genetic_code_elem.text:
                result['genetic_code'] = genetic_code_elem.text
            
            mitochondrial_genetic_code_elem = taxon.find('.//MitoGeneticCode/MGCId')
            if mitochondrial_genetic_code_elem is not None and mitochondrial_genetic_code_elem.text:
                result['mitochondrial_genetic_code'] = mitochondrial_genetic_code_elem.text
            
            # 提取其他名称
            other_names = taxon.find('OtherNames')
            if other_names is not None:
                # 俗名
                for common_name in other_names.findall('.//CommonName'):
                    if common_name.text:
                        result['common_names'].append(common_name.text)
                
                # 同义词
                for synonym in other_names.findall('.//Synonym'):
                    if synonym.text:
                        result['synonyms'].append(synonym.text)
            
            # 提取谱系信息
            lineage_elem = taxon.find('Lineage')
            if lineage_elem is not None and lineage_elem.text:
                result['lineage'] = [name.strip() for name in lineage_elem.text.split(';') if name.strip()]
            
        except Exception as e:
            result['status'] = 'parse_error'
            result['error'] = str(e)
        
        return result
    
    def batch_process(self, species_names: List[str], delay: float = 0.5) -> List[Dict]:
        """
        批量处理物种名称列表
        
        Args:
            species_names: 物种名称列表
            delay: 请求间隔时间（秒）
            
        Returns:
            处理结果列表
        """
        results = []
        
        print(f"开始批量处理 {len(species_names)} 个物种...")
        
        for i, name in enumerate(species_names, 1):
            print(f"[{i}/{len(species_names)}] 处理物种: {name}")
            
            # 搜索taxon ID
            search_result = self.search_taxon_by_name(name, delay)
            
            if search_result['status'] == 'success' and search_result['count'] > 0:
                # 获取第一个结果的详细信息
                taxon_id = search_result['taxon_ids'][0]
                details = self.fetch_taxon_details(taxon_id, delay)
                
                # 合并结果
                result = {
                    'input_name': name,
                    'search_count': search_result['count'],
                    'all_taxon_ids': search_result['taxon_ids'],
                    **details
                }
                
                if details['status'] == 'success':
                    print(f"  ✓ 成功: {details['scientific_name']} (ID: {details['taxon_id']})")
                else:
                    print(f"  ✗ 获取详细信息失败: {details.get('error', '未知错误')}")
                
            else:
                result = {
                    'input_name': name,
                    'search_count': 0,
                    'all_taxon_ids': [],
                    'status': search_result['status'],
                    'error': search_result.get('error', '未找到结果')
                }
                
                print(f"  ✗ 搜索失败: {result['error']}")
            
            results.append(result)
        
        self.results = results
        return results
    
    def export_results(self, filename: str = 'wheat_taxonomy_results', formats: List[str] = ['csv', 'json']) -> List[str]:
        """
        导出结果到文件
        
        Args:
            filename: 文件名前缀
            formats: 导出格式列表
            
        Returns:
            导出文件路径列表
        """
        if not self.results:
            print("没有结果可导出")
            return []
        
        exported_files = []
        
        for fmt in formats:
            if fmt == 'csv':
                # 导出CSV
                csv_file = os.path.join(self.output_dir, f"{filename}.csv")
                
                # 准备CSV数据
                csv_data = []
                for result in self.results:
                    csv_row = {
                        'input_name': result.get('input_name', ''),
                        'taxon_id': result.get('taxon_id', ''),
                        'scientific_name': result.get('scientific_name', ''),
                        'rank': result.get('rank', ''),
                        'parent_taxon_id': result.get('parent_taxon_id', ''),
                        'genetic_code': result.get('genetic_code', ''),
                        'mitochondrial_genetic_code': result.get('mitochondrial_genetic_code', ''),
                        'common_names': '; '.join(result.get('common_names', [])),
                        'synonyms': '; '.join(result.get('synonyms', [])),
                        'lineage': '; '.join(result.get('lineage', [])),
                        'search_count': result.get('search_count', 0),
                        'all_taxon_ids': '; '.join(result.get('all_taxon_ids', [])),
                        'status': result.get('status', '')
                    }
                    csv_data.append(csv_row)
                
                df = pd.DataFrame(csv_data)
                df.to_csv(csv_file, index=False, encoding='utf-8')
                exported_files.append(csv_file)
                print(f"CSV文件已保存: {csv_file}")
            
            elif fmt == 'json':
                # 导出JSON
                json_file = os.path.join(self.output_dir, f"{filename}.json")
                
                with open(json_file, 'w', encoding='utf-8') as f:
                    json.dump(self.results, f, ensure_ascii=False, indent=2)
                
                exported_files.append(json_file)
                print(f"JSON文件已保存: {json_file}")
        
        return exported_files
    
    def get_wheat_species_suggestions(self) -> List[str]:
        """
        获取小麦相关物种的建议列表
        
        Returns:
            小麦相关物种名称列表
        """
        return [
            "Triticum aestivum",          # 普通小麦
            "Triticum durum",             # 硬粒小麦
            "Triticum monococcum",        # 一粒小麦
            "Triticum turgidum",          # 圆锥小麦
            "Triticum spelta",            # 斯佩尔特小麦
            "Triticum compactum",         # 密穗小麦
            "Triticum carthlicum",        # 格鲁吉亚小麦
            "Triticum timopheevii",       # 提莫菲小麦
            "Aegilops tauschii",          # 节节麦
            "Aegilops speltoides",        # 拟斯佩尔特山羊草
            "Aegilops searsii",           # 西尔斯山羊草
            "Secale cereale",             # 黑麦
            "Hordeum vulgare",            # 大麦
            "Hordeum bulbosum",           # 球茎大麦
            "Elymus repens",              # 偃麦草
            "Agropyron cristatum",        # 冰草
            "Dasypyrum villosum",         # 簇毛麦
            "Pseudoroegneria spicata",    # 刺冰草
            "Thinopyrum elongatum"        # 长穗偃麦草
        ]


def main():
    """主函数 - 命令行接口"""
    parser = argparse.ArgumentParser(description='小麦种质Taxonomy信息提取工具')
    parser.add_argument('--input', '-i', type=str, help='输入文件路径（每行一个物种名称）')
    parser.add_argument('--output', '-o', type=str, default='./output', help='输出目录路径')
    parser.add_argument('--species', '-s', type=str, nargs='+', help='直接指定物种名称列表')
    parser.add_argument('--delay', '-d', type=float, default=0.5, help='请求间隔时间（秒）')
    parser.add_argument('--format', '-f', type=str, nargs='+', default=['csv', 'json'], 
                        choices=['csv', 'json'], help='导出格式')
    parser.add_argument('--examples', action='store_true', help='显示小麦相关物种的建议列表')
    
    args = parser.parse_args()
    
    # 创建提取器
    extractor = WheatGermplasmTaxonomyExtractor(output_dir=args.output)
    
    # 显示示例
    if args.examples:
        suggestions = extractor.get_wheat_species_suggestions()
        print("小麦相关物种建议列表:")
        for i, species in enumerate(suggestions, 1):
            print(f"{i:2d}. {species}")
        return
    
    # 获取物种名称列表
    species_names = []
    
    if args.species:
        species_names = args.species
    elif args.input:
        with open(args.input, 'r', encoding='utf-8') as f:
            species_names = [line.strip() for line in f if line.strip()]
    else:
        print("请提供物种名称列表（--species）或输入文件（--input）")
        print("使用 --examples 查看建议的小麦相关物种")
        return
    
    if not species_names:
        print("没有找到要处理的物种名称")
        return
    
    # 批量处理
    results = extractor.batch_process(species_names, delay=args.delay)
    
    # 导出结果
    exported_files = extractor.export_results(formats=args.format)
    
    # 显示统计信息
    successful_count = len([r for r in results if r.get('status') == 'success'])
    print(f"\n处理完成！")
    print(f"- 总数: {len(results)}")
    print(f"- 成功: {successful_count}")
    print(f"- 失败: {len(results) - successful_count}")
    print(f"- 导出文件: {exported_files}")


if __name__ == "__main__":
    main()
