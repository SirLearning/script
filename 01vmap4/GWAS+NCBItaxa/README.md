# 小麦种质Taxonomy信息提取工具

这是一个专门用于从NCBI Taxonomy数据库提取小麦种质taxonomy信息的Python工具。该工具可以根据输入的物种名称自动提取taxon ID、同义词、俗名、分类等级、谱系等详细信息。

## 功能特点

- **自动查询**: 连接NCBI Taxonomy数据库进行实时查询
- **智能匹配**: 支持学名、俗名和同义词的智能匹配
- **批量处理**: 支持大批量物种名称的自动化处理
- **多格式导出**: 支持CSV和JSON格式的结果导出
- **错误处理**: 完善的错误处理和重试机制
- **小麦专用**: 预置了小麦相关物种的建议列表

## 安装依赖

```bash
pip install pandas requests
```

## 使用方法

### 1. 基本用法

```bash
# 查询单个物种
python wheat_taxonomy_extractor.py --species "Triticum aestivum"

# 查询多个物种
python wheat_taxonomy_extractor.py --species "Triticum aestivum" "Triticum durum" "Aegilops tauschii"
```

### 2. 从文件批量处理

```bash
# 从文件读取物种名称列表
python wheat_taxonomy_extractor.py --input species_list.txt --output ./results
```

创建输入文件 `species_list.txt`:
```
Triticum aestivum
Triticum durum
Triticum monococcum
Aegilops tauschii
Secale cereale
```

### 3. 查看建议物种列表

```bash
python wheat_taxonomy_extractor.py --examples
```

### 4. 自定义参数

```bash
# 设置请求间隔和导出格式
python wheat_taxonomy_extractor.py \
    --species "Triticum aestivum" "Triticum durum" \
    --delay 1.0 \
    --format csv json \
    --output ./my_results
```

## 输出格式

### CSV格式
包含以下字段：
- `input_name`: 输入的物种名称
- `taxon_id`: NCBI taxon ID
- `scientific_name`: 标准学名
- `rank`: 分类等级
- `parent_taxon_id`: 父级taxon ID
- `genetic_code`: 遗传密码
- `mitochondrial_genetic_code`: 线粒体遗传密码
- `common_names`: 俗名（分号分隔）
- `synonyms`: 同义词（分号分隔）
- `lineage`: 完整谱系（分号分隔）
- `search_count`: 搜索结果数量
- `all_taxon_ids`: 所有找到的taxon ID
- `status`: 处理状态

### JSON格式
包含完整的结构化数据，适用于程序化处理。

## 预置小麦相关物种

工具预置了以下小麦相关物种：

1. **Triticum属（小麦属）**
   - Triticum aestivum（普通小麦）
   - Triticum durum（硬粒小麦）
   - Triticum monococcum（一粒小麦）
   - Triticum turgidum（圆锥小麦）
   - Triticum spelta（斯佩尔特小麦）
   - Triticum compactum（密穗小麦）

2. **Aegilops属（山羊草属）**
   - Aegilops tauschii（节节麦）
   - Aegilops speltoides（拟斯佩尔特山羊草）
   - Aegilops searsii（西尔斯山羊草）

3. **相关属**
   - Secale cereale（黑麦）
   - Hordeum vulgare（大麦）
   - Elymus repens（偃麦草）
   - Agropyron cristatum（冰草）

## 作为Python模块使用

```python
from wheat_taxonomy_extractor import WheatGermplasmTaxonomyExtractor

# 创建提取器
extractor = WheatGermplasmTaxonomyExtractor(output_dir='./results')

# 查询单个物种
result = extractor.search_taxon_by_name("Triticum aestivum")
details = extractor.fetch_taxon_details(result['taxon_ids'][0])

# 批量处理
species_list = ["Triticum aestivum", "Triticum durum", "Aegilops tauschii"]
results = extractor.batch_process(species_list)

# 导出结果
exported_files = extractor.export_results(formats=['csv', 'json'])
```

## 注意事项

1. **API限制**: 该工具使用NCBI的公共API，请遵守使用规范，避免过于频繁的请求
2. **网络连接**: 需要稳定的网络连接以访问NCBI数据库
3. **物种名称**: 建议使用标准的学名进行查询，以获得最准确的结果
4. **请求间隔**: 默认请求间隔为0.5秒，可根据需要调整

## 错误处理

工具包含完善的错误处理机制：
- 网络连接错误会自动重试
- 解析错误会记录在结果中
- 查询失败会在输出中标明状态

## 技术细节

- 使用NCBI E-utilities API进行查询
- XML解析提取详细的taxonomy信息
- 支持异步处理以提高效率
- 完整的日志记录和状态跟踪

## 示例输出

```csv
input_name,taxon_id,scientific_name,rank,common_names,synonyms,status
Triticum aestivum,4565,Triticum aestivum,species,"Canadian hard winter wheat; common wheat","Triticum aestivum subsp. aestivum; Triticum sativum",success
Triticum durum,4567,Triticum turgidum subsp. durum,subspecies,"","Triticum durum; Triticum durum ssp. durum",success
```

## 支持和贡献

如有问题或建议，请联系开发者。
