# FastCall2 Pipeline 故障排除指南

## 常见问题和解决方案

### 1. timeline.params.output_dir 错误

**问题**: 运行时出现"无法找到timeline.params.output_dir"错误

**原因**: 配置文件中的报告路径配置有问题

**解决方案**:
- 使用默认的 `nextflow.config` (已禁用报告功能)
- 或者使用 `nextflow.config.with_reports` 启用报告功能：
  ```bash
  nextflow run runFastCall2.nf -c nextflow.config.with_reports [其他参数]
  ```

### 2. Java 内存错误

**问题**: Java heap space 或 OutOfMemoryError

**解决方案**:
```bash
# 增加内存分配
nextflow run runFastCall2.nf --memory "200g" [其他参数]

# 减少线程数
nextflow run runFastCall2.nf --threads 16 [其他参数]
```

### 3. 找不到 BAM 文件

**问题**: Cannot find BAM files 或相关错误

**解决方案**:
1. 检查 taxaBamMap.txt 文件格式
2. 确保 BAM 文件路径正确
3. 确保 BAM 文件已建立索引 (.bai 文件)

### 4. 染色体命名不匹配

**问题**: 染色体名称在参考基因组和 BAM 文件中不一致

**解决方案**:
```bash
# 自定义染色体列表
nextflow run runFastCall2.nf --chromosomes "chr1A,chr1B,chr1D" [其他参数]
```

### 5. Nextflow 语法错误

**问题**: Nextflow script 语法错误

**解决方案**:
```bash
# 检查语法
nextflow run runFastCall2.nf --help

# 运行测试
test_pipeline.bat  # Windows
# 或
bash test_pipeline.sh  # Linux/Unix
```

### 6. 权限问题

**问题**: Permission denied 错误

**解决方案**:
```bash
# Windows (以管理员身份运行)
# Linux/Unix
chmod +x *.sh
sudo nextflow run runFastCall2.nf [参数]
```

### 7. 磁盘空间不足

**问题**: No space left on device

**解决方案**:
1. 清理临时文件
2. 更改工作目录到有更多空间的位置：
   ```bash
   nextflow run runFastCall2.nf -w /path/to/large/disk/work [其他参数]
   ```

### 8. TIGER jar 文件问题

**问题**: Could not find or load main class

**解决方案**:
1. 检查 TIGER jar 文件路径
2. 确保 Java 版本兼容
3. 测试 jar 文件：
   ```bash
   java -jar TIGER.jar
   ```

### 9. samtools 路径问题

**问题**: samtools: command not found

**解决方案**:
```bash
# 使用完整路径
nextflow run runFastCall2.nf --samtools_path /full/path/to/samtools [其他参数]

# 或添加到 PATH
export PATH=$PATH:/path/to/samtools/bin
```

### 10. 进程被杀死

**问题**: Process killed (exit status 137)

**原因**: 通常是内存不足

**解决方案**:
1. 减少并行处理的染色体数量
2. 增加系统内存
3. 使用交换文件

## 调试技巧

### 1. 查看详细日志
```bash
nextflow run runFastCall2.nf -with-report -with-timeline -with-trace [参数]
```

### 2. 保留工作目录
```bash
# 不自动清理工作目录，方便调试
nextflow run runFastCall2.nf -resume [参数]
```

### 3. 测试单个染色体
```bash
nextflow run runFastCall2.nf --chromosomes "1A" [其他参数]
```

### 4. 使用测试配置
```bash
nextflow run runFastCall2.nf -profile test [其他参数]
```

## 联系支持

如果问题仍然存在，请提供：
1. 完整的错误消息
2. 使用的命令
3. 系统信息 (操作系统、内存、Java版本)
4. nextflow.log 文件内容
