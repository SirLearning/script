awk [-F re] [parameter...] ['prog'] [-f progfile] [in_file...]
# -F re: re is a regular expression, which is used to specify the field separator. The default field separator is a space or a tab.
# parameter: The parameter is used to control the input and output of the program. The parameter can be one of the following:
#   -f: Read the program from the file progfile.
#   -F: Specify the field separator.
#   -v: Assign a value to a variable.
#   -W: Specify the maximum number of fields.
#   -w: Specify the maximum number of characters in a string.
#   -mf: Specify the maximum number of fields in a record.
#   -mr: Specify the maximum number of records in a file.
# -v var=value: Assign the value value to the variable var before the program begins.
# 'prog': The program is enclosed in single quotes. The program can be a series of patterns and corresponding actions.
# 'pattern {action}': The pattern is enclosed in braces, and the action is enclosed in braces. The action is performed when the pattern is matched.
# -f progfile: Read the program from the file progfile. This is useful when the program is long.
# in_file: The input file. If in_file is not specified, awk reads from the standard input.

# 处理记录（文本文件中的一行）、字段（记录中的一个单词或者一个字符串）和变量（用于存储数据的名称）。
# 内置变量：
#   $0表示整个记录，$1表示第一个字段，$2表示第二个字段，以此类推。分隔符默认是空格，内置变量FS表示该分隔符。
#   NF表示字段的数量，NR表示记录的数量，RS表示记录的分隔符，OFS表示输出字段的分隔符，ORS表示输出记录的分隔符。
#   FILENAME表示当前文件的名称，FNR表示当前文件的记录数。
#   BEGIN和END是两个特殊的模式，分别在处理开始和结束时执行。
# example 1
awk -F % 'NR==7,NR==15 {printf $1 $3 $7}'

# 内置函数：
#   length(s)：返回字符串s的长度。
# example 2
awk '{print"%03d%s\n",NR,$1}' myfile