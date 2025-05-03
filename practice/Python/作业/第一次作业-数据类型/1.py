#定义各输出字段的宽度
width = input("请输入价格表总字符宽度：")
price_width = 10
############补充代码#############
width = int(width)
#################################
item_width = width - price_width

############补充代码#############
#定义格式化方式
head = "{{:<{}}}".format(item_width) + "{{:>{}}}".format(price_width)
body = "{{:<{}}}".format(item_width) + "{{:>{}.2f}}".format(price_width)
#################################

#打印价格列表
print("="*width)
print(head.format("Item","Price"))
print("-"*width)
print(body.format("Apple",4.0))
print(body.format("Watermelon",2.0))
print(body.format("Peach",5.0))
print("="*width)

