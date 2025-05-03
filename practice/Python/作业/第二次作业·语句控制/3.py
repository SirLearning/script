e = 1
factNum = 1
for n in range(1,10**5):
    factNum = factNum * n
    e += 1/factNum
print(e)


