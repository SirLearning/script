def primeJudge(n):
    if n <= 1:
        return False
    for i in range(2, int(n**0.5)+1):
        if n%i == 0:
            return False
    return True

primes = []
num = 1000
while len(primes) < 20 and num >= 1:
    if primeJudge(num):
        primes.append(num)
    num -= 1
    
print(sum(primes))





