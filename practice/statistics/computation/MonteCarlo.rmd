---
title: "R Notebook"
output: html_notebook
---

# Examples

1. count $\sum\cos(x)dx$
```{r}
result <- vector("numeric")
while (length(result) < 10000) {
  y <- runif(1, 0, pi/2)
  u <- runif(1)
  if (u <= cos(y)) {
    result <- append(result, y)
  }
}
judge1 <- result <= 1
judge2 <- result >= 0
m <- sum(judge1*judge2)
```

2. count $\int_0^1\cos(x)dx$
```{r}
x <- result[judge1]
gx <- exp(x)*cos(x)
mean(gx)
```

3. count by `mosaicCalc`
```{r}
library(mosaicCalc)
g <- antiD(e^x*cos(x)~x)
g(1, e = exp(1)) - g(0, e = exp(1))
```

4. integrate $\int_0^3\cos(exp(x))dx$
```{r}
x <- runif(10000, 0, 3)
gx <- cos(exp(x))
3*mean(gx)

library(mosaicCalc)
g <- antiD(cos(exp(x))~x)
g(3)-g(0)
```




