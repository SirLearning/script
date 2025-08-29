multi.cov.test = function(data, ind, r) {
  n = nrow(data); p = ncol(data)-1; data = data[, 1:p]
  V = 0
  for (i in 1:r) {
    datai = data[ind==i, ]; ni = nrow(datai)
    V = V + (ni-1)*cov(datai)
  }
  det.V = 0
  for (i in 1:r) {
    datai = data[ind==i, ]; ni = nrow(datai)
    det.V = det.V + (ni-1)*log(det(cov(datai)))
  }
  M = (n-r)*log(det(V/(n-r))) - det.V
  nu = (2*p^2+3*p-1)*(r+1)/(6*(p+1)*(n-r))
  f = p*(p+1)*(r-1)/2; T = (1-nu)*M
  p.value = 1 - pchisq(T, f)
  return(p.value = p.value)
}
Y = factor(plant_height$country)
multi.cov.test(plant_height, ind = Y, r = 3)
