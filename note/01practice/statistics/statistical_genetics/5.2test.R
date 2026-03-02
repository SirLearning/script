alpha=.05; df=1;
c=qchisq(alpha, df, ncp = 0, lower.tail = F)

T=1.645; df=1;
prob=pchisq(T, df, ncp = 0, lower.tail = F)

lambda=4.5; alpha=.01; df=3;
c=pchisq(c, df=df, ncp=lambda, lower.tail = F):beta=1-power

n1=1000; lambda1=1.551; alpha=0.05; df=1;
n2=4000; lambda2=n2*(lambda1/n1);
c=qchisq(alpha, df, lower.tail = F);
power1=pchisq(c, df=df, ncp = lambda1, lower.tail = F);
power2=pchisq(c, df=df, ncp = lambda2, lower.tail = F)
