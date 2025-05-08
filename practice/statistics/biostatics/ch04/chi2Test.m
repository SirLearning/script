function [p, Q]= chi2test(x)
% Check the arguments.
if(nargin ~= 1), error('One and only one argument required!');
end
if(ndims(x) ~= 2), error('The argument (x) must be a 2d matrix!');
end
if(any(size(x) == 1)), error('The argument (x) must be a 2d matrix!'); end
if(any(~isreal(x))), error('All values of the argument (x) must be real values!'); end
% Calculate Q = sum( (a-np*)^2/(np*(1-p*)) )
s= size(x, 1);
r= size(x, 2);
np= sum(x, 2)/sum(sum(x)) * sum(x); % p=sum(x, 2)/sum(sum(x)) and n=sum(x)
Q= sum(sum((x-np).^2./(np)));
% Calculate cdf of chi-squared to Q. Degrees of freedom, v, is (r-1)*(s-1).
p= 1 - gammainc(Q/2, (r-1)*(s-1)/2);
End