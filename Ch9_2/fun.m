% Computes the differences between forward model prediction and data,
% normalized by the standard deviation for the slug test Example 9.1.
% from Parameter Estimation and Inverse Problems, 3rd edition, 2018
% by R. Aster, B. Borchers, C. Thurber
%
% fvec=fun(p)
%
%
% p is expected to be a five element vector
%   [a, n, m, thetaS, thetaR]
function fvec=fun(p)
% global variables, these are 
global theta
global h
global sigma
% Compute the function values.
fvec=zeros(length(theta),1);
for i=1:length(theta)
  fvec(i)= (p(5) + (p(4) - p(5))./((1+(p(1).*abs(h(i)).^p(2))).^p(3))-theta(i))./sigma;
end  
