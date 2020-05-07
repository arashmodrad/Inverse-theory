% Computes the Jacobian of fun for the slug test Example 9.1
% from Parameter Estimation and Inverse Problems, 3rd edition, 2018
% by R. Aster, B. Borch(ii)ers, C. Th(ii)urber
% 
% p is expected to be a five element vector
%   [a, n, m, thetaS, thetaR]
%
function J=jac(p)
% global variables, these are 
% theta, h
global theta;
global h;
global sigma;

% use known formula for the derivatives in the Jacobian
nn=length(theta);
J=zeros(nn,2);
for ii=1:nn
    % Recalculated with Weights
    J(ii,1) = -(p(2).*p(3).*abs(h(ii)).*(p(4)-p(5)).*(p(1).*abs(h(ii))).^(p(2)-1).*((p(1).*abs(h(ii))).^(p(2)) + 1).^(-p(3)-1))./sigma;
    J(ii,2) = -(p(3).*log(p(1).*abs(h(ii))).*(p(4)-p(5)).*(p(1).*abs(h(ii))).^(p(2)).*((p(1).*abs(h(ii))).^(p(2)) + 1).^(-p(3)-1))./sigma;
    J(ii,3) = -(log((p(1).*abs(h(ii))).^(p(2))+1).*(p(4)-p(5)).*((p(1).*abs(h(ii))).^(p(2))+1).^(-p(3)))./sigma;
    J(ii,4) = 1./(sigma.*((p(1).*abs(h(ii))).^(p(2))+1).^(p(3)));
    J(ii,5) = (1./sigma) - 1./(sigma.*((p(1).*abs(h(ii))).^(p(2))+1).^(p(3)));

% Orignal Form
%   J(ii,1)= (p(4)-p(5)).*(-p(1).*h(ii).^2.*p(3).*(p(2)).*abs(p(1).*h(ii)).^(p(2)-2)...
%       .*(abs(p(1).*h(ii)).^(p(2))+1).^(-p(3)-1));
%   J(ii,2)= (p(4)-p(5)).*(-p(3).*abs(p(1).*h(ii)).^(p(2)).*log(abs(p(1).*h(ii)))...
%   .*(abs(p(1).*h(ii)).^(p(2))+1).^(-p(3)-1));
%   J(ii,3)= (p(4)-p(5)).*(-(abs(p(1).*h(ii)).^(p(2))+1).^(-p(3)).*log(p(1).*h(ii)).^(p(2))+1);
%   J(ii,4) = 1./((1+abs(p(1).*h(ii)).^p(2)).^p(3));
%   J(ii,5) = 1 - 1./((1+abs(p(1).*h(ii)).^p(2)).^p(3));
end
