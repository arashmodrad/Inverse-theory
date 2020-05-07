%% Group Assignment Chapter 9
clear; close all; clc;
%% 
global theta;
global h;
global sigma;
[h] = xlsread('soil.xlsx','sheet1','B2:B29');
theta = xlsread('soil.xlsx','sheet1','C2:C29');
sigma = std(theta);
% sigma = 1;
h(1) = 0.00001;
% p0 = ([1,1,1,.6,.1]');
% Parameters from table
% m = 1/(exp(.182)-1)
% p0 = [exp(-2.076),exp(.182),1,.482,.090]';
% Guess with Weighting
p0 = ([.01,1,1,1,1]')./2;

tol = eps.^(.125);
maxiter = 1000;

[pstar, iter] = lm('fun', 'jac', p0, 0.0001, maxiter);

goodWill =  vanGenuchten(pstar);

% Posterior Standard Deviation Estimate
dof = length(h) - 5;
s = norm(goodWill - theta)./sqrt(dof);

figure();plot((goodWill),(h),'k','linewidth',2);
hold on; errorbar((theta),(h),(ones(length(theta),1).*s./2),'horizontal','or','markerfacecolor','r','markersize',5);
xlabel('Volumetric Water Content (%)')
ylabel('Pressure Head (cm w.e.)')
set(gca,'fontweight','bold')
% figure();plot(log(goodWill),log(h),'k','linewidth',2);
% hold on; errorbar(log(theta),log(h),(ones(length(theta),1).*s./2),'horizontal','or','markerfacecolor','r');

% Report the chi2 Value
chi2 = sum(((goodWill - theta)./s).^2);
% Report the p value
p = chi2cdf(chi2,dof);

% Compute the Covariance and Correlation Matrecies
J=jac(pstar);
covm = s.^2.*inv(J'*J);
sig = sqrt(diag(covm));
corrm = covm./(sig*sig');

%% Random Starting Guess
pstar = [];
for ii = 1:50
p0 = [rand, rand, rand, rand, rand]';

tol = eps.^(.125);
maxiter = 1000;

[tmp, iter(ii)] = lm('fun', 'jac', p0, tol, maxiter);
if isreal(tmp)
    pstar = [pstar,tmp];
end
    
% pstar(:,ii) = real(pstar(:,11));
goodWill =  vanGenuchten(pstar(:,ii));

% Posterior Standard Deviation Estimate
dof = length(h) - 5;
s(ii) = norm(goodWill - theta)./sqrt(dof);
end
hmm = abs(pstar);
%% Chi Squared Objective Function
n = 500;
testVals = [linspace(0,.25,n);linspace(0.5,1.5,n);linspace(0,1,n); linspace(0,1,n); linspace(0,.5,n)]';
% Loop over a and n
X1 = zeros(n);
for ii = 1:n
    for jj = 1:n
%        tmp = vanGenuchten([testVals(ii,1),testVals(jj,2),testVals(500,3),testVals(500,4),testVals(500,5)]);
       tmp = vanGenuchten([testVals(ii,1),testVals(jj,2),pstar(3:5)']);

       X1(jj,ii) = sum(((tmp-theta)./s).^2);
    end
end

% Loop over n and m
X2 = zeros(n);
for ii = 1:n
    for jj = 1:n
%        tmp = vanGenuchten([testVals(ii,1),testVals(jj,2),testVals(500,3),testVals(500,4),testVals(500,5)]);
       tmp = vanGenuchten([pstar(1),testVals(ii,2),testVals(jj,3),pstar(4:5)']);

       X2(jj,ii) = sum(((tmp-theta)./s).^2);
    end
end

% Loop over m and thetas
X3 = zeros(n);
for ii = 1:n
    for jj = 1:n
%        tmp = vanGenuchten([testVals(ii,1),testVals(jj,2),testVals(500,3),testVals(500,4),testVals(500,5)]);
       tmp = vanGenuchten([pstar(1:2)',testVals(ii,3),testVals(jj,4),pstar(5)']);

       X3(jj,ii) = sum(((tmp-theta)./s).^2);
    end
end

% Loop over thetas and thetar
X4 = zeros(n);
for ii = 1:n
    for jj = 1:n
%        tmp = vanGenuchten([testVals(ii,1),testVals(jj,2),testVals(500,3),testVals(500,4),testVals(500,5)]);
       tmp = vanGenuchten([pstar(1:3)',testVals(ii,4),testVals(jj,5)]);

       X4(jj,ii) = sum(((tmp-theta)./s).^2);
    end
end
c = 35;
% k = linspace(10,500,50);
k = linspace(10,500,10);

figure();
subplot(2,2,1)
% imagesc(testVals(:,1),testVals(:,2),X1);colormap(bone)
contour(testVals(:,1),testVals(:,2),X1,k);colormap(bone)
hold on;
contour(testVals(:,1),testVals(:,2),X1,c);colormap(bone)


xlabel('\alpha')
ylabel('n')
title('\chi^2 Objective Function at Optimal Values m, \theta_S, \theta_R')
set(gca,'fontweight','bold')


subplot(2,2,2)
% imagesc(testVals(:,2),testVals(:,3),X2);colormap(bone)
contour(testVals(:,2),testVals(:,3),X2,k);colormap(bone)
hold on;
contour(testVals(:,2),testVals(:,3),X2,c);colormap(bone)


xlabel('n')
ylabel('m')
title('\chi^2 Objective Function at Optimal Values \alpha, \theta_S, \theta_R')
set(gca,'fontweight','bold')


subplot(2,2,3)
% imagesc(testVals(:,3),testVals(:,4),X3);colormap(bone)
contour(testVals(:,3),testVals(:,4),X3,k);colormap(bone);
hold on;
contour(testVals(:,3),testVals(:,4),X3,c);colormap(bone)


xlabel('m')
ylabel('\theta_S')
title('\chi^2 Objective Function at Optimal Values \alpha, n, \theta_R')
set(gca,'fontweight','bold')


subplot(2,2,4)
% imagesc(testVals(:,4),testVals(:,5),X4);colormap(bone)
contour(testVals(:,4),testVals(:,5),X4,k);colormap(bone)
hold on;
contour(testVals(:,4),testVals(:,5),X4,c);colormap(bone)


xlabel('\theta_S')
ylabel('\theta_R')
title('\chi^2 Objective Function at Optimal Values \alpha, n, m')
set(gca,'fontweight','bold')

% Summed Objective
X5 = X1+X2+X3+X4;
figure();
imagesc(X5);colormap(bone)