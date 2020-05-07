% make sure we have a clean environment
clear
rand('state',0);
randn('state',0)
addpath 'C:\Users\snowfield\Desktop\Backup\Math\Inverse Theory\PEIP-code\Lib';

% Generate the x and y valuesy
t = [3.493500000000000;4.285300000000000;5.137400000000000;5.818100000000000;6.863200000000000;8.184100000000000];
x = [6;10.133300000000000;14.266700000000000;18.400000000000000;22.533300000000000;26.666700000000000];
ytrue = t;
N=length(x);
sig = 0.1;
sigma=sig*ones(N,1);
y = ytrue+sig*randn(size(x)).*ytrue;

% Build the matrix
G=[ones(size(x)) x];

% Apply the weighting
yw = y./sigma;
Gw = G./[sigma,sigma];

% Solve for the least-squares solution
disp('least-squares solution')
m2 = Gw\yw

% Solve for the 1-norm solution
disp('1-norm solution')
m1 = irls(Gw,yw,1.0e-5,1.0e-5,1,125)

% comute modeled
ymod1 = G * m1;
ymod2 = G * m2;

%% Plots
%  1. Data and fitted line.
%  2. residuals.
figure(1)
clf
plot(x,y,'ko');
hold on
p = plot(x,ymod1,'k',x,ymod2,'b--');
p(1).LineWidth = 2;
%%p(2).Marker = '*';
%%plot();
xlabel('x');
ylabel('t');
title('$\ell 1$ and $\ell 2$ Regression','interpreter','latex');
legend({'Data','$\ell 1$','$\ell 2$'},'Interpreter','latex','location','northwest')
set(gca,'fontweight','bold')

disp('Displaying Data and Linear Regression Line (fig 1)')


figure(2)
clf
plot(x,y-ymod1,'ko');hold on
plot(x,y-ymod2,'ro');
xlabel('x');
ylabel('r');
title('$\mathbf{r}_{\ell 1}$ and $\mathbf{r}_{\ell 2}$ Residuals','interpreter','latex')
legend({'$\mathbf{r}_{\ell 1}$','$\mathbf{r}_{\ell 2}$'},'Interpreter','latex','location','northwest')
set(gca,'fontweight','bold')



disp('Displaying Model Residuals vs. x (fig 2)');
%%
% Get the covariance matrix
ginv = inv(Gw'*Gw)*Gw';

disp('Covariance matrix')
covm = ginv*ginv'

% Find the parameter correlations
s=sqrt(diag(covm))

disp('correlation matrix')
r = covm./(s*s')

% Output covm and the eigenvalues/eigenvectors of covm.
disp('Covariance matrix for fitted parameters.')
covm

disp('Eigenvalues/eigenvectors of the covariance matrix');
[u,lam]=eig(inv(covm))
disp('95% confidence ellipsoid semiaxis lengths');
semi_axes = [sqrt(chi2inv(0.95,2)*(1./diag(lam)))]'
disp('95% confidence ellipsoid semiaxes')

[semi_axes(1)*u(:,1), semi_axes(2)*u(:,2)]

%
% Plot the 95% error ellipses for each pair of parameters
% Note that because we're doing pairs of parameters there are 2
% degrees of freedom in the Chi-square here, rather than 3.  
%
%generate a vector of angles from 0 to 2*pi
theta=(0:.01:2*pi)';
delta=sqrt(chi2inv(0.95,2));
%the radii in each direction from the center
r=zeros(length(theta),2);

figure(3)
clf

% compute the data for the m1, m2 ellipsoid.
C=covm((1:2),(1:2));
[u,lam]=eig(inv(C));
%calculate the x component of the ellipsoid for all angles
r(:,1)=(delta/sqrt(lam(1,1)))*u(1,1)*cos(theta)+(delta/sqrt(lam(2,2)))*u(1,2)*sin(theta);
%calculate the y component of the ellipsoid for all angles
r(:,2)=(delta/sqrt(lam(1,1)))*u(2,1)*cos(theta)+(delta/sqrt(lam(2,2)))*u(2,2)*sin(theta);

% plot the data for the m1, m2 ellipsoid
plot(m2(1)+r(:,1),m2(2)+r(:,2),'k');
fill(m2(1)+r(:,1),m2(2)+r(:,2),'k');

xlabel('m_1 (s)');
ylabel('m_2 (km)');
title('$\mathbf{m_{L2}}$ 95\% Confidence Region','interpreter','latex');
dim = [.2 .05 .3 .3];
str = {['$\mathbf{m_{1}} =$ ',num2str(m2(1)), ' $\pm$ ', num2str(semi_axes(1))],['$\mathbf{m_{2}} =$ ',num2str(m2(2)), ' $\pm$ ', num2str(semi_axes(2))]};
annotation('textbox',dim,'String',str,'FitBoxToText','on','interpreter','latex');
set(gca,'fontweight','bold')
%%
% Because there are 2 parameters to estimate, we have N-2 degrees
% of freedom.
%
dof = N-2;
disp(['Chi-square misfit for ',num2str(dof),' dof'])
chi2 = norm((y - G*m2)./sigma)^2

% Find the p-value for this data set
disp('chi-square p-value')
p = 1-chi2cdf(chi2,dof)

% Monte Carlo Section
y0 = G*m2; 

for nreal = 1:1000
  % Generate a trial data set of perturbed, weighted data
  ytrial = y0+sigma.*randn(N,1);
  ywtrial=ytrial./sigma;
  mmc(nreal,:)=(Gw\ywtrial)';
  chimc(nreal)= norm((G*mmc(nreal,:)'-ytrial)./sigma)^2;
end
% Calculate the Theoretical X^2 pdf
chi = chi2pdf(linspace(min(chimc),max(chimc),1000),dof);
% chi = chi2pdf(mmc,dof);


% Plot the histogram of chi squared values
figure(4)
clf
histogram(chimc,30,'FaceColor','k','normalization','pdf');
hold on;
plot(linspace(min(chimc),max(chimc),1000),chi,'k','linewidth',2)
ylabel('Probability Density');
xlabel('\chi^2')
title('Simulated and Theoretical \chi^2 Distribution for L2 Solution');
legend('Experimental PDF','Theoretical PDF')
set(gca,'fontweight','bold')

%%
% Monte Carlo Section for the 1-norm model
% The 1-norm estimated data
y01 = G*m1; 

nreal=1000;

disp('generating Monte Carlo realizations');
%generate the nreal Monte Carlo solutions
for j = 1:nreal
  % Generate a trial data set of perturbed, weighted data
  ytrial = y01+sigma.*randn(N,1);
  ywtrial=ytrial./sigma;
  
  % Store the 1-norm parameters and misfit
  mmc(j,:)=irls(Gw,ywtrial,1.0e-5,1.0e-5,1,500)';
  chimc(j)= norm((G*mmc(j,:)'-y01)./sigma)^2;
end

% figure out how many Monte Carlo points are in the %95 confidence region
disp('confidence interval inclusion check (should be about 0.95):') 
nnz(chimc<=chi2inv(.95, dof))/nreal

% Calculate the covariance of the parameters in the realizations
disp('Emperical covariance of m1 models');
covmemp=mmc-ones(nreal,1)*mean(mmc);
covmemp=(covmemp'*covmemp)/nreal

% Get the 1.96-sigma (95%) conf intervals
disp('95% parameter confidence intervals (m-, mest, m+) on 1-norm solution')
del = 1.96*sqrt(diag(covmemp));
[m1-del , m1 , m1+del]
figure(1);
dim = [.5 .05 .3 .3];
str = {['$\mathbf{m_{L1}}(1) =$ ',num2str(m1(1)), ' $\pm$ ', num2str(del(1))],['$\mathbf{m_{L1}}(2) =$ ',num2str(m1(2)), ' $\pm$ ', num2str(del(2))]};
annotation('textbox',dim,'String',str,'FitBoxToText','on','interpreter','latex');
%% contributions from each of the data points to the 1-norm misfit measure
 
for i = 1:N
    xnew = x;
    tnew = t;
    xnew(i) = [];
    tnew(i) = [];
    ytrueNew = tnew;
    Nnew=length(xnew);
    sigmanew=sig*ones(Nnew,1);
    ynew = ytrueNew+sig*randn(size(xnew)).*ytrueNew;

    % Build the matrix
    Gnew=[ones(size(xnew)) xnew];

    % Apply the weighting
    ywnew = ynew./sigmanew;
    Gwnew = Gnew./[sigmanew,sigmanew];

    % Solve for the least-squares solution
    disp('least-squares solution')
    m2new(i,:) = Gwnew\ywnew

    % Solve for the 1-norm solution
    disp('1-norm solution')
    m1new(i,:) = irls(Gwnew,ywnew,1.0e-5,1.0e-5,1,1)

    % comute modeled
    ymod1new(:,i) = Gnew * m1new(i,:)';
    ymod2new(:,i) = Gnew * m2new(i,:)';
    
    % compute residuals
    r1(:,i) = ynew - ymod1new(:,i);
    r2(:,i) = ynew - ymod2new(:,i);
    
    % compute norm1 of residuals
    nr1(1,i) = norm(r1(:,i),1);
    nr2(1,i) = norm(r2(:,i),1);

end

figure(5)
clf
plot(1:N,nr1,'ko');
hold on
plot(1:N,nr2,'ro');
xlabel('Data points');
ylabel('norm1 of Residuals');
title('$\mathbf{r}_{\ell 1}$ and $\mathbf{r}_{\ell 2}$ Residuals','interpreter','latex')
legend({'$\mathbf{r}_{\ell 1}$','$\mathbf{r}_{\ell 2}$'},'Interpreter','latex','location','northwest')
set(gca,'fontweight','bold')


%%%%%%%%% confidence interval and data outliers %%%%%%%%%%%%%%%%%%%

con1 = [mean(nr1) - (tinv(0.975,5)*std(nr1)/(6^0.5)), mean(nr1) + (tinv(0.975,5)*std(nr1)/(6^0.5))];
con2 = [mean(nr2) - (tinv(0.975,5)*std(nr2)/(6^0.5)), mean(nr2) + (tinv(0.975,5)*std(nr2)/(6^0.5))];    
for i = 1:6
    if (nr1(1,i) > con1(1,2) || nr1(1,i) < con1(1,1))
        disp(['Data rank number ', num2str(i), ' acording to norm1 of mL1 estimate is an outlier'])
    end
end
for i = 1:6
    if (nr2(1,i) > con2(1,2) || nr2(1,i) < con2(1,1))
        disp(['Data rank number ', num2str(i), ' acording to norm1 of mL2 estimate is an outlier'])
    end
end
