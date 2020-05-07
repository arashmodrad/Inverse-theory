% Example 11.4
% from Parameter Estimation and Inverse Problems, 3rd edition, 2018
% by R. Aster, B. Borchers, C. Thurber
%
% make sure we have a clean environment
clear
rand('state',0);
randn('state',0);

%Global variables for use by the mcmc function calls
global step;

addpath('C:\Users\snowfield\Desktop\Backup\Math\Inverse Theory\PEIP-code\Lib');
global theta;
global h;
global sigma;
[h] = xlsread('soil.xlsx','sheet1','B2:B29');
m = length(h);
theta = xlsread('soil.xlsx','sheet1','C2:C29');
dev = xlsread('soilstd.xlsx','sheet2','F2:F29');
dev = [NaN;dev];
sigma = inpaint_nans(dev,1);

disp('Note: This example will take several minutes to run on most computers')

% Generate the data set.
% mtrue=[.033,1.49,.33,0.86,.04]';
% mtrue=[.033,1.09,.33,0.86,.04]';
mtrue = [0.0080,1.01,.4553,.5434,.1503]';
thetaTrue=fun(mtrue);
% % sigma=0.01*ones(size(ytrue));
% y=thetaTrue+sigma.*randn(size(thetaTrue));

%set the MCMC parameters
%number of skips to reduce autocorrelation of models
skip=25000;
%burn-in steps
BURNIN=100000;
%number of posterior distribution samples
N=500000;
%MVN step size
step = 0.0035*ones(5,1);

% We assume flat priors here
%% initialize model at a random point on [0,1]
m0=[rand, rand + 1, rand,rand,rand]';

[mout,mMAP,pacc]=mcmc('lognormalprior','loglikelihood','generate','logproposal',m0,N);
 
%
% plot autocorrelations.
%

figure(1)
clf
%plot parameter correlations
laglen=25000;
% lags=(-laglen:laglen)';
% lags = [1,[2:100:laglen]]';
% acorr=zeros(2*laglen+1,5);

disp('Calculating autocorrelations for fig. 1')
for i=1:5
      [acorr(:,i),lags]=calc_corr(mout(i,:)',laglen);
end

for i=1:5
  subplot(5,1,i);
  bookfonts;
  plot([0 laglen],[0 0],'Color',[0.7 0.7 0.7],'LineWidth',3);
  hold on
%   plot(lags(laglen+1:200:laglen*2+1),acorr(laglen+1:200:laglen*2+1,i),'ko');
    plot(lags(round(length(acorr)./2) + 1:end),acorr(round(length(acorr)./2)+1:end,i),'ko');

  hold off
  ylabel(['A ( m_',num2str(i),')'])
  ylim([-0.5 1])
  if i~=5
    set(gca,'Xticklabel',[]);
  end
  bookfonts
  xlim([0,laglen]);
end
xlabel('Lag')
bookfonts

print -depsc2 c11MCMCmcorrbefore.eps
disp(['Displaying the autocorrelation of individual parameters before' ...
      ' thinning (fig. 1)']);

%  %track accepted new models to report acceptance rate
%  nacc=0;
%   
%   %sample the posterior here
%   for t = 2:N,
%       %calculate a candidate model
%    c=m(:,t-1)+randn(4,1).*step;
%       %calculate acceptance probability for c,m(:,t-1)
%       lnprobacc=getlnprobacc(m(:,t-1),c,x,y,sigma);
%       if log(rand) < lnprobacc
%           m(:,t)=c;
%           nacc=nacc+1;
%       else
%           m(:,t)=m(:,t-1);
%       end
%   end
disp(['Acceptance Rate: ',num2str(pacc)]);
  
%downsample results to reduce correlation
k=(BURNIN:skip:N);
 
mskip=mout(:,k);
  
%histogram results, and find the modes of the subdistributions as an
%estimte of the MAP model
disp(['m_map','  m_true'])
mMapMCMC = mMAP;
for ii = 1:5
%     hist(mout(ii,:), 100000);
[counts, center] = hist(mout(ii,:), 10000);

[~,ix] = max(counts);
mMAP(ii,1) = center(ix);%find(counts>0, 1, 'last');
end
% mMAP = mean(mout,2);
[mMAP,mtrue]
 
% estimate the 95% credible intervals
for i=1:5
%   msort=sort(mskip(i,:));
  msort=sort(mout(i,BURNIN:N));
    m2_5(i) = msort(round(2.5/100*length(mout(i,BURNIN:N))));
  m97_5(i) =  msort(round(97.5/100*length(mout(i,BURNIN:N))));
%   m2_5(i) = msort(round(2.5/100*length(mskip)));
%   m97_5(i) =  msort(round(97.5/100*length(mskip)));
  disp(['95% confidence interval for m', num2str(i),' is [', num2str(m2_5(i)),',', num2str(m97_5(i)),']'])
end

%plot a scatter plot and histogram of the posterior distribution
mlims=[0 0.025; 1 2 ; 0 1 ; 0.45 0.625; 0 .25];

figure(2)
clf
for i=1:5
  for j=1:5
    subplot(5,5,5*(i-1)+j)
    if i==j
      hist(mskip(i,:));
      h = findobj(gca,'Type','patch');
      set(h,'FaceColor','k')
      set(gca,'Yticklabel',[]);
      xlim(mlims(i,:));
%       bookfonts
    else
      plot(mskip(j,:),mskip(i,:),'k.','Markersize',6,'Color',[0.6 0.6 0.6]);
      hold on
      % plot the true answer as a large black dot
      plot(mtrue(j),mtrue(i),'k.','Markersize',24);
      % plot the accepted answers as gray dots
      plot(mMAP(j),mMAP(i),'b.','Markersize',24);
      % plot the 95% ci as a box
      plot([m2_5(j),m97_5(j)],[m2_5(i),m2_5(i)],'k-','LineWidth',1);
      plot([m2_5(j),m97_5(j)],[m97_5(i),m97_5(i)],'k-','LineWidth',1);
      plot([m2_5(j),m2_5(j)],[m2_5(i),m97_5(i)],'k-','LineWidth',1);
      plot([m97_5(j),m97_5(j)],[m2_5(i),m97_5(i)],'k-','LineWidth',1);
      xlim(mlims(j,:));
      ylim(mlims(i,:));
%       bookfonts
      hold off
    end
  end
end

subplot(5,5,1)
ylabel('m_1')
% bookfonts
subplot(5,5,6)
ylabel('m_2')
% bookfonts
subplot(5,5,11)
ylabel('m_3')
% bookfonts
subplot(5,5,16)
ylabel('m_4')
% bookfonts
subplot(5,5,21)
ylabel('m_5')
% bookfonts
subplot(5,5,21)
xlabel('m_1')
% bookfonts
subplot(5,5,22)
xlabel('m_2')
% bookfonts
subplot(5,5,23)
xlabel('m_3')
% bookfonts
subplot(5,5,24)
xlabel('m_4')
% bookfonts
subplot(5,5,25)
xlabel('m_5')
% bookfonts
print -depsc2 c11MCMCscat.eps
disp('Displaying the true parameters, thinned points, and 95% confidence intervals (fig. 2)');

%plot parameter sample histories
figure(3)
clf
for i=1:5
  subplot(5,1,i)
  plot([1 length(mskip)],[mtrue(i) mtrue(i)],'Color',[0.6 0.6 0.6],'LineWidth',3);
  hold on
  plot(mskip(i,:),'ko')
  hold off
  if i~=5
    set(gca,'Xticklabel',[]);
  end
  xlim([1 length(mskip)])
end
xlabel('Sample Number')
subplot(5,1,1)
ylabel('m_1')
% bookfonts
subplot(5,1,2)
ylabel('m_2')
% bookfonts
subplot(5,1,3)
ylabel('m_3')
% bookfonts
subplot(5,1,4)
ylabel('m_4')
% bookfonts
subplot(5,1,5)
ylabel('m_5')
% bookfonts
set(gca,'fontweight','bold')
print -depsc2 c11MCMCmhist.eps
disp('Displaying the true parameters and thinned sample histories (fig. 3)');

%plot parameter correlations
figure(4)
clf
laglen=round(length(mskip)./2);
lags=(-laglen:laglen)';
acorr=zeros(2*laglen+1,4);
for i=1:5
  acorr(:,i)=calc_corr_skip(mskip(i,:)',laglen);
  subplot(5,1,i);
  plot([0 laglen],[0 0],'Color',[0.7 0.7 0.7],'LineWidth',3);
  hold on
  plot(lags(laglen+1:laglen*2+1),acorr(laglen+1:laglen*2+1,i),'ko');
  hold off
  ylabel(['A ( m_',num2str(i),')'])
  ylim([-0.5 1])
%    bookfonts
  if i~=5
    set(gca,'Xticklabel',[]);
  end
end
xlabel('Lag')
% bookfonts
print -depsc2 c11MCMCmcorr.eps
disp('Displaying the autocorrelation of thinned parameter samples (fig. 4)');

% Q-Q Plot
figure(5);clf;
for ii = 1:5
    subplot(1,5,ii)
qqplot(mskip(ii,:))
title(['m_',num2str(ii)])
if ii == 1
    ylabel('Quantiles of Posterior Samples')
end
if ii == 3
    xlabel('Quantiles of Normal Distribution')
end
set(gca,'fontweight','bold')
end
