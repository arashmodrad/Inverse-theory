% Example 11.4
% from Parameter Estimation and Inverse Problems, 3rd edition, 2018
% by R. Aster, B. Borchers, C. Thurber
% lp=logprior(m)
%
% For this problem, our prior is uniform on m1=[0 2], m2=[-1 0], m3=[0 2],
% m4=[-1 0]
%
function lp=lognormalprior(m)
% mtrue = [.033,1.49,.33,0.86,.04];
mtrue = [0.0080,1.01,.4553,.5434,.1503]';
mstd = 0.1.*ones(5,1);
if (m(1)>mtrue(1) - 5.*mstd(1)) && (m(1)<mtrue(1)+5.*mstd(1)) && (m(2)>mtrue(2) - 5.*mstd(2))...
        && (m(2)<mtrue(2)+5.*mstd(2)) && (m(3)>mtrue(3) - 5.*mstd(3)) && (m(3)<mtrue(3)+5.*mstd(3))...
        && (m(4)>mtrue(4) - 5.*mstd(4)) && (m(4)<mtrue(4)+5.*mstd(4)) && (m(5)>mtrue(5) - 5.*mstd(5)) && (m(5)<mtrue(5)+5.*mstd(5))
  lp=1./(mstd.*m.*sqrt(2.*pi)).*exp(-(log(m) - mtrue).^2./(2.*mstd.^2));
else
  lp=-Inf;
end
