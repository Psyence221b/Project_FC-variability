% Robust linear regression.
%% Usage 1: [rob] = andlab_robustfit(x,y)
% 
% [rob] is the output of robustfit in struct format, including r-value calculated with re-weghted x and y:
% y_weighted = sqrt(rob.stats.w).*(y-rob.b(1)); % reconstruct the weighted x and y series that robustfit used
% x_weighted = sqrt(rob.stats.w).*x;
% rob.corrw = corr(x_weighted,y_weighted); % robust fit regression output (r-value)
% 
% the weigted data is in the structre, rob.weighted.x & rob.weighted.y
%
% To know other fields of the structure, please see 'help robustfit'
% http://www.mathworks.com/support/solutions/en/data/1-CMABGO/index.html?product=ML&solution=1-CMABGO
%
%% Usage 2: [rob, leastsq] = my_robustfit(x,y)
%
% [leastsq] is the output of least square regression % regular correlation output
%
% To know other fields of the structure, please see 'help regress'
%
%% Description example
% "We employed iterative reweighted least squares (the robustfit function from Matlab, Mathworks,
% Natick, MA, USA), given that standard Pearson correlation is very sensitive to even a few influential data points 
% (Wager et al.,2005; see also Wilcox, 2005)" - Choi, Padmala & Pessoa (2011, neuroimage)
%
% Modified by TH Lee (Dec, 2015) based on Moon (Jan, 2013)' robust fit
% Updated by THLee (Nov 2019)
% function
%
% All rights are Not reserved.
%==========================================================================

function [rob, leastsq, rob_corrw] = dsn_robustfit(x,y)
% simple regression
[leastsq_b,leastsq_bint,leastsq_r,leastsq_rint,leastsq_stats] = regress(y,[ones(size(x)) x]); 

[rob_b, rob_stats] = robustfit(x,y);

%  scatter(x,y,'filled'); grid on; hold on
%  plot(x,leastsq_b(1)+leastsq_b(2)*x,'r','LineWidth',2);
%  plot(x,rob_b(1)+rob_b(2)*x,'g','LineWidth',2)
%  legend('Data','Ordinary Least Squares','Robust Regression')

% reconstruct the weighted x and y series that robustfit used
yw=sqrt(rob_stats.w).*(y-rob_b(1));                 
xw=sqrt(rob_stats.w).*x;

% calculate correlation coefficient
rob_corrw = corr(xw,yw);                     

%% Output
leastsq.corrcoef = leastsq_b(2)/abs(leastsq_b(2)) * sqrt(leastsq_stats(1));
leastsq.b = leastsq_b;
leastsq.bint = leastsq_bint;
leastsq.r = leastsq_r;
leastsq.rint = leastsq_rint;
leastsq.stats = leastsq_stats;

rob.b = rob_b;
rob.corrw = rob_corrw;
rob.stats = rob_stats;
rob.weighted.x = xw;
rob.weighted.y = yw;


%% Plot regression line
%  scatter(x,y,'filled'); grid on; hold on
%  plot(x,leastsq_b(1)+leastsq_b(2)*x,'r','LineWidth',2);
%  plot(x,rob_b(1)+rob_b(2)*x,'g','LineWidth',2)
%  legend('Data','Ordinary Least Squares','Robust Regression')

