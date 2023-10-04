% ==========================
% Problem set solving the stochastic neoclassical growth model different
% ways
% ==========================
% clear the workspace
clear
% close all figures
close all
addpath('C:\Users\tbroer\Dropbox\Teaching\PSE\2021 Quantitative Macro\Problem sets\PS RBC')
addpath('C:\Users\tbroer\Dropbox\Teaching\PSE\2021 Quantitative Macro\Problem sets')

% ============
% parameters
% ============
theta=0.4 % capital share
beta = 0.99; % discount factor
rho = 0.95;   % persistence of TFP shock
gamma = 1.00000001; % CRRA coefficient (for 1 equals log, but need to replace the function, so set close to 1)
delta=0.1;

% ============
% options and convergence criteria
% ============
T=150; % periods for transition

%mean of capital non-stochastic steady state
kbar=[FILL THIS IN];
% initial level of capital in the transition
k_0=kbar*0.9;


% =====================
% 4. Log-linearization
% =====================

% some ratios as function of parameters
ybar=[FILL THIS IN];
cbar=[FILL THIS IN];
ckrat=[FILL THIS IN];
R=1/beta;
margprod=R-1+delta;

% a. write system as A E[y_t+1]+B y_t=0

% Euler equation
% -gamma*chat_t=(-gamma)*chat_t+1
% +beta*margprod*khat_t+1+beta*margprod*zbar*zhat_t+1

% becomes
% (-gamma)*chat_t-beta*margprod*khat_t+1  -  beta*margprod*zhat_t+1+gamma chat_t=0

% res constraint
% khat_t+1 - ckrat*chat_t - (margprod+(1-delta)) khat_t -kbar^theta*zbar*zhat_t=0

% order c,k,z
A=[[FILL THIS IN]];

B= [[FILL THIS IN]];

D = [FILL THIS IN];

% note that these are right-hand eigenvectors
[ [FILL THIS IN]]=eig(D);
aaa=[FILL THIS IN];

BKcond=[FILL THIS IN];

if BKcond~=1
    disp('BK conditions not satisfied')
else
    indic=find(abs(diag(lambda))>1);
    indic1=find(abs(diag(lambda))<=1);
    polfunc_temp=aaa(indic,:);
    % this gives the consumption policy function
    polfunc=[FILL THIS IN]
end


% calculate the deterministic transition
k_lin(1)=k_0;
for t=2:T
    c_lin(t-1)=[FILL THIS IN];
    k_lin(t)=[FILL THIS IN];
end


% ==============
% Figures
% ==============

% plot the transition
figure(2)
title('Simulated transition - deterministic')
subplot(2,1,1)
plot(k_lin,'r','Linewidth',2)
subplot(2,1,2)
plot(c_lin,'r','Linewidth',2)
set(h,'fontsize',12,'Interpreter','Latex');

