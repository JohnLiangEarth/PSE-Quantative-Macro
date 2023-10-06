% ==========================
% Question 2: Problem set solving the stochastic neoclassical growth model different
% ways
% ==========================
% clear the workspace
clear
% close all figures
close all
addpath('/Users/tubaseker/Desktop/QuantitativeMacro/PS2')

% ============
% parameters
% ============
alpha = 0.3; % capital share
beta = 0.99; % discount factor
rho = 0.95;   % persistence of TFP shock
sigma = 1.0001; % CRRA coefficient (for 1 equals log, but need to replace the function, so set close to 1)
delta=0.025;

% ============
% options and convergence criteria
% ============
T=150; % periods for transition

%mean of capital non-stochastic steady state
kbar = (alpha*beta/(1-beta*(1-delta)))^(1/(1-alpha)); %write feasibility constraint at the steady-state and EE. Then, replace c_bar from feasibiltiy constraint into the EE, and obtain k_bar. 
% initial level of capital in the transition
k_0=kbar*0.9;

%%
% =====================
% Log-linearization
% =====================

% some ratios as function of parameters
ybar = kbar^(alpha);
cbar = ybar-delta*kbar;
ckrat=cbar/kbar; %c/k ratio
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
A=[-sigma beta*alpha*(alpha-1)*kbar^(alpha)-1; 0 1];

B= [-sigma 0; -ckrat alpha*kbar^(alpha-1)+(1-delta)];

D = A\B; %inv(A)*B can be slower;

% note that these are right-hand eigenvectors
[ev, lambda]=eig(D); %Note: [V,D] = eig(A) produces a diagonal matrix D of eigenvalues and 
                            %a full matrix V whose columns are the corresponding eigenvectors  
                            %so that A*V = V*D
aaa=inv(ev);


BKcond=sum(abs(diag(lambda))>1);
%If all the eigenvalues are less than 1 in absolute value then the system is stable. 
% If all the eigenvalues are greater than 1 in absolute value then the system is unstable. 
% If at least one eigenvalue is less than 1 in absolute value the system is saddle-path stable.

if BKcond~=1
    disp('BK conditions not satisfied')
else
    indic=find(abs(diag(lambda))>1);
    indic1=find(abs(diag(lambda))<=1);
    polfunc_temp=aaa(indic,:);
    % this gives the consumption policy function
    polfunc = -polfunc_temp(1,2)/polfunc_temp(1,1); % policy function for chat(t+1)
    % Recall that we found policy function as -p12/p11*khat(t+1)
end


% calculate the deterministic transition
k_lin(1)=0.9*kbar; %k0
for t=2:T
    c_lin(t-1) = polfunc*((k_lin(t-1)-kbar)/kbar)*cbar+cbar;
    k_lin(t)=k_lin(t-1)^alpha+(1-delta)*k_lin(t-1)-c_lin(t-1);
end

%%
% ==============
% Figures
% ==============

subplot(2, 1, 1); % 2 rows, 1 column, first subplot
plot(k_lin,'Linewidth',2);
title('Simulated transition for capital stock','interpreter','latex','FontSize',15);
xlabel('Periods ($t$)','interpreter','latex','FontSize',16);
ylabel('Capital $k_t$','interpreter','latex','FontSize',16);
grid on;
specialPointK = kbar; % Replace with the y-coordinate of the special point
specialPointLabelK = 'k_*'; % The label for the special point
line([min(k_lin), 150], [specialPointK, specialPointK], 'Color', 'red', 'LineStyle', '--', 'LineWidth', 1);
text(0.5, specialPointK, specialPointLabelK, 'FontSize', 12, 'HorizontalAlignment', 'right', 'VerticalAlignment', 'middle');

% Create the second subplot (lower plot)
subplot(2, 1, 2); % 2 rows, 1 column, second subplot
plot(c_lin,'Color', 'black','LineWidth',2);
title('Simulated transition for consumption','interpreter','latex','FontSize',15);
xlabel('Periods ($t$)','interpreter','latex','FontSize',16);
ylabel('Consumption $c_t$','interpreter','latex','FontSize',16);
grid on;
specialPointC = cbar; % Replace with the y-coordinate of the special point
specialPointLabelC = 'c_*'; % The label for the special point
line([min(c_lin), 150], [specialPointC, specialPointC], 'Color', 'red', 'LineStyle', '--', 'LineWidth', 1);
text(0.5, specialPointC, specialPointLabelC, 'FontSize', 12, 'HorizontalAlignment', 'right', 'VerticalAlignment', 'middle');
grid on;
saveas(gcf,'loglinpath.png');
hold off; 

%%
% ==========================
% Question 3: Model with labor supply
% ==========================

lbar = 1/3; %Calibration target
kbar2 = kbar*lbar;
ybar2 = alpha*kbar2^(alpha)*lbar^(1-alpha);
cbar2 = ybar2-delta*kbar;
rbar = alpha*(ybar2/kbar2);
% Write system as A E[y_t+1]+B y_t=0

a = [kbar2/ybar2,0,0; 0,0,0; beta*(1-alpha)*rbar,-beta*rbar*(1-alpha), sigma];
b = [alpha, 1+alpha,-sigma; 0,0,sigma];


nk = 1; %nk: number of state variables ?2
[f,p] = solab(a,b,nk);

% Extract policies and capital law of motion
c_polfunc=f(1,:);
n_polfunc=f(2,:);
lom=p(1,:);

% Grids

kmax=1.1*kbar2;
kmin=kbar2*0.9;
if delta==1
    if linear==1
        kgrid=linspace(kmin,kmax,N);
    else
        temp=linspace(0,0.5,N).^5/0.5^5*(kmax-kmin);
        kgrid=kmin+temp;
    end
else
    kgrid=linspace(kmin ,kmax,N);
end

%Multiply each with xbar since policy functions are in log deviations:
for j = 1:M
    c_pol_lin(:,j) = cbar2*exp(c_polfunc(1)*log(kgrid/kbar2));
    n_pol_lin(:,j) = lbar*exp(n_polfunc(1)*log(kgrid/kbar2));
    k_pol_lin(:,j)=kbar2*exp(lom(1)*log(kgrid/kbar2));
end


