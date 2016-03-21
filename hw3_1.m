% Homework 3
% Lingxiao Zhao

% Q1 A truncated observation of wages
clc
clear

% Propose parameter values
beta = 1;
mu = 0;
sigma = 1;
gamma = 0.5;
phi = 15;
H = 1000; % number of simulation draws

% Simulation
t_min = 21;
t_max = 60;
T = unidrnd(40,H,1) + 20;
e = normrnd(mu,sigma,H,1);
w = beta * T + e;
r = gamma * T + phi;
for i = 1:H;
    if w(i) >= r(i)
        v(i) = w(i);
        a(i) = 1;
    else 
        v(i) = 0;
        a(i) = 0;
    end
end

% Estimate beta, mu, sigma without accounting fot the selection
lm_1 = fitlm(T,w,'linear');
anova(lm_1,'summary');
beta_hat_1 = lm_1.Coefficients.Estimate(2)
mu_hat_1 = lm_1.Coefficients.Estimate(1)
sigma_hat_1 = sqrt(var(lm_1.Residuals.Raw))

% Estimate beta, mu, sigma with the selection
lm_2=fitlm(T(a~=0),w(a~=0),'linear');
anova(lm_2,'summary');
beta_hat_2=lm_2.Coefficients.Estimate(2)
mu_hat_2=lm_2.Coefficients.Estimate(1)
sigma_hat_2=sqrt(var(lm_2.Residuals.Raw))

