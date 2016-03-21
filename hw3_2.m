% Homework 3
% Lingxiao Zhao

% Q2 A model of separation
clc
clear
% Propose parameter values for structural parameters
beta = 1;
mu = 0;
sigma = 1;
gamma = 0.5;
phi = 15;
H = 1000; % number of simulation draws
t_min = 21;
t_max = 60; 
u_z = 0;
sigma_z = 0.5;
z_l = -5;

% Simulation
t_min = 21;
t_max = 60;
T = unidrnd(40,H,1) + 20;
e = normrnd(mu,sigma,H,1);
w = beta * T + e;
r = gamma * T + phi;
z = ones(H,1);
    for i = 1:H
        z(i)=normrnd(0,0.5*T(i));
    end
for i=1:H;
        if w(i)>=r(i)
        v(i)=w(i);
        a(i)=i;
        else
        v(i)=0;
        a(i)=0;
        end 
end
indices_T=find(T<=40);
whos indices_T;
indices_z=find(z(indices_T)<=z_l);
whos indices_z;
lm_3=fitlm(T(indices_z),w(indices_z),'linear')
beta_hat_3=lm_3.Coefficients.Estimate(2)
mu_hat_3=lm_3.Coefficients.Estimate(1)
sigma_hat_3=sqrt(var(lm_3.Residuals.Raw))

for j=1:40
    count_number(j)=0;
    for i=1:H
            if T(i)== j+20
            count_number(j)=count_number(j)+1;
            end
    end
end

for j=1:40
    count(j)=0;
    for i=1:H
        if T(i)==j+20
            if w(i) >= r(i)
                count(j)=count(j)+1;
            else
                count(j)=count(j);
            end
        end 
    end
    prob(j)=count(j)/count_number(j);
end