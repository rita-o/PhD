%%%%%%% EQUATION 1 & 2 - Solve for several g's at the same time %%%%%%% 
close all
npoints_list = [10, 100, 500, 800];
figure(1)
figure(2)
figure(3)

for i=1:length(npoints_list)
    
%% Initiate (run only once to use the same G's all the time)
alpha   = 0.1;
beta    = 0.7;
npoints = npoints_list(i);
G_raw   = 0.7*ones(npoints,1);
V_raw   = 7*ones(npoints,1);
G       = awgn(G_raw,30); % Make gaussian noise on g's
%G       = normrnd(0.7,0.02,[npoints,1]);
v       = 7;

% Plot noise on g
figure(1)
subplot(2,1,1)
plot(1:npoints,G_raw)
hold on
plot(1:npoints,G)
title('G-ratio')
xlabel('Number of samples')
ylabel('G-ratio value')
legend('Clean signal','Noise signal')
subplot(2,1,2)
histogram(G,15)
ylabel('Number of samples')
xlabel('G-ratio value')
hold off

%% Solve
x0      = [0.3;0.1]; % First guess on M and theta
options = optimoptions('fsolve','Display','iter','FunctionTolerance',10^-6,'OptimalityTolerance',10^-6,'algorithm','Levenberg-Marquardt');
res1    = fsolve(@(x) myfun(x,alpha,beta,G,v),x0,options);
options = optimoptions('lsqnonlin','Display','iter','FunctionTolerance',10^-6,'OptimalityTolerance',10^-6,'algorithm','Levenberg-Marquardt');
res2    = lsqnonlin(@(x) myfun(x,alpha,beta,G,v),x0,[],[],options); % other option

M_with_noise        = res1(1);
Theta_with_noise    = res1(2);

sprintf('fsolver - Mode is: %.003f and Theta is: %.003f',M_with_noise,Theta_with_noise)

M_with_noise        = res2(1);
Theta_with_noise    = res2(2);

sprintf('lsqnonlin - Mode is: %.003f and Theta is: %.003f',M_with_noise,Theta_with_noise)

%% Get True Value
g   = 0.7;
v   = 6;
syms M theta
[~, ~,v_eqn_GammaLevel] = make_equations_v(M,theta,v,alpha,beta,1);
[~, ~,g_eqn_GammaLevel] = make_equations_g(M,theta,g,alpha,beta,1);

% Solve equations
[s1, s2]    = vpasolve([g_eqn_GammaLevel v_eqn_GammaLevel],[M theta],[0 1; 0 1]);

M_true              = s1;
Theta_true          = s2;

sprintf('Mode is: %.003f and Theta is: %.003f',M_true,Theta_true)

%% Calculate difference between true and noise values
diff_M          = M_true - M_with_noise;
diff_Theta      = Theta_true - Theta_with_noise;

sprintf('Difference in mode is: %.003f and in theta is: %.003f',diff_M,diff_Theta)

figure(2)
subplot(2,1,1)
scatter(npoints_list(i),abs(diff_M),'Filled')
xlim([0 800])
xlabel('N of points')
ylabel('Error')
title('Absolute diff in M')
hold on
subplot(2,1,2)
scatter(npoints_list(i),abs(diff_Theta),'Filled')
hold off
xlim([0 800])
xlabel('N of points')
ylabel('Error')
title('Absolute diff in Theta')
hold on


%% Display 

theta       = Theta_true;
mode        = M_true;
shape       = mode/theta+1;
X           = 0:0.01:4;
pdf_true    = gampdf(X,shape,theta);


theta       = Theta_with_noise;
mode        = M_with_noise;
shape       = mode/theta+1;
X           = 0:0.01:4;
pdf_noise   = gampdf(X,shape,theta);

figure(3)
subplot(1,4,i)
t   =  plot(X,pdf_true,'r');
hold on
n   = plot(X,pdf_noise,'b');
hold on
title(sprintf('%d points',npoints_list(i)))
xlabel('Axon radius (um)')
ylabel('pdf value')
legend([t, n],'True value','Noise value')
hold on


end
%% FUNCTION MATLAB - All together

function F = myfun(x,alpha,beta,G,v)
% Variables
M = x(1);
theta = x(2);

% Velocity equation
F = [-v + ((2*5.5)/beta)*theta^(1-alpha) ...
    * gamma(M/theta+1-alpha)/gamma(M/theta+1) ... 
    * (M/theta+1-alpha)]; % Velocity function

% G ratio equation
for i=1:length(G)
    g = G(i);
    F  = [F;
         -g^2 + 2^(2*alpha)*beta^2*theta^(2*alpha) ...
        * gamma(M/theta+1)/gamma(M/theta+1-2*alpha) ... % Without approximmation
        * (M + 3*theta + 2*(theta^2)/M) ...
        / (M + (3 - 4*alpha)*theta + (2+4*alpha^2-6*alpha)*(theta^2)/M)];
end

end

%% FUNCTION MATLAB - Equation 1 (g ratio)

function [g_eqn_EasyLevel_raw, g_eqn_HardLevel,g_eqn_GammaLevel] = make_equations_g(M,theta,g,alpha,beta,A)
g_eqn_EasyLevel_raw     = g^2 == 2^(2*alpha)*beta^2*theta^(2*alpha) ...
    * A*(M/theta)^(2*alpha) ... % Lower limit and factor A
    * (M + 3*theta + 2*(theta^2)/M) ...
    / (M + (3 - 4*alpha)*theta + (2+4*alpha^2-6*alpha)*(theta^2)/M);
g_eqn_HardLevel         = g^2 == 2^(2*alpha)*beta^2*theta^(2*alpha) ...
    * (M/theta+1/2-alpha)^(2*alpha) ... % Complex limit without factor
    * (M + 3*theta + 2*(theta^2)/M) ...
    / (M + (3 - 4*alpha)*theta + (2+4*alpha^2-6*alpha)*(theta^2)/M);
g_eqn_GammaLevel        = g^2 == 2^(2*alpha)*beta^2*theta^(2*alpha) ...
    * gamma(M/theta+1)/gamma(M/theta+1-2*alpha) ... % Without approximmation
    * (M + 3*theta + 2*(theta^2)/M) ...
    / (M + (3 - 4*alpha)*theta + (2+4*alpha^2-6*alpha)*(theta^2)/M);
end

%%  FUNCTION MATLAB - Equation 2 (velocity)

function [v_eqn_EasyLevel, v_eqn_HardLevel,v_eqn_GammaLevel] = make_equations_v(M,theta,v,alpha,beta,C)
v_eqn_EasyLevel     = v == ((2*5.5)/beta)*theta^(1-alpha) ...
    * C*(M/theta)^(-alpha) ... % Lower limit and factor A
    * (M/theta+1-alpha);
v_eqn_HardLevel         = v == ((2*5.5)/beta)*theta^(1-alpha) ...
    * (M/theta+1/2-alpha/2)^(-alpha) ... % Complex limit without factor
    * (M/theta+1-alpha);
v_eqn_GammaLevel        = v == ((2*5.5)/beta)*theta^(1-alpha) ...
    * gamma(M/theta+1-alpha)/gamma(M/theta+1) ... % Without approximmation
    * (M/theta+1-alpha);
end
