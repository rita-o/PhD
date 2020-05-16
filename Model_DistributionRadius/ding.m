% Find alpha and beta

d = 0.1:0.1:3;

figure()
alpha   = 0.2;
beta    = 0.8/((0.69)^alpha);
g = beta*(d).^alpha;
plot(d,g,'r')
hold on
alpha   = 0.1;
beta    = 0.8/((0.69)^alpha);
g = beta*(d).^alpha;
plot(d,g,'b')
alpha   = 0.15;
beta    = 0.8/((0.69)^alpha);
g = beta*(d).^alpha;
plot(d,g,'g')

xlabel('diameter (um)')
ylabel('gratio')

alpha=0.1;
beta=0.7;

%% Find Constant A
mode=0.1;
R = [];
figure()
for theta=0.02:0.01:0.7
    LL      = (mode/theta)^(2*alpha);
    UL      = (mode/theta + 1)^(2*alpha);
    TV      = gamma(mode/theta+1)/gamma(mode/theta+1-2*alpha);
    GV = (mode/theta+0.5-alpha)^(2*alpha);
    ratio1   = TV/LL;
    R = [R; ratio1];
    scatter(theta,LL,'r')
    scatter(theta,TV,'b')
    scatter(theta,UL,'g')
    scatter(theta,GV,'c')
    scatter(theta,ratio1,'k')

    hold on
end
ratio1 = mean(R);

%% Find Constant C
mode=0.4;
R = [];
figure()
for theta=0.02:0.01:0.2
    UL      = (mode/theta)^(-alpha);
    LL      = (mode/theta + 1)^(-alpha);
    TV      = gamma(mode/theta+1-alpha)/gamma(mode/theta+1);

    ratio2   = TV/UL;
    R = [R ratio2];
    scatter(theta,LL,'r')
    scatter(theta,TV,'b')
    scatter(theta,UL,'g')
    scatter(theta,ratio2,'k')

    hold on
end
ratio2 = mean(R);

clearvars mode

%% 

%% Test with different g 
G = [0.65, 0.70, 0.72,0.8];
ls = {'-.','-','--',':'};
alpha=0.2;
beta=0.7;

figure()

for i=1:length(G)
g=G(i);
v=10;
syms M theta
eqn1 = (g^2*(2+4*alpha^2-6*alpha)/(2^(2*alpha)*theta^(2*alpha)*beta^2*(M/theta+1/2-alpha)^(2*alpha)*M)-2/M)*theta^2 ...
    +((g^2/(2^(2*alpha)*theta^(2*alpha)*beta^2*(M/theta+1/2-alpha)^(2*alpha)))*(3-4*alpha)-3)*theta ...
    + (g^2/(2^(2*alpha)*theta^(2*alpha)*beta^2*(M/theta+1/2-alpha)^(2*alpha)*M^(-1))-M) ==0;
eqn2 = v == 11/(beta*2^(2*alpha))*(M/theta+0.5-alpha/2)^alpha*theta^(1-alpha)*(M/theta+1-alpha);
[S1, S2] = vpasolve([eqn1 eqn2],[M theta],[0.2 1; 0 1]);
sprintf('g-ratio: %.2f, mode: %.2f, theta: %.2f',g,S1,S2)


% fimplicit3(eqn1,'FaceColor','y')
% hold on
% fimplicit3(eqn2,'FaceColor','r')
% xlim([0.1,1])
% ylim([0,1])
% xlabel('M')
% ylabel('theta')

fimplicit(eqn1,'Color','b','LineStyle',ls{i})
hold on
fimplicit(eqn2,'Color','r','LineStyle',ls{i})
hold on
xlim([0.1,1])
ylim([0,1])
xlabel('M')
ylabel('theta')
end
legend('g ratio fun','vel fun')

%% Use g from data
alpha=0.1;
A = 2^(2*alpha)*1.05; %
C = 2^(-alpha)*0.98; %
g_map = gratio.*tract;
R = 5.5*2;

for x=1:96
    for y=1:106
        for z=1:70
            if tract(x,y,z)== 1
                g=g_map(x,y,z);
                v=8; 
                 syms M theta
%                 eqn1 = ((A*(beta^2)*M^(2*alpha))*(M+3*((((v*beta*M^alpha)/(R*C))-M)/(1-alpha))+(2/M)*((((v*beta*M^alpha)/(R*C))-M)/(1-alpha))^2)/(M+(3-4*alpha)*((((v*beta*M^alpha)/(R*C))-M)/(1-alpha))+(((2+4*alpha^2-6*alpha)/M)*((((v*beta*M^alpha)/(R*C))-M)/(1-alpha))^2)))-g^2 ==0;
%                 V = vpasolve(eqn1,M,[0.1;0.6]);
%                 Mode=double(V);
%                 if isempty(Mode)
%                     sprintf('Problem with g-ratio: %.2f',g)
%                 else
%                 theta = ((beta*v*Mode^alpha)/(R*C) - Mode)/(1-alpha);
                eqn1 = (g^2*(2+4*alpha^2-6*alpha)/(2^(2*alpha)*theta^(2*alpha)*beta^2*(M/theta+1/2-alpha)^(2*alpha)*M)-2/M)*theta^2 ...
    +((g^2/(2^(2*alpha)*theta^(2*alpha)*beta^2*(M/theta+1/2-alpha)^(2*alpha)))*(3-4*alpha)-3)*theta ...
    + (g^2/(2^(2*alpha)*theta^(2*alpha)*beta^2*(M/theta+1/2-alpha)^(2*alpha)*M^(-1))-M) ==0;
                eqn2 = v == 11/(beta*2^(2*alpha))*(M/theta+0.5-alpha/2)^alpha*theta^(1-alpha)*(M/theta+1-alpha);
                [S1, S2] = vpasolve([eqn1 eqn2],[M theta],[0.2 0.6; 0 1]);
                Mode = S1;
                theta = S2;
                mean = Mode + theta;
                sprintf('g-ratio: %.2f, mode: %.2f, theta: %.2f, mean: %.2f',g,Mode,theta, mean)
                %end
            end
        end 
    end
end

%% Equation parts
M=0.1:0.01:0.9;
g=0.72;
v=4;
beta=0.7;
C =0.98;
pt1 = (M.^2)*(v*g^2)*(-4*alpha^2+8*alpha-1)/(11*C*A*beta);
pt2 = - (M.^(2*alpha+1))*alpha*(alpha+1);
pt3 = (M.^(3*alpha))*v*beta*(3*alpha+1)/(11*C);
pt4 = (M.*g^2*alpha*(alpha-1))/(A*beta^2);
pt5 = (M.^(2*alpha-1))*(g^2/(A*beta^2))*((v*beta/(11*C))^2)*(2+4*alpha^2-6*alpha);
pt6 = -(M.^(4*alpha-1))*2*((v*beta)/(11*C))^2;
total = pt1+pt2+pt3+pt4+pt5+pt6;
total2 = pt3 + pt5 + pt6;

close all
figure()
plot(M,pt1);
hold on
plot(M,pt2);
hold on
plot(M,pt3);
hold on
plot(M,pt4);
hold on
plot(M,pt5);
hold on
plot(M,pt6);
hold on
plot(M,total,'k','LineWidth',2);
% plot(M,total2,'k--','LineWidth',2);
hold on


g=0.71;
v=8;
pt1 = (M.^2)*(v*g^2)*(-4*alpha^2+8*alpha-1)/(11*C*A*beta);
pt2 = - (M.^(2*alpha+1))*alpha*(alpha+1);
pt3 = (M.^(3*alpha))*v*beta*(3*alpha+1)/(11*C);
pt4 = (M.*g^2*alpha*(alpha-1))/(A*beta^2);
pt5 = (M.^(2*alpha-1))*(g^2/(A*beta^2))*((v*beta/(11*C))^2)*(2+4*alpha^2-6*alpha);
pt6 = -(M.^(4*alpha-1))*2*((v*beta)/(11*C))^2;
total = pt1+pt2+pt3+pt4+pt5+pt6;
plot(M,total,'k--','LineWidth',2);
legend()

%% Two equations similar?

alpha=0.1;
beta=0.7;
theta=0.1;
g=0.6;
syms M
eqn1 = M*2*gamma(M/0.1)==3;
eqn1 = gamma(M)==1.7725;
syms M theta
eqn1 = 0 == g^2 - A*beta^2*M^(2*alpha)*(M+3*theta+2*theta^2/M)/(M+(3-4*alpha)*theta+(2+4*alpha^2-6*alpha)*theta^2/M);

eqn2 = v == (11/beta)*theta(1-alpha) ...
    * 0.98*(M/theta)^(alpha) ... % Lower limit and factor A
    * (M/theta+1-alpha);
[S1, S2] = vpasolve([eqn1 eqn2],[M theta],[0 1; 0 1]);
syms v
M =0.22;
theta=0.1;
eqn = theta==((v.*beta*M^alpha/(11*C))-M)/(1-alpha);
V = vpasolve(eqn,v)
v=5.18;
M=0.22;
syms theta
eqn = theta==((v.*beta*M^alpha/(11*C))-M)/(1-alpha);
V = vpasolve(eqn,theta)
theta=0.1;
syms M
eqn1 = (M^2)*(v.*g.^2)*(-4*alpha^2+8*alpha-1)/(11*C*A*beta) - (M^(2*alpha+1))*alpha*(alpha+1) + (M^(3*alpha))*v.*beta*(3*alpha+1)/(11*C) + (M*g.^2*alpha*(alpha-1))/(A*beta^2) + (M^(2*alpha-1))*(g.^2/(A*beta^2))*(v.*beta/(11*C))^2*(2+4*alpha^2-6*alpha) - (M^(4*alpha-1))*2*(v.*beta/(11*C))^2 == 0;
V = vpasolve(eqn1,M,[0;1])

g = 0.7
xt = @(M) (-g + sqrt(A*beta^2*M^(2*alpha)*(M+3*theta+2*theta^2/M)/(M+(3-4*alpha)*theta+(2+4*alpha^2-6*alpha)*theta^2/M)));
figure()
fplot(xt,[0,1])
g = 0.72
xt = @(M) (-g + sqrt(A*beta^2*M^(2*alpha)*(M+3*theta+2*theta^2/M)/(M+(3-4*alpha)*theta+(2+4*alpha^2-6*alpha)*theta^2/M)));
hold on
fplot(xt,[0,1])
beta=0.7;
xt2 = @(M) ((M^2)*(v.*g.^2)*(-4*alpha^2+8*alpha-1)/(11*C*A*beta) - (M^(2*alpha+1))*alpha*(alpha+1) + (M^(3*alpha))*v.*beta*(3*alpha+1)/(11*C) + (M*g.^2*alpha*(alpha-1))/(A*beta^2) + (M^(2*alpha-1))*(g.^2/(A*beta^2))*(v.*beta/(11*C))^2*(2+4*alpha^2-6*alpha) - (M^(4*alpha-1))*2*(v.*beta/(11*C))^2);
hold on
fplot(xt2,[0,1],'r')

figure
for theta=0.1:0.01:0.4
M=0.3;
g = sqrt(A*beta^2*M^(2*alpha)*(M+3*theta+2*theta^2/M)/(M+(3-4*alpha)*theta+(2+4*alpha^2-6*alpha)*theta^2/M));
scatter(theta,g)
hold on
end

theta=0.1;
M=0.22;
k=M/theta+1;
g=sqrt(A*beta^2*((k-1)*theta)^(2*alpha)*(k+1)*k/((k-2*alpha+1)*(k-2*alpha)))

%% 

A = 2^(2*alpha)*1; %
C = 2^(-alpha)*1; %
g=0.70;
v=5;
R = 5.5*2;
syms M  
xt = @(M) (((A*(beta^2)*M^(2*alpha))*(M+3*((((v*beta*M^alpha)/(R*C))-M)/(1-alpha))+(2/M)*((((v*beta*M^alpha)/(R*C))-M)/(1-alpha))^2)/(M+(3-4*alpha)*((((v*beta*M^alpha)/(R*C))-M)/(1-alpha))+(((2+4*alpha^2-6*alpha)/M)*((((v*beta*M^alpha)/(R*C))-M)/(1-alpha))^2)))-g^2);
figure()
fplot(xt,[0,1])
xlim([0 1])
ylim([-0.5 0.6])
xlabel('Mode')
ylabel('Function of mode')

eqn1= ((A*(beta^2)*M^(2*alpha))*(M+3*((((v*beta*M^alpha)/(11*C))-M)/(1-alpha))+(2/M)*((((v*beta*M^alpha)/(11*C))-M)/(1-alpha))^2)/(M+(3-4*alpha)*((((v*beta*M^alpha)/(11*C))-M)/(1-alpha))+(((2+4*alpha^2-6*alpha)/M)*((((v*beta*M^alpha)/(11*C))-M)/(1-alpha))^2))-g^2) ==0;
V = vpasolve(eqn1,M,[0.2;0.5])

alpha=0.2
g=0.6
v=8
value=[];
for M=0.1:0.1:1
    value=[value; (((A*(beta^2)*M^(2*alpha))*(M+3*((((v*beta*M^alpha)/(R*C))-M)/(1-alpha))+(2/M)*((((v*beta*M^alpha)/(R*C))-M)/(1-alpha))^2)/(M+(3-4*alpha)*((((v*beta*M^alpha)/(R*C))-M)/(1-alpha))+(((2+4*alpha^2-6*alpha)/M)*((((v*beta*M^alpha)/(R*C))-M)/(1-alpha))^2)))-g^2);];
end
figure()
plot(0.1:0.1:1,value)

Mode=0.3;
alpha=0.2;
beta=0.8;
v=8;
theta = ((beta*v*Mode^alpha)/(5.5*2*C) - Mode)/(1-alpha);
X = 0:0.01:3;
pdf = gampdf(X,Mode/theta+1,theta);
figure()
plot(X,pdf)

g=0.65;
v=6;
value=[];
for Mode=0.1:0.1:1
    value=[value; ((beta*v*Mode^alpha)/(5.5*2*C) - Mode)/(1-alpha)];
end
figure()
plot(0.1:0.1:1,value)

Mode=double(V);

%%
figure()
for M=0.1:0.01:0.6
theta=((v.*beta*M^alpha/(11*C))-M)/(1-alpha)
scatter(M,theta)
hold on
end

%%
g=0.70;
M=0.3;
alpha=0.1;
beta=0.8;
A= 1.05*2^(2*alpha);
C = 2^(-alpha)*0.98;
% a = g^2*(2+4*alpha^2-6*alpha)/(A*beta^2*M^(2*alpha+1))-2/M;
% b = (g^2/(A*beta^2*M^(2*alpha)))*(3-4*alpha)-3;
% c = (g^2/(A*beta^2*M^(2*alpha-1))-M);
% theta1 = (-b + sqrt(b^2-4*a*c))/(2*a)
% theta2 = (-b - sqrt(b^2-4*a*c))/(2*a)
v=8;
syms M theta
eqn1 = (g^2*(2+4*alpha^2-6*alpha)/(A*beta^2*M^(2*alpha+1))-2/M)*theta^2 ...
    +((g^2/(A*beta^2*M^(2*alpha)))*(3-4*alpha)-3)*theta ...
    + (g^2/(A*beta^2*M^(2*alpha-1))-M) ==0;
eqn2 = (v*beta*M^alpha)/(11*C) - M - (1-alpha)*theta ==0;
figure()
%fimplicit3([eqn1 eqn2])
fimplicit(eqn1,'Color','b')
hold on
fimplicit(eqn2,'Color','r')
hold on
xlim([0.1,1])
ylim([0,1])
xlabel('M')
ylabel('theta')

[S1, S2] = vpasolve([eqn1 eqn2],[M theta],[0.2 1; 0 1]);
%%
M=0.3;
g=0.7;
beta=0.7;
syms theta
eqn1 = (g^2*(2+4*alpha^2-6*alpha)/(A*beta^2*M^(2*alpha+1))-2/M)*theta^2 ...
    +((g^2/(A*beta^2*M^(2*alpha)))*(3-4*alpha)-3)*theta ...
    + (g^2/(A*beta^2*M^(2*alpha-1))-M) ==0;
S = vpasolve(eqn1, theta,[0,1])
theta=0.5;
v=8;
eqn2 = (v*beta*M^alpha)/(11*C) - M - (1-alpha)*theta ==0;

% fun1 = @(M,theta) ((g^2*(2+4*alpha^2-6*alpha)/(A*beta^2*M^(2*alpha+1))-2/M)*theta^2 ...
%     +((g^2/(A*beta^2*M^(2*alpha)))*(3-4*alpha)-3)*theta ...
%     + (g^2/(A*beta^2*M^(2*alpha-1))-M));
% figure()
% zhandle = ezsurf(fun1,[0.1,1])
% fun2 = @(M,theta) ((v*beta*M^alpha)/(11*C) - M - (1-alpha)*theta)
% hold on
% zhandle = ezsurf(fun2,[0.1,1])
% 
% [S1, S2] = vpasolve([eqn1 eqn2],[M theta],[0.2 1; 0 1]);
% vpa(S1)
% vpa(S2)