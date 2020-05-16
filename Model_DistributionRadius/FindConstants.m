%% Find alpha and beta

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
for theta=0.02:0.01:0.5
    LL      = (mode/theta)^(2*alpha);
    UL      = (mode/theta + 1)^(2*alpha);
    TV      = gamma(mode/theta+1)/gamma(mode/theta+1-2*alpha);
 
    ratio1   = TV/LL;
    R = [R; ratio1];
    
    scatter(theta,LL,'r')
    scatter(theta,TV,'b')
    scatter(theta,UL,'g')
    scatter(theta,ratio1,'k')

    hold on
end

ratio1 = mean(R);
xlabel('Theta')
ylabel('Value of limit')
title(sprintf('Ratio between 2 gammas for mode = %.1f - Equation 1 (g)',mode))
legend('Ratio TV/LL','Lower Limit (LL)','True Value (TV)','Upper Limit (UL)')

%% Find Constant A (chaging mode and theta)

R =[];

figure()
for mode=0.1:0.1:0.4
for theta=0.02:0.01:0.5
    LL      = (mode/theta)^(2*alpha);
    UL      = (mode/theta + 1)^(2*alpha);
    TV      = gamma(mode/theta+1)/gamma(mode/theta+1-2*alpha);
 
    ratio1   = TV/LL;
    R = [R; ratio1];
    
    scatter3(mode,theta,ratio1,'r')
    hold on
end
end 

ratio1 = mean(R);
xlabel('Mode')
ylabel('Theta')
zlabel('Ratio')
title('Ratio between 2 gammas')

%%

alpha = 0.1;
A = 1.02;
lgn = cell(12,1);
k = 1;

figure()
for A=0.98:0.02:1.2
syms mode theta
eq1 = (mode/theta)^(2*alpha)*A  == gamma(mode/theta+1)/gamma(mode/theta+1-2*alpha);
symvar(eq1);
fimplicit(eq1);
hold on
lgn{k} = strcat('Bias=',num2str(A));
k = k+1;
end

xlim([0 1])
ylim([0 1])
xlabel('M')
ylabel('theta')
legend(lgn)

%% Find Constant C
mode=0.1;
R = [];
figure()
for theta=0.02:0.01:0.5
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
xlabel('Theta')
ylabel('Value of limit')
title(sprintf('Ratio between 2 gammas for mode = %.1f - Equation 2 (v)',mode))
legend('Ratio TV/LL','Lower Limit (LL)','True Value (TV)','Upper Limit (UL)')
