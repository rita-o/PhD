close all
clc
figure()
% Display gamma distribution
theta       = 0.2;
mode        = 0.3;
shape       = mode/theta+1;
X           = (0:0.01:4);
pdf         = gampdf(X,shape,theta);

subplot(2,1,1)
t   =  plot(2*X,pdf,'r');
title('Distribution of axons in fiber')
xlabel('Axon diameter (um)')
ylabel('pdf value')

% Display gamma distribution
theta       = 0.2;
mode        = 0.3;
shape       = mode/theta+1;
X           = (0:0.01:4);
pdf         = gampdf(X,shape,theta);

subplot(2,1,2)
t   =  plot(X,pdf,'r');
title('Distribution of axons in fiber')
xlabel('Axon radius (um)')
ylabel('pdf value')
%%
close all
figure

  
    %fimplicit(@(M,theta) -v + ((2*5.5)/beta)*theta^(1-alpha) ...
   % * gamma(M/theta+1-alpha)/gamma(M/theta+1) ... % Without approximmation
    %* (M/theta+1-alpha,'Color',[0, 0.4470, 0.7410]);
    fimplicit(@(M,theta) myfun(M,theta,v,alpha,beta));
     warning('off')
    %res     = lsqnonlin(@(M,theta) make_equations_inc_gv(M,theta,alpha,beta,G,V),x0,[0 0],[1 1],options); % other option

    hold on
    xlim([0,1])
    ylim([0,1])
    xlabel('M')
    ylabel('theta')
    hold on
    syms M theta
   % fimplicit(v_eqn_GammaLevel,'Color','r')
    
    
%%

function F = myfun(M,theta,v,alpha,beta)

F = [v - ((2*5.5)/beta)*theta^(1-alpha) ...
    * gamma(M/theta+1-alpha)/gamma(M/theta+1) ... % Without approximmation
    * (M/theta+1-alpha)];
end
