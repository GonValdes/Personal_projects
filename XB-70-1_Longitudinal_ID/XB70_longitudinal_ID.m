%% Identifying longitudinal pitching moment derivatives of the XB-70-1
% Study the effectiveness and the sensitivity to different inputs of
% ordinary least squares method for identifying pitching moment
% derivatives.

% Part of the Aircraft system Identification module labs in Cranfield University.
% Data extracted from Technical Report NASA TN D-4578

% Code written by Gonzalo Valdés, Feb 2020

clc
clear all
close all

% Aircraft parameters and flight conditions
S = 6297.8;% ft2
c_bar = 78.53;% ft
Iy = 10*10^6;% slug/ft2
V0 = 1936;% ft/s
q_bar = 424;% dynamic pressure lb/ft2
alpha0 = 6.2*pi/180; % trim alpha rad

T = 0:0.01:30;% time vector

% Noise definition
N = 0.1*[0.1 0.1 0.002 0.001 0.001 0.001 1];% [u w q qdot alpha alphadot eta]

% Input definition
%Different inputs are tested for better identification
U(:,1) = -[zeros(1,100) ones(1,50) -ones(1,50) ones(1,50) zeros(1,2751)];% 1 2 1
U(:,2) = -[zeros(1,100) ones(1,50) -ones(1,50) zeros(1,2801)]; % doublet
U(:,3) = -[zeros(1,100) ones(1,300) -ones(1,200) ones(1,100) ...
    -ones(1,100) ones(1,100) zeros(1,2101)];% 1 2 1 1 1

% Aircraft model: XB-70 60,000ft Mach 2.00 
% Modelled as transfer functions.
q_delE = zpk([0 -0.00174 -0.175],[-0.0743 -0.0181],-1.62)...
    *tf(1,[1 2*0.145*1.93 1.93^2]); %elevator to pitch rate
u_delE = zpk(-226,[-0.0743 -0.0181],1.50)...
    *tf([1 2*0.832*0.164 0.164^2],[1 2*0.145*1.93 1.93^2]); %forward speed to pitch rate
w_delE = zpk([-226 -0.0196 -0.01],[-0.0743 -0.0181],-13.8)...
    *tf(1,[1 2*0.145*1.93 1.93^2]);%vertical speed to pitch rate

% Run simulations and add measurement and input noise
for i = 1:size(U,2)
    us = (2*U(:,i)+N(7)*randn(numel(T),1))*pi/180;% define input in rads
    
    q = lsim(q_delE,us,T)+N(3)*randn(numel(T),1);% rad
    qdot = gradient(q,T);
    
    u = lsim(u_delE,us,T)+N(1)*randn(numel(T),1)+V0*cos(alpha0);% ft/s
    w = lsim(w_delE,us,T)+N(2)*randn(numel(T),1)+V0*sin(alpha0);% ft/s
    
    alpha = atan(w./u) + N(5)*randn(numel(T),1);%AoA (rad)
    alphadot = gradient(alpha,T) + N(6)*randn(numel(T),1);%dAoA/dt (rad/s)
    % Write response data to Excel file
    XY{i} = [T' us q qdot alpha alphadot u w];
%     sheet_name = ['Data_' num2str(i)];
%     xlswrite('Flight_data_XB70',XY{i},sheet_name);
end

% Visualise simulation data
h = figure;
for i=1:size(U,2)
    % plot response
    subplot(2,size(U,2),i)
    plot(T,XY{i}(:,5)*180/pi)
    %     axis([0 max(T) -15 15])
    xlabel('Time (s)'); ylabel('q (deg/s)')
    % plot elevator deflection
    subplot(2,size(U,2),i+size(U,2))
    plot(T,XY{i}(:,2)*180/pi)
    %     axis([0 max(T) -20 20])
    xlabel('Time (s)'); ylabel('\eta (degrees)')
end

%% Identification routine
% H&J data: CMa ~ -0.2, Cmq ~ -1, Cmde ~ -0.5
% Cm*B = dq/dt = B(Cm0 + Cma*da + C*Cmq*q + Cmde*de)

B = q_bar*S*c_bar/Iy;
C = c_bar/(2*V0);

for i=1:size(U,2)
    %XY{i} = [T' U(:,i) q qdot alpha alphadot];
    %Regressors: x = B[1 da C*q u] 
    x = B*[ones(numel(T),1) (XY{i}(:,5)-alpha0) C*XY{i}(:,3) XY{i}(:,2)];
    z = XY{i}(:,4);% qdot
    theta_est(:,i) = inv(x'*x)*x'*z;
    
    qdot_est = B*[ones(numel(T),1) (XY{i}(:,5)-alpha0) C*XY{i}(:,3) XY{i}(:,2)]...
        *theta_est(:,i);
    v = z-qdot_est;
    XY{i}(:,end+1) = qdot_est;
    
    sigma2(:,i) = sum(v'*v)/(size(x,1)-size(x,2));%fit error
    cov_est{i} = (inv(x'*x)*x')*sigma2(:,i)*(x*inv(x'*x));% covariance matrix
    var_est(:,i) = sigma2(i)*diag(inv(x'*x));
end

%Validation run. Simulate quality of estimated parameters by simulating
%response to the validation input not previously used.
uv = -[zeros(1,100) ones(1,100) zeros(1,2801)]'+N(7)*randn(numel(T),1);%Pulse input
q = lsim(q_delE,uv,T);% rad
qdot = gradient(q,T);% rad/s

u = lsim(u_delE,uv,T)+V0*cos(alpha0);% ft/s
w = lsim(w_delE,uv,T)+V0*sin(alpha0);% ft/s
alpha = atan(w./u);%rad
% estimated qdot
qdot_est = [ones(size(alpha)) B*(alpha-alpha0) B*C*q B*uv]*theta_est(:,3);%rad/s

%Plot data
figure
subplot(3,1,1)
plot(T,qdot/B,'r','Linewidth',1);
hold on
plot(T,qdot_est/B,'k--','Linewidth',1)
legend('Estimated','Model')
xlabel('Time (s)')
ylabel('C_m')
subplot(3,1,2)
plot(T,uv,'k','Linewidth',1)
xlabel('Time (s)')
ylabel('Input: \eta (deg)')
subplot(3,1,3)
plot(T,(qdot-qdot_est)/B,'k','Linewidth',1)
xlabel('Time (s)')
ylabel('Residual')

%% Figure with inputs

figure
subplot(2,2,1)
plot(T,U(:,1),'k','Linewidth',2); hold on
axis([0 max(T)/2 -2 2])
ylabel('\eta (deg)'); xlabel('Time (s)')
title('Test input 1')
subplot(2,2,2)
plot(T,U(:,2),'k','Linewidth',2);
axis([0 max(T)/2 -2 2])
ylabel('\eta (deg)'); xlabel('Time (s)')
title('Test input 2')
subplot(2,2,3)
plot(T,U(:,3),'k','Linewidth',2); 
axis([0 max(T)/2 -2 2])
ylabel('\eta (deg)'); xlabel('Time (s)')
title('Test input 3')
subplot(2,2,4)
plot(T,uv,'k','Linewidth',2); 
axis([0 max(T)/2 -2 2])
ylabel('\eta (deg)'); xlabel('Time (s)')
title('Validation input')


%% Parameter estimates and variances
color = ['ko';'r^';'bs'];
figure
for i=1:3
    fig_n{i} = plot(1:size(theta_est,1),theta_est(:,i),color(i,:),'MarkerFaceColor',color(i,1))
    hold on
    parx = [[1:size(theta_est,1)]+0.1;[1:size(theta_est,1)]-0.1];
    %Upper confidence interval prediction
    plot(parx,[[theta_est(:,i)+2*sqrt(sigma2(:,i))]'; [theta_est(:,i)+2*sqrt(sigma2(:,i))]'],[color(i,1) '-'],'MarkerFaceColor','b')
    %Lower confidence interval prediction
    plot(parx,[[theta_est(:,i)-2*sqrt(sigma2(:,i))]'; [theta_est(:,i)-2*sqrt(sigma2(:,i))]'],[color(i,1) '-'],'MarkerFaceColor','b')
end
set(gca,'XLim',[0 5])
ylabel('Value')
xlabel('Cm0  Cma  Cmq  Cmde')
legend([fig_n{1},fig_n{2},fig_n{3}],{'1-2-1','Doublet','1-2-1-1-1'})
