%% MASS SPRING DAMPER IDENTIFICATION
%Identification of the parameters of a mass spring damper.
clc
clear all
close all

% define system properties
m = 1; % mass, kg
k = 10; % stiffness, N/m
lambda = 0.75; % damping

sim_step = 0.02;
t = 0:sim_step:10; % time vector
v = [t' 0.1*randn(numel(t),3)]; % measurement noise
u = 10*[zeros(1,length(0:sim_step:1)) ones(1,length(1.02:sim_step:t(end)))]; % step 

% run simulation
model = 'mass_spring_damper.slx';
[t,x,y]= sim(model,[t(1) t(end)]);

%Plot

y_name={'d2x/dt2 (m/s^2)','dx/dt (m/s)','x (m)','u (N)'};
figure
for i=1:4
    subplot(4,1,i)
    plot(t,y(:,i))
    ylabel(y_name{i});xlabel('Time (s)'); grid
end

% plot system phase plane diagram
figure
plot(x(:,1),x(:,2))
ylabel('x(m)')
xlabel('dx/dt(m/s)')

%% Linear regression
% Perform a linear regression with normalized equation

%First define the model 
% m*d2x/dt2 + z(damping)*dx/dt + kx = u => eq: y = m*x_dot_dot = u - z*x_dot - k*x
% Model parameters, defined as: theta = [1 z k] 

N=length(t);
X = [y(:,4) -y(:,2) -y(:,3)];%Define the regressors [u -x_dot -x]
z = m*y(:,1); %Define the measured output m*x_dot_dot 

%Minimum of least squares cost function occurs at theta_est = (X'X)^-1 X'z
D = inv(X'*X);
theta_est = D*X'*z;%estimation of theta
y_est = X*theta_est; %output estimation
res = z - y_est;%residual 

sigma2=sqrt(res'*res/(N-3));%variance
%Define confidence interval for theta. 95% of theta lying at 2 std deviations
theta_min = theta_est - 2*sqrt(sigma2*diag(D));
theta_max = theta_est + 2*sqrt(sigma2*diag(D));
y_min=zeros(N,1);y_max=zeros(N,1);
%Obtain the 95% confidence interval for the estimated output
for i=1:N
    y_min(i) = y_est(i) - 2*sqrt(sigma2*X(i,:)*D*X(i,:)');
    y_max(i) = y_est(i) + 2*sqrt(sigma2*X(i,:)*D*X(i,:)');
end

%Plot
figure
plot(t,[z,y_est,y_min,y_max],'LineWidth',0.5)
legend('Experimental Data','Estimated Data','Lower Confidence Bound','Upper Confidence Bound');xlabel('Time (s)'); grid
title('d2x/dt2 (m/s^2)');

%To compare the effect of noise visualize the estimates of theta and the
%'real' values. Measurement noise should negatively affect estimates of the estimated values
disp(strcat('Real values are: ',num2str([m,lambda,k])))
disp(strcat('Estimated values are: ',num2str(theta_est')))

%% OPTIMIZATION
%Obtain cost function topology for mass-spring-damper. Objective is the
%reduction of  variance in the system’s response to random noise.
%Cost function defined as: 
% J = var(X) + var(X_dot) +g(lambda) + g(m) where g(z)=z^2 + exp(1-0.1z)

clear all

%Domain of the different parameters
lambda_set = 0:0.02:1;
m_set = 0.1:0.02:1;
k=1;

%Initialize parameters
sim_step = 0.02;
t = 0:sim_step:10; % time vector
v = [t' 0*randn(numel(t),3)]; %No noise
rng(1); % set random number seed
u = randn(size(t));
model = 'mass_spring_damper.slx';

%Calculate every point in the domain
tic;
for i=1:numel(lambda_set)
    for j=1:numel(m_set)
        m=m_set(j); lambda= lambda_set(i);
        [~,X,~] = sim(model);
        J(i,j) = var(X(:,1))+var(X(:,2))+lambda_set(i)^2+exp(1-0.1*lambda_set(i))+m_set(j)^2+exp(1-0.1*m_set(j));
    end
end
time = toc;

% find minima in the whole domain
[a,b]=min(J(:));
[R,C] = ind2sub(size(J),b);
lambda_0 = lambda_set(R);
m_0 = m_set(C);

% plot contours of cost function value
[X,Y]=meshgrid(lambda_set,m_set);
figure
subplot(1,2,1)
plot(lambda_set,lambda_set.^2+exp(1-0.1*lambda_set)); grid on
xlabel('\lambda')
ylabel('g(\lambda)')

subplot(1,2,2)
contour(X,Y,J',50); hold on; grid on
plot(lambda_0,m_0,'ro','MarkerFaceColor','r')
xlabel('\lambda'); ylabel('m');
text(lambda_0+0.02,m_0,['minima J = ' num2str(J(R,C))])
text(lambda_0+0.02,m_0+0.04,['\lambda = ' num2str(lambda_0)])
text(lambda_0+0.02,m_0+0.08,['m = ' num2str(m_0)])

