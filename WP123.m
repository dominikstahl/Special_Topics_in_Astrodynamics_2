clear;
clc;
close all;

%% Constants and Parameters

m_s = 1.9891e30; % Mass of Sun [kg]
m_m = 641.85e21; % Mass of Mars [kg]
%m_m=6.0477e24;%Earth
a = 227939200; % Mars semi-major axis [km] % dimensional distance unit [km]
AU = 149597870.7; % [km]
global mu; mu = m_m/(m_s+m_m); % Parameter of cr3pb [-]
x_s = -mu; % x-coordinate of Sun in xyz-system [-]
x_m = 1-mu; % x-coordinate of Mars in xyz-system [-]
global beta; beta = .1; % lightness number
%G = 6.67430e-20; % gravitational constant [km^3/kgs^2]
T_mars = (a/AU)^1.5; % T Mars [Earthyears]
t_unity = T_mars/(2*pi); % dimensional time unit [Earthyears]
t_units = t_unity*365.25*24*60*60; % dimensional time unit [s]
%clear G omega;

%% Task 1.1a

% Finding L1:
new_x_L1 = 1-2*mu; % initial guess for L1
while true
    x = new_x_L1;
    f_of_x = compute_Ux(x,0,0);
    f_prime_of_x = compute_Uxx(x,0,0);
    new_x_L1 = x-f_of_x/f_prime_of_x;
    if(abs(compute_Ux(new_x_L1,0,0))<1e-15)
        break;
    end
end
x_L1 = new_x_L1;
clear x f_of_x f_prime_of_x new_x_L1

% Finding L3:
new_x_L3 = -1; % initial guess for L1
while true
    x = new_x_L3;
    f_of_x = compute_Ux(x,0,0);
    f_prime_of_x = compute_Uxx(x,0,0);
    new_x_L3 = x-f_of_x/f_prime_of_x;
    if(abs(compute_Ux(new_x_L3,0,0))<1e-15)
        break;
    end
end
x_L3 = new_x_L3;
clear x f_of_x f_prime_of_x new_x_L3

%% Task 1.1b

% According to Wakker eq. (3.64)

% L1
f_of_L1 = x_L1-(1-mu)*(mu+x_L1)/abs(mu+x_L1)^3+mu*(1-mu-x_L1)/abs(1-mu-x_L1)^3;

% L3
f_of_L3 = x_L3-(1-mu)*(mu+x_L3)/abs(mu+x_L3)^3+mu*(1-mu-x_L3)/abs(1-mu-x_L3)^3;

%% Task 1.2a

shift = .005;

% L1
x_shifted_L1 = x_L1-shift;
r1 = abs(x_shifted_L1-x_s);
needed_a_at_shifted_L1 = abs(compute_Ux(x_shifted_L1,0,0)); % [m/s^2]
beta_shifted_L1 = needed_a_at_shifted_L1*r1^2/(1-mu);
clear r1

% L3
x_shifted_L3 = x_L3+shift;
r1 = abs(x_shifted_L3-x_s);
needed_a_at_shifted_L3 = abs(compute_Ux(x_shifted_L3,0,0)); % [m/s^2]
beta_shifted_L3 = needed_a_at_shifted_L3*r1^2/(1-mu);
clear r1

%% Task 1.2c

beta_for_rest_of_code = beta; % beta to be reseted in the end

% shifted L1
beta = beta_shifted_L1;
global alpha; alpha = 0;
global delta; delta = 0;
figure
for x=-10:10
    for y=-10:10
        [~,s]=ode45(@compute_State_Derivative_with_Sailing,[0 .5],[x_shifted_L1+x*1e-5;y*1e-5;0;0;0;0],odeset('RelTol',1e-12,'AbsTol',1e-12));
        plot(s(1,1),s(1,2),'.','linewidth',2,'markersize',15,'color',[226 26 26]/255); hold on;
        plot(s(:,1),s(:,2),'linewidth',2,'color',[0 166 214]/255,'displayname','unstable manifold'); hold on;
    end 
end
grid on;
xlabel('$x$','interpreter','latex');
ylabel('$y$','interpreter','latex');
axis equal;

% shifted L3
beta = beta_shifted_L3;
figure
for x=-10:10
    for y=-10:10
        [~,s]=ode45(@compute_State_Derivative_with_Sailing,[0 .5],[x_shifted_L3+x*1e-5;y*1e-5;0;0;0;0],odeset('RelTol',1e-12,'AbsTol',1e-12));
        plot(s(1,1),s(1,2),'.','linewidth',2,'markersize',15,'color',[226 26 26]/255); hold on;
        plot(s(:,1),s(:,2),'linewidth',2,'color',[0 166 214]/255,'displayname','unstable manifold'); hold on;
    end 
end
grid on;
xlabel('$x$','interpreter','latex');
ylabel('$y$','interpreter','latex');
axis equal;

beta = beta_for_rest_of_code;
clear beta_for_rest_of_code

%% Task 1.2d

solar_radiation_pressure = 4.5e-6; % [Pa] [N/m^2]
mass = 10; % [kg]

% L1
force1 = needed_a_at_shifted_L1*mass; % [kgm/s^2] [N]
area1 = force1/solar_radiation_pressure; % [m^2]
sidelength1 = sqrt(area1); % [m]
clear force1 area1

% L1
force3 = needed_a_at_shifted_L3*mass; % [kgm/s^2] [N]
area3 = force3/solar_radiation_pressure; % [m^2]
sidelength3 = sqrt(area3); % [m]
clear force3 area3

clear solar_radiation_pressure mass

%% Task 2.1a

% L1
A_L1 = [0,0,1,0;0,0,0,1;-compute_Uxx(x_L1,0,0),-compute_Uxy(x_L1,0,0),0,2;-compute_Uxy(x_L1,0,0),-compute_Uyy(x_L1,0,0),-2,0];
[vec1,eig1] = eig(A_L1);
eig1 = diag(eig1);

% L3
A_L3 = [0,0,1,0;0,0,0,1;-compute_Uxx(x_L3,0,0),-compute_Uxy(x_L3,0,0),0,2;-compute_Uxy(x_L3,0,0),-compute_Uyy(x_L3,0,0),-2,0];
[vec3,eig3] = eig(A_L3);
eig3 = diag(eig3);

%% Task 2.1b

% L1
[r1,r2] = get_r1r2(x_L1,0,0);
K = (1-mu)/r1^3+mu/r2^3;
Alpha = 1-1/2*K;
Beta_squared = (2*K+1)*(K-1);
Lambda1 = -Alpha+sqrt(Alpha^2+Beta_squared);
Lambda2 = -Alpha-sqrt(Alpha^2+Beta_squared);
lambda1_L1 = -sqrt(Lambda1);
lambda2_L1 = sqrt(Lambda1);
lambda3_L1 = -sqrt(Lambda2);
lambda4_L1 = sqrt(Lambda2);

% L1
[r1,r2] = get_r1r2(x_L3,0,0);
K = (1-mu)/r1^3+mu/r2^3;
Alpha = 1-1/2*K;
Beta_squared = (2*K+1)*(K-1);
Lambda1 = -Alpha+sqrt(Alpha^2+Beta_squared);
Lambda2 = -Alpha-sqrt(Alpha^2+Beta_squared);
lambda1_L3 = -sqrt(Lambda1);
lambda2_L3 = sqrt(Lambda1);
lambda3_L3 = -sqrt(Lambda2);
lambda4_L3 = sqrt(Lambda2);

clear K r1 r2 Alpha Beta_squared Lambda1 Lambda2

%% Task 2.2a

maxtime = 5/t_unity;

figure;
initialcondition = [x_L1;0;0;0;0;0]+1e-5*[vec1(1:2,2);0;vec1(3:4,2);0];
[~,s]=ode45(@compute_State_Derivative,[0 maxtime],initialcondition,odeset('RelTol',1e-12,'AbsTol',1e-12));
plot(s(:,1),s(:,2),'linewidth',2,'color',[0 166 214]/255,'displayname','unstable manifold'); hold on;
initialcondition = [x_L1;0;0;0;0;0]-1e-5*[vec1(1:2,2);0;vec1(3:4,2);0];
[~,s]=ode45(@compute_State_Derivative,[0 maxtime],initialcondition,odeset('RelTol',1e-12,'AbsTol',1e-12));
plot(s(:,1),s(:,2),'--','linewidth',2,'color',[0 166 214]/255,'displayname','unstable manifold'); hold on;
initialcondition = [x_L1;0;0;0;0;0]+1e-5*[vec1(1:2,1);0;vec1(3:4,1);0];
[~,s]=ode45(@compute_State_Derivative,[maxtime 0],initialcondition,odeset('RelTol',1e-12,'AbsTol',1e-12));
plot(s(:,1),s(:,2),'linewidth',2,'color',[230 70 22]/255,'displayname','stable manifold'); hold on;
initialcondition = [x_L1;0;0;0;0;0]-1e-5*[vec1(1:2,1);0;vec1(3:4,1);0];
[~,s]=ode45(@compute_State_Derivative,[maxtime 0],initialcondition,odeset('RelTol',1e-12,'AbsTol',1e-12));
plot(s(:,1),s(:,2),'--','linewidth',2,'color',[230 70 22]/255,'displayname','stable manifold'); hold on;
plot(-mu,0,'kp','linewidth',3,'markersize',10,'displayname','Sun');
plot(1-mu,0,'k.','linewidth',2,'markersize',30,'displayname','Mars');
plot(x_L1,0,'k+','linewidth',2,'markersize',10,'displayname','L1');
axis equal; grid on;
legend('interpreter','latex');
xlabel('$x$','interpreter','latex');
ylabel('$y$','interpreter','latex');
clear s;

%% Task 2.2b

maxtime = 27.25/t_unity;

figure;
initialcondition = [x_L1;0;0;0;0;0]+1e-5*[vec1(1:2,2);0;vec1(3:4,2);0];
[~,s]=ode45(@compute_State_Derivative,[0 maxtime],initialcondition,odeset('RelTol',1e-12,'AbsTol',1e-12));
plot(s(:,1),s(:,2),'linewidth',2,'color',[0 166 214]/255,'displayname','unstable manifold L1'); hold on;
initialcondition = [x_L3;0;0;0;0;0]+1e-5*[vec3(1:2,3);0;vec3(3:4,3);0];
[~,s]=ode45(@compute_State_Derivative,[0 maxtime],initialcondition,odeset('RelTol',1e-12,'AbsTol',1e-12));
plot(s(:,1),s(:,2),'linewidth',2,'color',[230 70 22]/255,'displayname','stable manifold L3'); hold on;
plot(-mu,0,'kp','linewidth',3,'markersize',10,'displayname','Sun');
plot(1-mu,0,'k.','linewidth',2,'markersize',30,'displayname','Mars');
plot(x_L1,0,'k+','linewidth',2,'markersize',10,'displayname','L1');
plot(x_L3,0,'kx','linewidth',2,'markersize',10,'displayname','L3');
axis equal; grid on;
legend('interpreter','latex');
xlabel('$x$','interpreter','latex');
ylabel('$y$','interpreter','latex');
clear s;

%% Task 2.3a

maxtime = 5/t_unity;

global delta; delta = 0; % clock angle
global alpha; % cone angle
alpha_deg = -70:2.5:70;
least_norm = 3*ones(size(alpha_deg));
best_position_error = zeros(size(alpha_deg));
best_velocity_error = zeros(size(alpha_deg));
time1 = zeros(size(alpha_deg));
time3 = zeros(size(alpha_deg));
for alpha_index=1:length(alpha_deg)
    alpha = alpha_deg(alpha_index)*pi/180;
    initialcondition = [x_L1;0;0;0;0;0];
    [t1,s1]=ode45(@compute_State_Derivative_with_Sailing,[0 maxtime],initialcondition,odeset('RelTol',1e-12,'AbsTol',1e-12));
    initialcondition = [x_L3;0;0;0;0;0];
    [t3,s3]=ode45(@compute_State_Derivative_with_Sailing,[0 -maxtime],initialcondition,odeset('RelTol',1e-12,'AbsTol',1e-12));
    i = 1;
    while i<=length(t1)
        j = 1;
        while j<=length(t3)
            if(norm(s1(i,:)-s3(j,:))<least_norm(alpha_index))
                least_norm(alpha_index) = norm(s1(i,:)-s3(j,:));
                best_position_error(alpha_index) = norm(s1(i,1:3)-s3(j,1:3));
                best_velocity_error(alpha_index) = norm(s1(i,4:6)-s3(j,4:6));
                time1(alpha_index) = t1(i);
                time3(alpha_index) = t3(j);
                j = j+1;
            elseif(norm(s1(i,:)-s3(j,:))<2*least_norm)
                j = j+1;
            elseif(norm(s1(i,:)-s3(j,:))<10*least_norm)
                j = j+10;
            else
                j = j+100;
            end
        end
        i = i+1;
    end
end
[~,best_alpha_index] = min(least_norm);

figure;
semilogy(alpha_deg,least_norm,'x','linewidth',2,'color',[0 166 214]/255,'displayname',''); hold on; grid on;
xlabel('cone angle $\alpha$','interpreter','latex');
ylabel('minimum Euclidean-norm error','interpreter','latex');

figure;
alpha = alpha_deg(best_alpha_index)*pi/180;
initialcondition = [x_L1;0;0;0;0;0];
[~,s1]=ode45(@compute_State_Derivative_with_Sailing,[0 time1(best_alpha_index)],initialcondition,odeset('RelTol',1e-12,'AbsTol',1e-12));
plot(s1(:,1),s1(:,2),'linewidth',2,'color',[0 166 214]/255,'displayname','unstable manifold L1'); hold on;
initialcondition = [x_L3;0;0;0;0;0];
[~,s3]=ode45(@compute_State_Derivative_with_Sailing,[0 time3(best_alpha_index)],initialcondition,odeset('RelTol',1e-12,'AbsTol',1e-12));
plot(s3(:,1),s3(:,2),'linewidth',2,'color',[230 70 22]/255,'displayname','stable manifold L3'); hold on;
plot(-mu,0,'kp','linewidth',3,'markersize',10,'displayname','Sun');
plot(1-mu,0,'k.','linewidth',2,'markersize',30,'displayname','Mars');
plot(x_L1,0,'k+','linewidth',2,'markersize',10,'displayname','L1');
plot(x_L3,0,'kx','linewidth',2,'markersize',10,'displayname','L3');
xlabel('$x$','interpreter','latex');
ylabel('$y$','interpreter','latex');
axis equal; grid on;
legend('interpreter','latex');
position_error = best_position_error(best_alpha_index)*a; % [km]
velocity_error = best_velocity_error(best_alpha_index)*a/t_units; % [km/s]
time_unstable = time1(best_alpha_index)*t_unity; % [Earth years]
time_stable = -time3(best_alpha_index)*t_unity; % [Earth years]

%% Task 3.1

numbering = [['1st'];['2nd'];['3rd'];['4th'];['5th'];['6th'];['7th'];['8th'];['9th']];

% L1
initial_State_guess = [0.994251260979694;0;0;0.00989534972689311]; % From exercise
figure;
for q = 1:10
    initial_Phi = eye(4);
    initialcondition = [initial_State_guess;initial_Phi(:,1);initial_Phi(:,2);initial_Phi(:,3);initial_Phi(:,4)];
    [t,s]=ode45(@compute_State_and_Phi_Derivative,[0 2],initialcondition,odeset('RelTol',1e-12,'AbsTol',1e-12));
    for i=1:length(t)
        if(s(i,2)<0)
            % find the t where y coordinate is 0 (linear interpolation)
            t_end = (t(i)*abs(s(i-1,2))+t(i-1)*abs(s(i,2)))/(abs(s(i-1,2))+abs(s(i,2)));
            break;
        end
    end
    [t,s]=ode45(@compute_State_and_Phi_Derivative,[0 t_end],initialcondition,odeset('RelTol',1e-12,'AbsTol',1e-12));
    plot(s(:,1),s(:,2),'linewidth',2,'displayname',[numbering(q,:),' try']); hold on;
    % Update guess:
    current_Derivative = compute_State_and_Phi_Derivative(0,s(end,:)'); % state derivative at t1/2
    Phi34 = s(end,19); % Phi_34 at t1/2
    Phi24 = s(end,18); % Phi_24 at t1/2
    xdotdot = current_Derivative(3); % xdotdot at t1/2
    ydot = s(end,4); % ydot at t1/2
    xdot = s(end,3); % xdot at t1/2
    dydot = -(Phi34-xdotdot/ydot*Phi24)^(-1)*xdot;
    if(abs(xdot)<1e-8)
        break;
    end
    initial_State_guess = initial_State_guess+[0;0;0;dydot];
end
T1 = 2*t_end; % final period L1
initial1 = initial_State_guess; % final initial guess L1
plot(1-mu,0,'k.','linewidth',2,'markersize',30,'displayname','Mars');
plot(x_L1,0,'k+','linewidth',2,'markersize',10,'displayname','L1');
axis equal; grid on;
legend('interpreter','latex');
xlabel('$x$','interpreter','latex');
ylabel('$y$','interpreter','latex');

% L3
initial_State_guess = [-1.00100081972868;0;0;0.00217804834346359]; % From exercise
figure;
for q = 1:10
    initial_Phi = eye(4);
    initialcondition = [initial_State_guess;initial_Phi(:,1);initial_Phi(:,2);initial_Phi(:,3);initial_Phi(:,4)];
    [t,s]=ode45(@compute_State_and_Phi_Derivative,[0 5],initialcondition,odeset('RelTol',1e-12,'AbsTol',1e-12));
    for i=1:length(t)
        if(s(i,2)<0)
            % find the t where y coordinate is 0 (linear interpolation)
            t_end = (t(i)*abs(s(i-1,2))+t(i-1)*abs(s(i,2)))/(abs(s(i-1,2))+abs(s(i,2)));
            break;
        end
    end
    [t,s]=ode45(@compute_State_and_Phi_Derivative,[0 t_end],initialcondition,odeset('RelTol',1e-12,'AbsTol',1e-12));
    plot(s(:,1),s(:,2),'linewidth',2,'displayname',[numbering(q,:),' try']); hold on;
    % Update guess:
    current_Derivative = compute_State_and_Phi_Derivative(0,s(end,:)'); % state derivative at t1/2
    Phi34 = s(end,19); % Phi_34 at t1/2
    Phi24 = s(end,18); % Phi_24 at t1/2
    xdotdot = current_Derivative(3); % xdotdot at t1/2
    ydot = s(end,4); % ydot at t1/2
    xdot = s(end,3); % xdot at t1/2
    dydot = -(Phi34-xdotdot/ydot*Phi24)^(-1)*xdot;
    if(abs(xdot)<1e-8)
        break;
    end
    initial_State_guess = initial_State_guess+[0;0;0;dydot];
end
T3 = 2*t_end; % final period L3
initial3 = initial_State_guess; % final initial guess L3
plot(x_L3,0,'kx','linewidth',2,'markersize',10,'displayname','L3');
axis equal; grid on;
legend('interpreter','latex');
xlabel('$x$','interpreter','latex');
ylabel('$y$','interpreter','latex');

%% Task 3.2 and Task 3.3a

figure;
% L1
initial_Phi = eye(4);
initialcondition = [initial1;initial_Phi(:,1);initial_Phi(:,2);initial_Phi(:,3);initial_Phi(:,4)];
[nominalt,nominalorbit]=ode45(@compute_State_and_Phi_Derivative,[0 T1],initialcondition,odeset('RelTol',1e-12,'AbsTol',1e-12)); % Orbit to determine 10 starting conditions
initialconditions1 = cell(10,1);
for o=0:9
    start_time = o/10*T1;
    start_state = interpolate_State_at_Time(nominalt,nominalorbit,start_time)';
    initialconditions1{o+1} = start_state;
    initialcondition = [start_state(1:4);initial_Phi(:,1);initial_Phi(:,2);initial_Phi(:,3);initial_Phi(:,4)];
    [t,state]=ode45(@compute_State_and_Phi_Derivative,[0 T1],initialcondition,odeset('RelTol',1e-12,'AbsTol',1e-12));
    M = [state(end,5:8)',state(end,9:12)',state(end,13:16)',state(end,17:20)'];
    [vec,eigval] = eig(M);
    if(o==0)
        eigenvalues1 = diag(eigval); % for report
    end
    [~,indexmax] = max(real(diag(eigval)));
    [~,indexmin] = min(real(diag(eigval)));
    start_state_perturbed = start_state+1e-5*[vec(:,indexmax);zeros(16,1)];
    initialcondition = [start_state_perturbed(1:4);initial_Phi(:,1);initial_Phi(:,2);initial_Phi(:,3);initial_Phi(:,4)];
    [t,state]=ode45(@compute_State_and_Phi_Derivative,[0 maxtime],initialcondition,odeset('RelTol',1e-12,'AbsTol',1e-12));
    if(o==0)
        plot(state(:,1),state(:,2),'linewidth',2,'color',[0 166 214]/255,'displayname','unstable manifold'); hold on;
    else
        plot(state(:,1),state(:,2),'linewidth',2,'color',[0 166 214]/255,'handlevisibility','off'); hold on;
    end
    start_state_perturbed = start_state-1e-5*[vec(:,indexmax);zeros(16,1)];
    initialcondition = [start_state_perturbed(1:4);initial_Phi(:,1);initial_Phi(:,2);initial_Phi(:,3);initial_Phi(:,4)];
    [t,state]=ode45(@compute_State_and_Phi_Derivative,[0 maxtime],initialcondition,odeset('RelTol',1e-12,'AbsTol',1e-12));
    plot(state(:,1),state(:,2),'linewidth',2,'color',[0 166 214]/255,'handlevisibility','off'); hold on;
    start_state_perturbed = start_state+1e-5*[vec(:,indexmin);zeros(16,1)];
    initialcondition = [start_state_perturbed(1:4);initial_Phi(:,1);initial_Phi(:,2);initial_Phi(:,3);initial_Phi(:,4)];
    [t,state]=ode45(@compute_State_and_Phi_Derivative,[0 -maxtime],initialcondition,odeset('RelTol',1e-12,'AbsTol',1e-12));
    plot(state(:,1),state(:,2),'linewidth',2,'color',[230 70 22]/255,'handlevisibility','off'); hold on;
    start_state_perturbed = start_state-1e-5*[vec(:,indexmin);zeros(16,1)];
    initialcondition = [start_state_perturbed(1:4);initial_Phi(:,1);initial_Phi(:,2);initial_Phi(:,3);initial_Phi(:,4)];
    [t,state]=ode45(@compute_State_and_Phi_Derivative,[0 -maxtime],initialcondition,odeset('RelTol',1e-12,'AbsTol',1e-12));
    if(o==0)
        plot(state(:,1),state(:,2),'linewidth',2,'color',[230 70 22]/255,'displayname','stable manifold'); hold on;
    else
        plot(state(:,1),state(:,2),'linewidth',2,'color',[230 70 22]/255,'handlevisibility','off'); hold on;
    end
end
plot(nominalorbit(:,1),nominalorbit(:,2),'k','linewidth',2,'displayname','Lyapunov orbit'); hold on;
plot(-mu,0,'kp','linewidth',3,'markersize',10,'displayname','Sun');
plot(1-mu,0,'k.','linewidth',2,'markersize',30,'displayname','Mars');
plot(x_L1,0,'k+','linewidth',2,'markersize',10,'displayname','L1');
plot(x_L3,0,'kx','linewidth',2,'markersize',10,'displayname','L3');
axis equal; grid on;
legend('interpreter','latex');
xlabel('$x$','interpreter','latex');
ylabel('$y$','interpreter','latex');

% L3 (without figure)
initial_Phi = eye(4);
initialcondition = [initial3;initial_Phi(:,1);initial_Phi(:,2);initial_Phi(:,3);initial_Phi(:,4)];
[nominalt,nominalorbit]=ode45(@compute_State_and_Phi_Derivative,[0 T3],initialcondition,odeset('RelTol',1e-12,'AbsTol',1e-12)); % Orbit to determine 10 starting conditions
initialconditions3 = cell(10,1);
for o=0:9
    start_time = o/10*T3;
    start_state = interpolate_State_at_Time(nominalt,nominalorbit,start_time)';
    initialconditions3{o+1} = start_state;
    initialcondition = [start_state(1:4);initial_Phi(:,1);initial_Phi(:,2);initial_Phi(:,3);initial_Phi(:,4)];
    [t,state]=ode45(@compute_State_and_Phi_Derivative,[0 T3],initialcondition,odeset('RelTol',1e-12,'AbsTol',1e-12));
    M = [state(end,5:8)',state(end,9:12)',state(end,13:16)',state(end,17:20)'];
    [vec,eigval] = eig(M);
    if(o==0)
        eigenvalues3 = diag(eigval); % for report
    end
    [~,indexmax] = max(real(diag(eigval)));
    [~,indexmin] = min(real(diag(eigval)));
end

%% Task 3.3b

maxtime = 5/t_unity;

global delta; delta = 0; % clock angle
global alpha; % cone angle
alpha_deg = -70:2.5:70;
least_norm = 3*ones(size(alpha_deg));
best_position_error = zeros(size(alpha_deg));
best_velocity_error = zeros(size(alpha_deg));
time1 = zeros(size(alpha_deg));
time3 = zeros(size(alpha_deg));
best_o1 = zeros(size(alpha_deg));
best_o3 = zeros(size(alpha_deg));
for alpha_index=1:length(alpha_deg)
    alpha = alpha_deg(alpha_index)*pi/180;
    times1 = cell(10,1);
    states1 = cell(10,1);
    times3 = cell(10,1);
    states3 = cell(10,1);
    for o = 1:10 % iterates over the ten places at the Lyapunov orbits
        initialcondition = [initialconditions1{o}(1:2);0;initialconditions1{o}(3:4);0];
        [t1,s1]=ode45(@compute_State_Derivative_with_Sailing,[0 maxtime],initialcondition,odeset('RelTol',1e-12,'AbsTol',1e-12));
        times1{o} = t1;
        states1{o} = s1;
        initialcondition = [initialconditions3{o}(1:2);0;initialconditions3{o}(3:4);0];
        [t3,s3]=ode45(@compute_State_Derivative_with_Sailing,[0 -maxtime],initialcondition,odeset('RelTol',1e-12,'AbsTol',1e-12));
        times3{o} = t3;
        states3{o} = s3;
    end
    for o1=1:10
        for o3=1:10
            i = 1;
            while i<=length(times1{o1})
                j = 1;
                while j<=length(times3{o3})
                    if(norm(states1{o1}(i,:)-states3{o3}(j,:))<least_norm(alpha_index))
                        least_norm(alpha_index) = norm(states1{o1}(i,:)-states3{o3}(j,:));
                        best_position_error(alpha_index) = norm(states1{o1}(i,1:3)-states3{o3}(j,1:3));
                        best_velocity_error(alpha_index) = norm(states1{o1}(i,4:6)-states3{o3}(j,4:6));
                        time1(alpha_index) = times1{o1}(i);
                        time3(alpha_index) = times3{o3}(j);
                        best_o1(alpha_index) = o1;
                        best_o3(alpha_index) = o3;
                        j = j+1;
                    elseif(norm(states1{o1}(i,:)-states3{o3}(j,:))<2*least_norm)
                        j = j+1;
                    elseif(norm(states1{o1}(i,:)-states3{o3}(j,:))<10*least_norm)
                        j = j+10;
                    else
                        j = j+100;
                    end
                end
                i = i+1;
            end
        end
    end
end
[~,best_alpha_index] = min(least_norm);

%%

fig = figure;
set(fig,'defaultAxesColorOrder',[[226 26 26]/255; [0 166 214]/255]);
yyaxis left;
semilogy(alpha_deg,least_norm,'x','linewidth',2,'color',[226 26 26]/255,'handlevisibility','off'); hold on; grid on;
xlabel('cone angle $\alpha$','interpreter','latex');
ylabel('minimum Euclidean-norm error','interpreter','latex');
yyaxis right;
plot(alpha_deg,best_o1,'x','linewidth',2,'color',[0 166 214]/255,'displayname','L1'); hold on; grid on;
plot(alpha_deg,best_o3,'o','linewidth',2,'color',[0 166 214]/255,'displayname','L3'); hold on; grid on;
ylabel('Index of point along Lyapunov orbit','interpreter','latex');
xlim([-70 70]);
xline(alpha_deg(best_alpha_index),'handlevisibility','off');
xline(-alpha_deg(best_alpha_index),'displayname','best $\alpha$');
legend('interpreter','latex');

figure;
alpha = alpha_deg(best_alpha_index)*pi/180;
initialcondition = [initialconditions1{best_o1(best_alpha_index)}(1:2);0;initialconditions1{best_o1(best_alpha_index)}(3:4);0];
[t1,s1]=ode45(@compute_State_Derivative_with_Sailing,[0 time1(best_alpha_index)],initialcondition,odeset('RelTol',1e-12,'AbsTol',1e-12));
plot(s1(:,1),s1(:,2),'linewidth',2,'color',[0 166 214]/255,'displayname','unstable manifold L1'); hold on;
initialcondition = [initialconditions3{best_o3(best_alpha_index)}(1:2);0;initialconditions3{best_o3(best_alpha_index)}(3:4);0];
[t3,s3]=ode45(@compute_State_Derivative_with_Sailing,[0 time3(best_alpha_index)],initialcondition,odeset('RelTol',1e-12,'AbsTol',1e-12));
plot(s3(:,1),s3(:,2),'linewidth',2,'color',[230 70 22]/255,'displayname','stable manifold L3'); hold on;
[~,Lyapunov]=ode45(@compute_State_Derivative,[0 T1],[initial1(1:2);0;initial1(3:4);0],odeset('RelTol',1e-12,'AbsTol',1e-12));
plot(Lyapunov(:,1),Lyapunov(:,2),'k','linewidth',2,'displayname','Lyapunov orbit L1'); hold on;
[~,Lyapunov]=ode45(@compute_State_Derivative,[0 T3],[initial3(1:2);0;initial3(3:4);0],odeset('RelTol',1e-12,'AbsTol',1e-12));
plot(Lyapunov(:,1),Lyapunov(:,2),'k','linewidth',2,'displayname','Lyapunov orbit L3'); hold on;
plot(-mu,0,'kp','linewidth',3,'markersize',10,'displayname','Sun');
plot(1-mu,0,'k.','linewidth',2,'markersize',30,'displayname','Mars');
plot(x_L1,0,'k+','linewidth',2,'markersize',10,'displayname','L1');
plot(x_L3,0,'kx','linewidth',2,'markersize',10,'displayname','L3');
axis equal; grid on;
legend('interpreter','latex');
xlabel('$x$','interpreter','latex');
ylabel('$y$','interpreter','latex');
position_error = best_position_error(best_alpha_index)*a; % [km]
velocity_error = best_velocity_error(best_alpha_index)*a/t_units; % [km/s]
time_unstable = time1(best_alpha_index)*t_unity; % [Earth years]
time_stable = -time3(best_alpha_index)*t_unity; % [Earth years]

%% Functions

% First derivative of the potential U wrt x
function Ux = compute_Ux(x,y,z)
global mu;
[r1,r2] = get_r1r2(x,y,z);
Ux = -x+(1-mu)*(x+mu)/r1^3+mu*(x-1+mu)/r2^3;
end

% Second derivative of the potential U wrt x
function Uxx = compute_Uxx(x,y,z)
global mu;
[r1,r2] = get_r1r2(x,y,z);
Uxx = -1+(1-mu)*(-3*(x+mu)^2/r1^5+1/r1^3)+mu*(-3*(x-1+mu)^2/r2^5+1/r2^3);
end

% First derivative of the potential U wrt y
function Uy = compute_Uy(x,y,z)
global mu;
[r1,r2] = get_r1r2(x,y,z);
Uy = -y+(1-mu)*y/r1^3+mu*y/r2^3;
end

% Second derivative of the potential U wrt y
function Uyy = compute_Uyy(x,y,z)
global mu;
[r1,r2] = get_r1r2(x,y,z);
Uyy = -1+(1-mu)*(-3*y^2/r1^5+1/r1^3)+mu*(-3*y^2/r2^5+1/r2^3);
end

% Mixed derivative of the potential U wrt x and y
function Uxy = compute_Uxy(x,y,z)
global mu;
[r1,r2] = get_r1r2(x,y,z);
Uxy = -3*y*((1-mu)*(x+mu)/r1^5+mu*(x-1+mu)/r2^5);
end

% for 6 dimensional state, without sailing, without state transition matrix
function State_Derivative = compute_State_Derivative(t,State)
x = State(1);
y = State(2);
z = State(3);
potential_acceleration = [compute_Ux(x,y,z);compute_Uy(x,y,z);compute_Uz(x,y,z)];
coriolis_acceleration = -2*cross([0;0;1],State(4:6));
State_Derivative = [State(4:6);-potential_acceleration+coriolis_acceleration];
end

% for 4 dimensional state, without sailing, with state transition matrix.
% State = [x;y;xdot;ydot;Phi(:,1);Phi(:,2);Phi(:,3);Phi(:,4)]
function Derivative = compute_State_and_Phi_Derivative(t,State)
x = State(1);
y = State(2);
xdot = State(3);
ydot = State(4);
potential_acceleration = [compute_Ux(x,y,0);compute_Uy(x,y,0)];
coriolis_acceleration = -2*cross([0;0;1],[State(3);State(4);0]);
coriolis_acceleration = coriolis_acceleration(1:2);
State_Derivative = [xdot;ydot;-potential_acceleration+coriolis_acceleration];

A = [0,0,1,0;0,0,0,1;-compute_Uxx(x,y,0),-compute_Uxy(x,y,0),0,2;-compute_Uxy(x,y,0),-compute_Uyy(x,y,0),-2,0];
Phi = [State(5:8),State(9:12),State(13:16),State(17:20)];
Phi_Derivative = A*Phi;
Derivative = [State_Derivative;Phi_Derivative(:,1);Phi_Derivative(:,2);Phi_Derivative(:,3);Phi_Derivative(:,4)];
end

% for 6 dimensional state, with sailing, without state transition matrix
function State_Derivative = compute_State_Derivative_with_Sailing(t,State)
global beta mu delta alpha;
x = State(1);
y = State(2);
z = State(3);
[r1,~] = get_r1r2(x,y,z);
r1_norm = [x+mu;y;z]/r1;
theta_norm = cross([0;0;1],r1_norm)/norm(cross([0;0;1],r1_norm));
eta_norm = cross(r1_norm,theta_norm);
n_norm = [r1_norm,theta_norm,eta_norm]*[cos(alpha);sin(alpha)*sin(delta);sin(alpha)*cos(delta)];
potential_acceleration = [compute_Ux(x,y,z);compute_Uy(x,y,z);compute_Uz(x,y,z)];
coriolis_acceleration = -2*cross([0;0;1],State(4:6));
sailing_acceleration = beta*(1-mu)/r1^2*dot(r1_norm,n_norm)^2*n_norm;
State_Derivative = [State(4:6);-potential_acceleration+coriolis_acceleration+sailing_acceleration];
end

% First derivative of the potential U wrt z
function Uz = compute_Uz(x,y,z)
global mu;
[r1,r2] = get_r1r2(x,y,z);
Uz = (1-mu)*z/r1^3+mu*z/r2^3;
end

function [r1,r2] = get_r1r2(x,y,z)
global mu;
r1 = sqrt((x+mu)^2+y^2+z^2);
r2 = sqrt((x-1+mu)^2+y^2+z^2);
end

% linear interpolation
function state = interpolate_State_at_Time(times_vec,states_mat,time)
if(time<=times_vec(1))
    state = states_mat(1,:);
else
    i = 1;
    while(time>times_vec(i))
        i = i+1;
    end
    state = (states_mat(i,:)*abs(times_vec(i-1)-time)+states_mat(i-1,:)*abs(times_vec(i)-time))/(abs(times_vec(i-1)-time)+abs(times_vec(i)-time));
end
end