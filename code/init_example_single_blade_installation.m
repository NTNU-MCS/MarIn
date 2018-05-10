clc
clear

%% Wind field
meanWindSpeed = 10;
load(['WSP10_ClassC_Seed48_1000s.mat']);

%% Coordinate
load('bladePar.mat');
g = 9.81;
installation_height = -90;
m_yoke = 20e3;

%% Crane tip
posCraneTip = [0;0;-110];

%% Lift wire and sling stiffness
E = 2.1e11;
D = 0.25;
A = pi/4*D^2;
EA = E*A;

gamma = 0.1;

lw1_init = 10;
k_w1 = gamma * EA/lw1_init;

liftSpeed = 0;
%% Hock
m_h = 1e3;

elong_mainlw = (bladePar.totalMass+m_h+m_yoke)*g/k_w1;

p_h_init = posCraneTip + [0;0;lw1_init + elong_mainlw];
p_h_dot_init = [0;0;0];

%% Blade
p_init_n = [0;0;-90];
Theta_init_n = [0;-90*pi/180;0];

eta_init_n = [p_init_n;Theta_init_n];
nu_init_b = zeros(6,1);

R_b2n = Rzxy(Theta_init_n(1),Theta_init_n(2),Theta_init_n(3));

%% Slings
posSling_b = [bladePar.posCOG + [-2;-4.5;0] , bladePar.posCOG + [-2;4.5;0]];

for i =1:length(posSling_b(1,:))
    posSling_n(:,i) = eta_init_n(1:3)+R_b2n*(posSling_b(:,i)-bladePar.posCOG);
end

k_w2 = 1e8;
k_w3 = k_w2;

elong_lw2 = (bladePar.totalMass+m_yoke)*g/(sqrt(3)*k_w2);
elong_lw3 = elong_lw2;

lw2 = norm(p_h_init- posSling_n(:,1))-elong_lw2;
lw3 = norm(p_h_init- posSling_n(:,2))-elong_lw3;

%% Taglines
posTuggerline_b = [bladePar.posCOG+[0;-4.5;0],bladePar.posCOG+[0;4.5;0]];

for i =1:length(posTuggerline_b(1,:))
    posTuggerline_n(:,i) = eta_init_n(1:3)+R_b2n*(posTuggerline_b(:,i)-bladePar.posCOG);
    posTuggerlineBase_n(:,i) = posTuggerline_n(:,i) - [10;0;0];
end

lt1_init = norm(posTuggerline_n(:,1)- posTuggerlineBase_n(:,1));
lt2_init = norm(posTuggerline_n(:,2)- posTuggerlineBase_n(:,2));

k_t1 = 1e7;
k_t2 = k_t1;

for i =1:length(posSling_b(1,:))
    pLiftwire_init_n(:,i) = eta_init_n(1:3)+R_b2n*(posSling_b(:,i)-bladePar.posCOG);
end

d_w1 = k_w1/10e3;
d_w2 = k_w2/10e3;
d_w3 = d_w2;
d_t1 = 0; %k_t1/1e4;
d_t2 = d_t1;

%% Simulation
tend = 1000;
stepTime = 100;