% Piecewise Jerk Path Problem
clc;
clear;
%% Parameters

global Nstep; % Num of Points
global s_tol; % S Length
global Nstep3;
global ref_v;
global delta_s;

global v_ini;a

global Maxspeed Accmax Jerklimit;

global safe_zone;

% QP Problem settings
global w1 w2 w3;


Nstep = 51;
s_tol = 20;
Nstep3 = 3*Nstep;
delta_s = s_tol/(Nstep-1);

v_ini = 0.01;

w1 = 20;
w2 = 10;
w3 = 2;

ref_v = 0; % ref velocity;

Maxspeed = 1.5;
Accmax = 2;
Jerklimit =5*delta_s;

safe_zone = 0.5;

%% Compute l_bound up and down
l_upbound = zeros(Nstep,1);
l_downbound = zeros(Nstep,1);

s_sample = linspace(0,s_tol,Nstep);

% generate bounds

% for i = 1:Nstep
%     l_upbound(i) = -sqrt(35^2-s_sample(i)^2) + 35;
%     l_downbound(i) =-sqrt(39^2 - s_sample(i)^2) + 35;
% end

for i = 1:Nstep
    if i < Nstep/3
        l_upbound(i) = 5;
    elseif i<2*Nstep/3 && i >= Nstep/3
        l_upbound(i) = 2;
    else
        l_upbound(i) = 6;
    end
end

for i = 1:Nstep
    if i < Nstep/3-2
        l_downbound(i) = 1;
    elseif i<2*Nstep/3+2 && i >= Nstep/3-2
        l_downbound(i) = -2;
    else
        l_downbound(i) = 3;
    end
end


%% Get ref_x and ref_v

ref_x = zeros(Nstep,1);
ref_x = 0.3*l_upbound + 0.7*l_downbound;

x0 = zeros(Nstep3,1);
x0(1:Nstep) = ref_x;

%% Constrcuct QP constraints Matrix
% P matrix in dense format 
P = zeros(Nstep3,Nstep3);
% cost vector c matrix
c = zeros(Nstep3,1);

for i = 1:Nstep
    P(i,i) = 2*w1;
    c(i) = -2*ref_x(i)*w1;
end

for i = Nstep+1:2*Nstep
    P(i,i) = 2*w2;
    c(i) = -2*ref_v*w2;
end

for i = 2*Nstep+1:Nstep3
    P(i,i) = 2*w3;
    %c(i) = 0;
end

% Aeq
Aeq = zeros(2*Nstep+1,Nstep3);
% Part 1
for i = 1:Nstep-1
    
    j = i + Nstep;
    Aeq(i,j) = -1;
    Aeq(i,j+1) = 1;
    
    jj = i + 2*Nstep;
    Aeq(i,jj) = -delta_s/2.0;
    Aeq(i,jj+1) = -delta_s/2.0;
end
% Part 2
for i = Nstep:2*Nstep-2
    
    j = i - Nstep + 1;
    Aeq(i,j) = -1;
    Aeq(i,j+1) = 1;
    
    jj = i + 1;
    Aeq(i,jj) = -delta_s;
    
    jjj = i+Nstep+1;
    Aeq(i,jjj) = -delta_s*delta_s/3.0;
    Aeq(i,jjj+1) = -delta_s*delta_s/6.0;
  
end
%Part 3
Aeq(2*Nstep-1,1) = 1;
Aeq(2*Nstep,Nstep+1) = 1;
Aeq(2*Nstep+1,2*Nstep+1) = 1;

% Beq
beq = zeros(2*Nstep+1,1);

beq(2*Nstep-1) = ref_x(1);
beq(2*Nstep) = v_ini;
beq(2*Nstep+1) = 0;

% A

A = zeros(8*Nstep-2,Nstep3);
% A
for i = 1:Nstep3
    A(i,i) = 1;
end

for i = Nstep3+1:4*Nstep-1
    j = i - Nstep;
    A(i,j) = -1;
    A(i,j+1) = 1;
end

% -A
for i = 4*Nstep:7*Nstep-1
    j = i - 4*Nstep + 1;
    A(i,j) = -1;
end

for i = 7*Nstep:8*Nstep-2
    j = i - 5*Nstep + 1;
    A(i,j) = 1;
    A(i,j+1) = -1;
end


% ul

ul = zeros(8*Nstep-2,1); % lb < A < ub both in ul!
for i = 1:Nstep
    ul(i) = l_upbound(i)+safe_zone;
    ul(i+4*Nstep-1) = -l_downbound(i)-safe_zone;
end

for i =Nstep+1:2*Nstep
    ul(i) = Maxspeed;
    ul(i+4*Nstep-1) = Maxspeed;
end

for i = 2*Nstep+1:3*Nstep
    ul(i) = Accmax;
    ul(i+4*Nstep-1) = Accmax; 
end

for i = 3*Nstep+1:4*Nstep-1
    ul(i) = Jerklimit;
    ul(i+4*Nstep-1) = Jerklimit;
end

opts = optimoptions('quadprog','Algorithm','active-set','Display','iter','Maxiterations',200);
%opts = optimoptions(opts,'Maxiterations',200);
[x,fval,exitflag,output,lambda] = quadprog(P,c,A,ul,Aeq,beq,[],[],x0,opts);




plot(s_sample,l_upbound);
hold on;
plot(s_sample,l_downbound);
hold on;
plot(s_sample,ref_x,'*');
hold on;
plot(s_sample,x(1:Nstep));



