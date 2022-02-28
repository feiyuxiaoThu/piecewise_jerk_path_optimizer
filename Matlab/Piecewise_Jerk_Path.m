% Piecewise Jerk Path Problem
clc;
clear all;
%% Parameters

global Nstep; % Num of Points
global s_tol; % S Length
global Nstep3;
global ref_v;
global delta_s;

global v_ini;

global Maxspeed Accmax Jerklimit;

% QP Problem settings
global w1 w2 w3;


Nstep = 11;
s_tol = 30;
Nstep3 = 3*Nstep;
delta_s = s_tol/(Nstep-1);

v_ini = 15;

w1 = 20;
w2 = 10;
w3 = 2;

ref_v = 20; % ref velocity;

Maxspeed = 33.3;
Accmax = 5;
Jerklimit =50*delta_s;

%% Compute l_bound up and down
l_upbound = zeros(Nstep,1);
l_downbound = zeros(Nstep,1);

s_sample = linspace(0,s_tol,Nstep);

for i = 1:Nstep
    l_upbound(i) = -sqrt(35^2-s_sample(i)^2) + 35;
    l_downbound(i) = -sqrt(39^2 - s_sample(i)^2) + 35;
end


%% Get ref_x and ref_v

ref_x = zeros(Nstep,1);
ref_x = 0.3*l_upbound + 0.7*l_downbound;

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
Aeq(2*Nstep-1,1) = 1.0;
Aeq(2*Nstep,Nstep+1) = 1.0;
Aeq(2*Nstep+1,2*Nstep+1) =1.0;

% Beq
beq = zeros(2*Nstep+1,1);

beq(2*Nstep) = v_ini;

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
    ul(i) = l_upbound(i);
    ul(i+4*Nstep-1) = -l_downbound(i);
end

for i =Nstep+1:2*Nstep
    ul(i) = Maxspeed;
end

for i = 2*Nstep+1:3*Nstep
    ul(i) = Accmax;
    ul(i+4*Nstep-1) = Accmax; 
end

for i = 3*Nstep+1:4*Nstep-1
    ul(i) = Jerklimit;
    ul(i+4*Nstep-1) = Jerklimit;
end

[x,fval,exitflag,output,lambda] = quadprog(P,c,A,ul,Aeq,beq)





plot(s_sample,l_upbound);
hold on;
plot(s_sample,l_downbound);
hold on;
plot(s_sample,ref_x,'*');



