% Piecewise Jerk speed Problem
clc;
clear;
%% Parameters
tic
global Nstep; % Num of Points
global s_tol; % S Length
global Nstep3;
global ref_v;
global delta_s;

global v_ini 

global Maxspeed Accmax Jerklimit;


% QP Problem settings
global w1 w2 w3;

w1 = 1;
w2 = 0.1;
w3 = 1e2;


Maxspeed = 30;
Accmax = 4;
Jerklimit = 0.5*delta_s;


global st_seg_num % num of st curves 1 2 3
st_seg_num = 3;

global ST_C st_seg
ST_C = [0 10 0 0 0 0 -25 15 0  0 0 0 -15 14 0  0 0 0];
st_seg = [50 5;
    125 10;
    195 15];


s_tol = st_seg(3,2);

ref_v = ST_C(14); % ref velocity;

delta_s = 0.5;

Nstep = s_tol/delta_s + 1;

Nstep3 = 3*Nstep;


global s_ini  a_ini s_end v_end a_end
s_ini = 0;
v_ini = ST_C(2);
a_ini = 0;
s_end = st_seg(3,1);
v_end = ST_C(14); % target velocity
a_end = 0;

global veh_front veh_behind
veh_front = [80 10];
veh_behind = [-60 15];

%% Compute l_bound up and down
l_upbound = zeros(Nstep,1);
l_downbound = zeros(Nstep,1);

s_sample = linspace(0,s_tol,Nstep);

for i = 1:Nstep
    l_upbound(i) = veh_front(1) + veh_front(2)*s_sample(i);
    l_downbound(i) = veh_behind(1) + veh_behind(2)*s_sample(i);
end



%% Get ref_x and ref_v

ref_x = zeros(Nstep,1);

for i = 1:Nstep
    if s_sample(i) <= st_seg(1,2)
        ref_x(i) = ST_C(1) + ST_C(2)*s_sample(i);
    elseif s_sample(i) <= st_seg(2,2)
        ref_x(i) = ST_C(7) + ST_C(8)*s_sample(i);
    else
        ref_x(i) = ST_C(13) + ST_C(14)*s_sample(i);
    end
end

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
    ul(i) = l_upbound(i);
    ul(i+4*Nstep-1) = -l_downbound(i);
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
toc
opts = optimoptions('quadprog','Algorithm','active-set','Display','iter','Maxiterations',200);
%opts = optimoptions(opts,'Maxiterations',200);

[x,fval,exitflag,output,lambda] = quadprog(P,c,A,ul,Aeq,beq,[],[],x0,opts);

vel_o = zeros(Nstep,1);
vel_o(1) = v_ini;
acc_o = zeros(Nstep,1);
acc_o(1) = a_ini;
jerk_o = zeros(Nstep,1);
jerk = zeros(Nstep,1);
for i=1:Nstep-1
    vel_o(i+1) = (ref_x(i+1) - ref_x(i))/delta_s;
    if abs(vel_o(i+1) - vel_o(i)) > 0.0001
        acc_o(i) = (vel_o(i+1) - vel_o(i))/delta_s;
    else
        acc_o(i) = 0;
    end
    if abs(acc_o(i+1) - acc_o(i)) > 0.0001
        jerk_o(i) = (acc_o(i+1) - acc_o(i))/delta_s;
    else
        jerk_o(i) = 0;
    end
    jerk(i) = (x(2*Nstep+i+1) - x(2*Nstep+i))/delta_s;
end

figure();
subplot(3,1,1);
plot(s_sample,ref_x);
hold on;
plot(s_sample,x(1:Nstep),'*');
hold on;
plot(s_sample,l_upbound);
hold on;
plot(s_sample,l_downbound);
legend('original','op','upbound','lowbound','Location','NorthWest');
title('s-t')
subplot(3,1,2);
plot(s_sample,acc_o);
hold on;
plot(s_sample,x(2*Nstep+1:3*Nstep),'+');
legend('original','op','Location','NorthWest');
title('Acc')

subplot(3,1,3);
plot(s_sample,jerk_o);
hold on;
plot(s_sample,jerk,'*');
title('Jerk')
legend('original','op','Location','NorthWest');
% figure();
% plot(s_sample,x(2*Nstep+1:3*Nstep),'+');


