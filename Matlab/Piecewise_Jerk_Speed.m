% ST QP
%% Parameters
clear;
clc;
global st_seg_num % num of st curves 1 2 3
st_seg_num = 2;

% QP Weights
global w1 w2 w3
w1 = 0.0005; % error with DP ST results
w2 = 1; % acc
w3 = 1e3; %  jerk

global ST_Results
ST_Results = [0 10 0 0 0 0 50 5; 
    -25 15 0  0 0 0 125 10];

global s_ini v_ini a_ini s_end v_end a_end
s_ini = 0;
v_ini = ST_Results(1,2);
a_ini = 0;
s_end = ST_Results(2,7);
v_end = 15; % target velocity
a_end = 0;

global veh_front veh_behind
veh_front = [40 10];
veh_behind = [-80 15];

% PLot
t_sample = 0:1:10;
s_front = veh_front(1) + veh_front(2)*t_sample;
s_behind = veh_behind(1) + veh_behind(2)*t_sample;



%% Target function
H0 = zeros(6,6,st_seg_num); % distance (s-s_dp)^2
H1 = zeros(6,6,st_seg_num); % acceleration
H2 = zeros(6,6,st_seg_num); % jerk

f1 = zeros(6,st_seg_num); % part 1 of f
f2 = zeros(6,st_seg_num);

for i = 1:st_seg_num
    T = ST_Results(i,8); % T_seg
    H0(:,:,i) = 2*[T T^2/2 T^3/3 T^4/4 T^5/5 T^6/6;
        T^2/2 T^3/3 T^4/4 T^5/5 T^6/6 T^7/7;
        T^3/3 T^4/4 T^5/5 T^6/6 T^7/7 T^8/8;
        T^4/4 T^5/5 T^6/6 T^7/7 T^8/8 T^9/9;
        T^5/5 T^6/6 T^7/7 T^8/8 T^9/9 T^10/10;
        T^6/6 T^7/7 T^8/8 T^9/9 T^10/10 T^11/11];
    
    H1(:,:,i) = 2*[0 0 0 0 0 0;
        0 0 0 0 0 0;
        0 0  4*T 6*T^2 8*T^3 10*T^4;
        0 0 6*T^2 12*T^3 18*T^4 24*T^5;
        0 0 8*T^3 18*T^4 144*T^5/5 40*T^6;
        0 0 10*T^4 24*T^5 40*T^6 400*T^7/7]; % acc
    
    H2(:,:,i) = 2* [0,0,0,0,0,0;
        0,0,0,0,0,0;
        0,0,0,0,0,0;
        0,0,0,36*T,72*T^2,120*T^3;
        0,0,0,72*T^2,192*T^3,360*T^4;
        0,0,0,120*T^3,360*T^4,720*T^5];
    
    f1(:,i) = [T,T^2/2,T^3/3,T^4/4,T^5/5,T^6/6];
    f2(:,i) = [T^2/2,T^3/3,T^4/4,T^5/5,T^6/6,T^7/7];
end


H_seg1 = w1*H0(:,:,1) + w2*H1(:,:,1) + w3*H2(:,:,1);
H_seg2 = w1*(H0(:,:,2)-H0(:,:,1)) + w2*(H1(:,:,2)-H1(:,:,1)) + w3*(H2(:,:,2)-H2(:,:,1));

H_final = blkdiag(H_seg1,H_seg2);

f_seg1 = ST_Results(1,1)*f1(:,1) + ST_Results(1,2)*f2(:,1);
f_seg2 = ST_Results(2,1)*(f1(:,2)-f1(:,1)) + ST_Results(2,2)*(f2(:,2)-f2(:,1));

f_final = -2*w1*vertcat(f_seg1,f_seg2);


%% constraints
% two segs for example

t_seg = ST_Results(1,8);

% segs equal
Aeq1 = [1 t_seg t_seg^2 t_seg^3 t_seg^4 t_seg^5 -1 -t_seg -t_seg^2 -t_seg^3 -t_seg^4 -t_seg^5]; % s_left = s_right
Aeq2 = [0 1 2*t_seg 3*t_seg^2 4*t_seg^3 5*t_seg^4 0 -1 -2*t_seg -3*t_seg^2 -4*t_seg^3 -5*t_seg^4]; % v_left = v_right
Aeq3 = [0 0 2 6*t_seg 12*t_seg^2 20*t_seg^3 0 0 -2 -6*t_seg -12*t_seg^2 -20*t_seg^3];


Aseg = vertcat(Aeq1,Aeq2,Aeq3);
bseg = [0;0;0];

% initial and final conditions

%t_seg = 10; %ST_Results(2,8); % final time

t_seg1 = ST_Results(2,8);
Aini = [1 0 0 0 0 0;
    0 1 0 0 0 0;
    0 0 2 0 0 0];
Aend = [
    0 1 2*t_seg1 3*t_seg1^2 4*t_seg1^3 5*t_seg1^4;
    0 0 2 6*t_seg1 12*t_seg1^2 20*t_seg1^3]; %1 t_seg1 t_seg1^2 t_seg1^3 t_seg1^4 t_seg1^5;



%A_time = horzcat(Aini,zeros(3,6));
A_time = blkdiag(Aini,Aend);
b_time = [s_ini;v_ini;a_ini;v_end;a_end];

Aeq = vertcat(Aseg,A_time);
beq = vertcat(bseg,b_time);

% upper and lower bounds for s
% time_sample
A_sample = zeros(length(t_sample),12);
A_sample_speed = zeros(length(t_sample),12); % speed > 0 ds > 0
for i = 1:length(t_sample)
    T = t_sample(i);
    if t_sample(i) < ST_Results(1,8)
        A_sample(i,:) = [1 T T^2 T^3 T^4 T^5 0 0 0 0 0 0];
        A_sample_speed(i,:) = [0 1 2*T 3*T^2 4*T^3 5*T^4 0 0 0 0 0 0];
    else
        A_sample(i,:) = [0 0 0 0 0 0 1 T T^2 T^3 T^4 T^5];
        A_sample_speed(i,:) = [0 0 0 0 0 0 0 1 2*T 3*T^2 4*T^3 5*T^4];
    end
end

A = vertcat(A_sample, -A_sample);

ub = vertcat(s_front',-s_behind');


options = optimoptions('quadprog','Display','iter');
C_new = quadprog(H_final,f_final,A,ub,Aeq,beq,[],[],[],options);

s_ego = zeros(1,length(t_sample));
for i = 1:length(t_sample)
    if t_sample(i) <= 5 %ST_Results(1,8)
        s_ego(i) = ST_Results(1,1) + ST_Results(1,2)*t_sample(i);
    else
        s_ego(i) = ST_Results(2,1) + ST_Results(2,2)*t_sample(i);
    end
end

s_ego_new = zeros(1,length(t_sample));
for i = 1:length(t_sample)
    T = t_sample(i);
    if T <= 5 %ST_Results(1,8)
        s_ego_new(i) = C_new(1) + C_new(2)*T + C_new(3)*T^2 + C_new(4)*T^3 + C_new(5)*T^4 + C_new(6)*T^6;
    else
        s_ego_new(i) = C_new(7) + C_new(8)*T + C_new(9)*T^2 + C_new(10)*T^3 + C_new(11)*T^4 + C_new(12)*T^6;
    end
end


figure(1);
plot(t_sample,s_front);
hold on;
plot(t_sample,s_behind);
hold on;
plot(t_sample,s_ego,'o');
hold on;
plot(t_sample,s_ego_new);



