clear all
clc

%% TUGAS BESAR SKM
% Nama  : Gede Haris Widiarta
% NIM   : 1102174038
% Title : Analysis on Cart Inverted Pendulum System
%%
% Define Parameter System
M = 1.5;                %Cart Mass (kg)
m = 0.3;                %Pendulum Mass (kg)
l = 0.4;                %Length of Rod (m)
g = 9.8;                %Gravity (m/s^2)
I = 0.099;              %Inertia Moment (kg/s^2)
b = 0.05;               %Cart Friction Coefficient (Ns/m)
bt = 0.005;             %Pendulum Damping Coefficient (Nms/rad)
p  = I*(M+m)+M*m*l^2;   %Denumerator
%%
% Define State Space Matrix
A = [       0            1     0          0       ;
      (m*g*l*(M+m))/p    0     0       (m*l*b)/p  ;
            0            0     0          1       ;
    -g*(m^2)*(l^2)/p     0     0  -(b*(I+m*l^2))/p;]

B = [0 ; -(m*l)/p ; 0 ; (I+m*l^2)/p];
C = [0 0 1 0 ; 1 0 0 0];
D = [0];
system = ss(A,B,C,D)
%%
% Open-Loop Response
disp('Poles Sistem: ')
poles = eig(A)
subplot(2, 2, 1:4)
hold on
pzmap(system)
grid on
title('Poles of Cart-Inverted Pendulum System')
%%
% Controllability
Matrix_Co = ctrb(A,B)
Rank_Co   = rank(ctrb(A,B))
if Rank_Co == Matrix_Co;
    fprintf('System is Controllable\n')
else
    fprintf('System is Uncontrollable\n')
end
%%
% Observeability
Matrix_Ob = obsv(A,C)
Rank_Ob   = rank(obsv(A,C))
if Rank_Ob == Matrix_Ob;
    fprintf('System is Observeable\n')
else
    fprintf('System is Unobserveable\n')
end
%%
% Pole Placement
pp  =[-11 -9 -1+2j -1-2j]
Kpp = place(A,B,pp)
Acp = A-B*Kpp;
sys_pp = ss(Acp,B,C,D);
% Plot
t = 0:0.01:10;
r =0.2*ones(size(t));
figure(2);clf
[y,t,x]=lsim(sys_pp,r,t);
[AX,H1,H2] = plotyy(t,y(:,1),t,y(:,2),'plot');
set(get(AX(1),'Ylabel'),'String','Cart Position (m)')
set(get(AX(2),'Ylabel'),'String','Pendulum Angle (rad)')
title('Pole Placement')
grid
%%
% Linear Quadratic Regulator (LQR)
% Define LQR Parameter
Q  = diag([40 40 400 0]);
R  = 1;
                
% Calculate LQR Gain
K  = lqr(A,B,Q,R);
Ac = A-B*K; %control matrix
sys_lqr1 = ss(Ac,B,C,D)
% Plot
t = 0:0.01:5;
r =0.2*ones(size(t));
figure(3);clf
[y,t,x]=lsim(sys_lqr1,r,t);
[AX,H1,H2] = plotyy(t,y(:,1),t,y(:,2),'plot');
set(get(AX(1),'Ylabel'),'String','Cart Position (m)')
set(get(AX(2),'Ylabel'),'String','Pendulum Angle (rad)')
title('LQR Response')
grid

% LQR dengan Feedforward Gain (Kr)
Cn = [0 0 1 0];                    %modifikasi matriks C, karena input ref. hanya diaplikasikan pada posisi cart
Kr = -inv(Cn*(inv(A-B*K))*B);      %menghitung Kr
sys_lqr2 = ss(Ac,B*Kr,C,D)
%Plot
t = 0:0.01:5;
r =0.2*ones(size(t));
figure(4);clf
[y,t,x]=lsim(sys_lqr2,r,t);
[AX,H1,H2] = plotyy(t,y(:,1),t,y(:,2),'plot');
set(get(AX(1),'Ylabel'),'String','Cart Position (m)')
set(get(AX(2),'Ylabel'),'String','Pendulum Angle (rad)')
title('LQR with Kr')
grid
%%
% State Estimator (Observer)
ob = obsv(sys_lqr1);                
observability = rank(ob)
ctr_poles     = eig(Ac)             
obsr_poles    = [-25 -26 -27 -28]   
L  = place(A',C',obsr_poles)'       

% Observer-Based State-Feedback Control
Aco = [(A-B*K) (B*K);zeros(size(A)) (A-L*C)];
Bco = [B;zeros(size(B))];
Cco = [C zeros(size(C))];
Dco = [0;0];
sys_ob = ss(Aco,Bco,Cco,Dco);
%Plot
t = 0:0.01:5;
r = 0.2*ones(size(t));
figure(5);clf
[y,t,x]=lsim(sys_ob,r,t);
[AX,H1,H2] = plotyy(t,y(:,1),t,y(:,2),'plot');
set(get(AX(1),'Ylabel'),'String','Cart Position (m)')
set(get(AX(2),'Ylabel'),'String','Pendulum Angle (rad)')
title('Observer-Based State-Feedback Control')
grid

% Observer-Based State-Feedback Control + Kr
Acl  = [(A-B*K) (B*K);zeros(size(A)) (A-L*C)];
Bcl  = [B;zeros(size(B))];
Ccl  = [C zeros(size(C))];
Dcl  = [0;0];
Ccln = [Cn zeros(size(Cn))];
Kr2  = -inv(Ccln*(Acl\Bcl))
Bclt = [B*Kr2;zeros(size(B))]; 
sys_ob2 = ss(Acl,Bclt,Ccl,Dcl);
% Plot
t = 0:0.01:5;
r = 0.2*ones(size(t));
figure(6);clf
[y,t,x]=lsim(sys_ob2,r,t);
[AX,H1,H2] = plotyy(t,y(:,1),t,y(:,2),'plot');
set(get(AX(1),'Ylabel'),'String','Cart Position (m)')
set(get(AX(2),'Ylabel'),'String','Pendulum Angle (rad)')
title('Observer-Based State-Feedback Control + Kr')
grid