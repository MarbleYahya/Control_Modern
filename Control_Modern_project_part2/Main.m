clear;
close;
clc;

s = tf('s');
syms z
%parameters
mb = 0.11; %mass of ball(kg)
Jb = 1.76*10^-5; %rotational moment of inertia of the ball(kgm^2)
Rb = 0.02; %radius of the ball(m)
Jpx = 0.5; %rotational moment of inertia of the plate(kgm^2)
g = 9.8; %acceleration of the gravity


%states
syms x1 y1 x1_dot y1_dot %position and velocity of ball in x-axis and y-axis
syms a1 b1 a1_dot b1_dot %plate angular deflection and velocity in x-axis and y-axis

syms a1_dot2 b1_dot2
x = [x1,x1_dot,a1,a1_dot,y1,y1_dot,b1,b1_dot];
u = [a1_dot2 , b1_dot2];
B_parameter = mb/(mb + (Jb/Rb^2));

%state matrices(non linear)
A1 = [x(2),B_parameter*(x(1)*x(4)^2 + x(4)*x(5)*x(8) -g*sin(x(3))),x(4),0];
A2 = [x(6),B_parameter*(x(5)*x(8)^2 + x(1)*x(4)*x(8) -g*sin(x(7))),x(8),0];
A = [A1,A2]';
B = [0,0;0,0;0,0;1,0;0,0;0,0;0,0;0,1];
C = [1,0,0,0,0,0,0,0;0,0,0,0,1,0,0,0];

%linearization
x_balance = [0,0,0,0,0,0,0,0];
A_linear = jacobian(A,x);
A_numerical = double(subs(A_linear,[x1,x1_dot,a1,a1_dot,y1,y1_dot,b1,b1_dot],[0,0,0,0,0,0,0,0]));
B_numerical = B;
C_numerical = C;
D_numerical = [0,0;0,0];

%Q3:jordan matrix
[V,J] = jordan(A_numerical);
disp(J);

%Q4:transfer function
Tt = C_numerical*inv(s*eye(8) - A_numerical)*B_numerical;
%display transfer function
disp(Tt);

%Q4:pole-zero map
figure("Name", "Zero-Pole Plot for T(1,1)");
pzmap(Tt(1,1));
figure("Name", "Zero-Pole Plot for T(1,2)");
pzmap(Tt(1,2));
figure("Name", "Zero-Pole Plot for T(2,1)");
pzmap(Tt(2,1));
figure("Name", "Zero-Pole Plot for T(2,2)");
pzmap(Tt(2,2));
 
%Q5 : internal stability
Q = eye(size(A_numerical));
try
    P = lyap(A_numerical',-Q);  %system is not stable
    disp('system is stable according to lyapanov');
catch
    disp('system is not stable according to lyapanov');
end

eig_sys = eig(A_numerical);
if all(real(eig_sys) < 0)
    disp('system is stable according to eigen values');
else
    disp('system is not stable according to eigen values');
end

%Q6 : state transfer matrix

x0 = [2;0;0;1;0;0;0;1];
%time vector
t = linspace(0,5,100);
%inputs
u1 = ones(size(t));
u2 = ones(size(t));
u_signal = [u1;u2];
u_signal = reshape(u_signal , [] , 2);
%state space system
sys = ss(A_numerical,B_numerical,C_numerical,D_numerical);
%state response
[y_sys,t,x_sys] = lsim(sys ,u_signal , t , x0);

%plot state space response
figure;
plot(t,x_sys);
title('State Response');
xlabel('Time(s)');
ylabel('State values');
legend('x1','x2','x3','x4','x5','x6','x7','x8');
%plot system output response
figure;
plot(t,y_sys);
title('System Output Response');
xlabel('Time(s)');
ylabel('Outputs');
legend('y1','y2');

%Display state transfer matrix
syms t_syms
disp(exp(A_numerical*t_syms));

%display state response
disp(x_sys);

%display system output response
disp(y_sys);

%Q7 : 
[V,J] = jordan(A_numerical);
W = ones(8,1); %vector without zeros --- inv(V) * M = W --> M =V*W
initial_values = V*W;
disp('initial values');
disp(initial_values);

%Q10 :
C_control = ctrb(A_numerical,B_numerical);
O_observe = obsv(A_numerical,C_numerical);
disp('ctrb');
disp(C_control);
disp('obsv');
disp(O_observe);

% Rank of Controllability and Observability matrices
rank_C = rank(C_control);
disp(rank_C);
rank_O = rank(O_observe);
disp(rank_O);
% Number of states in the system (dimension of A matrix)
n_states = size(A_numerical, 1);

% Check controllability
if rank_C == n_states
    disp('The system is fully controllable');
else
    disp('The system is NOT fully controllable');
    % Identifying uncontrollable poles
    uncontrollable_poles = eig(A_numerical(~all(C_control, 2), :));
    disp('Uncontrollable poles:');
    disp(uncontrollable_poles);
end

% Check observability
if rank_O == n_states
    disp('The system is fully observable');
else
    disp('The system is NOT fully observable');
    % Identifying unobservable poles
    unobservable_poles = eig(A_numerical(~all(O_observe, 2), :));
    disp('Unobservable poles:');
    disp(unobservable_poles);
end


% Check if the system is detectable (if unobservable modes are stable)
if rank_O < n_states
    detectable = all(real(unobservable_poles) < 0);
    if detectable
        disp('The system is detectable (unobservable modes are stable)');
    else
        disp('The system is NOT detectable');
    end
end

% Check if the system is stabilizable (if uncontrollable modes are stable)
if rank_C < n_states
    stabilizable = all(real(uncontrollable_poles) < 0);
    if stabilizable
        disp('The system is stabilizable (uncontrollable modes are stable)');
    else
        disp('The system is NOT stabilizable');
    end
end
%% kalman
A_N=A_numerical;
B_N=B_numerical;
C_N=C_numerical;
D_N=D_numerical;
cc=ctrb(A_N,B_N); %check controlpaziri matris asli
%rref(cc); %check satri moghadamati controlpaziri matris asli
%rank(cc); %check rank controlpaziri matris asli
T=cc(1:8,1:8);%tashkil matris T az bordarhai mostaghel ctrb(A_N,B_N)

A_bar=T\A_N*T;
B_bar=T\B_N;
C_bar=C*T;
A_bar_C=A_bar(1:4,1:4);
B_bar_C=B_bar(1:4,1:2);
C_bar_C=C_bar(1:2,1:4);
%rref(ctrb(A_bar_C,B_bar_C))
%rank(ctrb(A_bar_C,B_bar_C))
ob_c=obsv(A_bar_C,C_bar_C);
%rref(obsv(A_bar_C,C_bar_C))
%rank(obsv(A_bar_C,C_bar_C))
U_C=[1 0 0 0;
     0 1 0 0;
     0 0 1 0;
     0 0 0 1];
A_bar_CC=U_C*A_bar_C/U_C;
B_bar_CC=U_C*B_bar_C;
C_bar_CC=C_bar_C/(U_C);

%control napaziri
A_bar_C_bar=A_bar(5:8,5:8);
B_bar_C_bar=B_bar(5:8,1:2);
C_bar_C_bar=C_bar(1:2,5:8);
%ctr_c_bar=ctrb(A_bar_C_bar,B_bar_C_bar) %check barai control napaziri
%matris
%rref(ctr_c_bar) %check barai control napaziri
%matris
%rank(ctr_c_bar) %check barai control napaziri
%matris
ob_c_bar=obsv(A_bar_C_bar,C_bar_C_bar);
%rref(ob_c_bar)
%rank(ob_c_bar)
U_C_bar=[1 0 0 0;
     0 1 0 0;
     0 0 1 0;
     0 0 0 1];
A_bar_C_bar_bar=U_C_bar*A_bar_C_bar/U_C_bar;
B_bar_C_bar_bar=U_C_bar*B_bar_C_bar;
C_bar_C_bar_bar=C_bar_C_bar/(U_C_bar);




%Q2 ----

sys = ss(A_numerical,B_numerical,C_numerical,D_numerical);
irreducable_system = minreal(sys);

%Q3 --- 

wn = 20;% rad/secs ,, according to paper
Ts = 0.4;% setteling time(secs) ,, according to paper
zeta = 4/(Ts*wn);
PO = exp((-pi*zeta)/(sqrt(1-zeta^2))); %overshoot

%desired poles ------
s1 = -zeta*wn + i*wn*sqrt(1-zeta^2);  
s2 = -zeta*wn - i*wn*sqrt(1-zeta^2);

%------
%note : system is completely controllable otherwise the uncontrollable
%poles should remain
%note 2: reapeated poles should not be more than B matrix's rank

%other poles----
s_p = real(s1)*5;
poles = [s1,s_p,s_p-1,s_p-2,s2,s_p-3,s_p-4,s_p-5];

%find state feedback values
K_statefeedback= place(A_numerical,B_numerical,poles);


%% small test


A_closed = A_numerical - B_numerical * K_statefeedback;

sys_closed = ss(A_closed, B_numerical, C_numerical, D_numerical);

step_info = stepinfo(sys_closed);

disp('Step Response Information:');
disp(step_info);

%% Static Pre-compensator
A_cl = A_numerical-B_numerical*K_statefeedback;
G_CL_0 = -C_numerical/(A_cl)*B_numerical ; % -C*inv(A_cl)*B
inv_G_CL_0 = inv(G_CL_0);

%% Dynamic Pre-compensator
n = length(A_numerical);
    [~,p] = size(B_numerical);
    [l,~] = size(C_numerical);
     D = zeros(l,p);
%      
% syms k_1 k_2 k_3 k_4 k_5 k_6 k_7 k_8 k_9 k_10 k_11 k_12 k_13 k_14 k_15 k_16 k_q1 k_q2     
% 
%     k_q = [k_q1 ; k_q2];
% 
Ai = [A_numerical zeros(n,l); -C_numerical zeros(l,l)];
Bi = [B_numerical;  zeros(l,p)];
% %syms k_1 k_2 k_3 k_4 k_5 k_6 k_7 k_8 k_9 k_10 k_11 k_12 k_13 k_14 k_15 k_16 k_q1 k_q2
% K_d = [k_1 k_2 k_3 k_4 k_5 k_6 k_7 k_8 ; k_9 k_10 k_11 k_12 k_13 k_14 k_15 k_16];
% K_D_S = [k_1 k_2 k_3 k_4 k_5 k_6 k_7 k_8 k_q1 ; k_9 k_10 k_11 k_12 k_13 k_14 k_15 k_16 k_q2];
% [~,p] = size(B_numerical);
% if (rank([B_numerical A_numerical; zeros(l,p) -C_numerical]) == length(Ai))  % Controllability Condition for Dynamic P.Comp.
%     disp ('Augmented system with Integral State is controllable');
%     A_CL_D = [A_numerical-B_numerical*K_d B_numerical*k_q;
%               -C_numerical 0];%zeros(2, 2)
%          % system order 
% 
% 
%     CE_D = vpa(det(s*eye(n+1)-A_CL_D),3);
%     vpa(collect (CE_D,s),3)
%     coeff_D = coeffs(CE_D,s);
%     coeff_D = vpa(fliplr(coeff_D),3);
%     alpha_D = (s+s1)*(s+s2)*(s+s_p)*(s+s_p-1)*(s+s_p-2)*(s+s_p-3)*(s+s_p-4)*(s+s_p-5); % Desired Poles
%     coeff_alpha = coeffs(alpha_D,s);
%     coeff_alpha = vpa(fliplr(coeff_alpha),3);
%     for i = 1:length(A_CL_D)
%         equns(i)=coeff_D(i+1)-coeff_alpha(i+1);
%     end
%     K_D = solve(equns,[k_1 k_2 k_3 k_4 k_5 k_6 k_7 k_8 k_q1 ; k_9 k_10 k_11 k_12 k_13 k_14 k_15 k_16 k_q2]);
%     K_d = double([K_D.k_1 K_D.k_2 K_D.k_3 K_D.k_4 K_D.k_5 K_D.k_6 K_D.k_7 K_D.k_8 K_D.k_9 K_D.k_10 K_D.k_11 K_D.k_12 K_D.k_13 K_D.k_14 K_D.k_15 K_D.k_16]);
%     k_q1 = double(K_D.k_q1);
%     k_q2 = double(K_D.k_q2);
%     k_q = [k_q1 ; k_q2];
% else
%     disp ('Augmented system with Integral State is not controllable')
% end

% % Design parameters for LQR
% Q_state = eye(n);          % State weighting matrix
% Q_integral = eye(q) * 100; % Integral state weighting matrix
% Q_aug = blkdiag(Q_state, Q_integral);
% R = eye(p);                % Control weighting matrix

Q = diag([100,100,100,100,100,100,100,100,10,10]);
R = diag([1,1]);
% [K,S,P] = lqr(A,B,Q,R,N) calculates the optimal gain matrix K, the solution S of the associated
% algebraic Riccati equation and the closed-loop poles P using the continuous-time state-space
% matrices A and B.
% S — Solution of the associated algebraic Riccati equation
% P — Poles of the closed-loop system
[K_LQR,~,P_LQR] = lqr(Ai,Bi,Q,R,0);
K_d_LQR = K_LQR(:,1:n);
K_q_LQR = K_LQR(:,n+1:n+p); % use -K_q_LQR in Simulink 
