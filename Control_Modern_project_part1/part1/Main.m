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
T = C_numerical*inv(s*eye(8) - A_numerical)*B_numerical;
%display transfer function
disp(T);

%Q4:pole-zero map
figure("Name", "Zero-Pole Plot for T(1,1)");
pzmap(T(1,1));
figure("Name", "Zero-Pole Plot for T(1,2)");
pzmap(T(1,2));
figure("Name", "Zero-Pole Plot for T(2,1)");
pzmap(T(2,1));
figure("Name", "Zero-Pole Plot for T(2,2)");
pzmap(T(2,2));
 
%Q5 : internal stability
Q = eye(size(A_numerical));
try
    P = lyap(A_numerical',-Q);  %system is not stable
    disp('system is stable according to lyapanov');
catch
    disp('system is not stable according to lyapanov');
end

eig_sys = eig(A_numerical);
if all(eig_sys > 0)
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






