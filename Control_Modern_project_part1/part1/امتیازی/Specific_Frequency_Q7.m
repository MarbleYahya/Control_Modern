%Q7
% Define the system matrices A, B, C, D
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


%linearization
x_balance = [0,0,0,0,0,0,0,0];
A_linear = jacobian(A,x);
A_numerical = double(subs(A_linear,[x1,x1_dot,a1,a1_dot,y1,y1_dot,b1,b1_dot],[0,0,0,0,0,0,0,0]));
A = A_numerical;  % Linearized state matrix
B = [0, 0; 0, 0; 0, 0; 1, 0; 0, 0; 0, 0; 0, 0; 0, 1];  % Given B matrix
C = [1, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 1, 0, 0, 0];  % Given C matrix
D = zeros(size(C, 1), size(B, 2));  % D matrix is zero

%The frequency we want to suppress (in radians per second)
omega_0 = 2 * pi * 1;  % Example: Suppress frequency at 1 Hz (omega = 2*pi*1)

% Compute the eigenvalues and eigenvectors of A
[eigenvectors, eigenvalues_matrix] = eig(A); 
eigenvalues = diag(eigenvalues_matrix);  % Extract eigenvalues as a vector


disp('Eigenvalues of A:');
disp(eigenvalues);

% Identify the eigenvalue closest to the specified frequency (omega_0)
% We want to find the eigenvalue whose imaginary part is closest to omega_0
[~, idx] = min(abs(imag(eigenvalues) - omega_0));

% Get the corresponding eigenvector for the selected frequency
v = eigenvectors(:, idx);

% The initial condition should be orthogonal to this eigenvector to avoid
% exciting this mode. We compute an initial condition that is orthogonal.
x0 = null(v');  % Finding a vector orthogonal to the chosen eigenvector

% If the null space gives more than one vector, we choose the first one
if size(x0, 2) > 1
    x0 = x0(:, 1);  % Take the first vector
end

% Normalize the initial condition for better numerical stability
x0 = x0 / norm(x0);

disp('Initial condition to suppress the specific frequency:');
disp(x0);
