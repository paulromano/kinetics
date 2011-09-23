% Senior Design I: MANE-4380
% Code by Paul Kollath-Romano

% Response to Control Rod Insertion
% Point Kinetics Equation - Piecewise Continuous Approximation

function insertion(v)

% Reactor Constants
beta_eff = 0.00765;
lambda = [3.0100 1.1400 0.3010 0.1110 0.0305 0.0124];
beta = 1/beta_eff*[0.041 0.115 0.396 0.196 0.219 0.033];
L = 1e-2;
% beta_eff = 0.007;
% lambda = [0.0127 0.0317 0.155 0.311 1.4 3.87];
% beta = [0.000266 0.001491 0.001316 0.002849 0.00896 0.000182];
% L = 2e-5;

% Determine reactivity function
c = [-0.48158745700046   0.90405552235442 ...
      0.07713328005310  -1.00502376863748];
R = @(t) c(1) + c(2)*sin(c(3)*(36-v*t/60) + c(4));
h = 0.01; %Timestep

% Source-free equilibrium initial condition
iter = 100;
x = zeros(7,iter);
x(1,1) = 1;
for i = 1:6
    x(i+1,1) = beta(i)/(lambda(i)*L);
end

% Roots of Inhour Equation
mat = zeros(7);
mat(1,1) = (0.003-beta_eff)/L;
for i = 1:6;
    mat(1,i+1) = lambda(i);
    mat(i+1,1) = beta(i)/L;
    mat(i+1,i+1) = -lambda(i);
end
ev = eig(mat);

% Eigenvalue matrix and Eigenvectors
D = eye(7);
U = ones(7);
V = ones(7);
for i = 1:7
    D(i,i) = exp(ev(i));
    nu = 1;
    for j = 1:6
        U(j+1,i) = (beta(j)/L)/(lambda(j)+ev(i));
        V(i,j+1) = lambda(j)/(lambda(j)+ev(i));
        nu = nu + (beta(j)*lambda(j)/L)/(lambda(j) + ev(i))^2;
    end
    V(i,:) = V(i,:)/nu;
end