function y = solvePKE(lambda, beta, beta_sum, L, target, h, ...
f_case, init_cond, v, time_insert, rho_ex)

% Determine number of delayed groups
m = length(lambda) + 1;

x = init_cond; time = 0;
f_hat = zeros(m,1);
d_hat= zeros(m,m);
big_d = zeros(m,m);
i = 1;
iterations = target / h;
result = zeros(m+1,iterations);

% Control rod positioning
c = [0.24079920041519   1.55078931457927];
rho = @(h) c(1) - c(2)*(1-h/36+sin(2*pi*h/36));
height = fzero(@(h) rho(h) - rho_ex, 20);

% Begin time dependent iterations
d = zeros(m+1,1);
while time <= target
    time = time + h;

    % Calculate the values of the reactivity and source at the midpoint
    mid_time = (time + time+h)/2;
%     if time < 720
%         p = rho(3*mid_time/60)*beta_sum;
%     elseif time >= 720 && time < time_insert
%         p = rho_ex*beta_sum;
    if time < time_insert
        p = rho_ex*beta_sum;
    else
        p = rho(height - v*(mid_time - time_insert)/60)*beta_sum;
    end
    source = f(f_case, mid_time);

    % Caculate the roots to the inhour equation
    d = inhour(lambda, L, beta, p, d);

    % Calculate the eigenvectors and the inverse of the matrix of eigenvectors
    Y = ev2(lambda, L, beta, d);
    Y_inv = ev_inv(lambda, L, beta, d);

    % Construct matrices for computation
    for k = 1:m
        d_hat(k,k) = exp(d(k)*h);
        big_d(k,k) = d(k);
    end
    f_hat(1) = source;
    big_d_inv = zeros(m,m);
    for k =1:m
        big_d_inv(k,k) = 1/big_d(k,k);
    end

    % Compute next time step
    x = (Y * d_hat * Y_inv)*(x + (Y*big_d_inv*Y_inv*f_hat)) - ...
    (Y*big_d_inv*Y_inv*f_hat);

    % Store results in a matrix
    result(1,i) = time;
    for j = 1:m
        result(j+1,i) = x(j);
    end
    %%Update counters
    i = i + 1;
end
y=result;

function y = ev2(lambda,L, beta, evals)
%This is a simple function that calculates the eigenvectors using the
%appropriate forms.
m = length(lambda) + 1;
evects = zeros(m,m);
for i = 1:m
    for j = 1:m
        if i == 1
            evects(i,j) = 1;
        end
        if i~= 1
            mu = beta(i-1)/L;
            evects(i,j) = mu / (lambda(i-1) + evals(j));
        end
    end
end
y = evects;

function y = ev_inv(lambda, L, beta, evals)
% This function returns the inverse of the matrix of eigenvalues
% based on some computations provieded in Aboanber and Nahla.
m = length(lambda) + 1;
mu = zeros(1,m-1);
for i = 1:m-1
    mu(i) = beta(i)/L;
end
normfact = zeros(m,1);
for k = 1:m
    sum = 0;
    for i = 1:m-1
        temp= mu(i)*lambda(i);
        temp2 = (lambda(i) + evals(k))^2;
        temp3 = temp/temp2;
        sum = sum+temp3;
    end
    normfact(k) = 1 / (sum+1);
end
result = zeros(m,m);
for i = 1:m
    for j = 1:m
        if i == 1
            result(i,j) = 1*normfact(j);
        end
        if i~=1
            result(i,j) = (lambda(i-1) / (lambda(i-1) + evals(j)))*normfact(j);
        end
    end
end
y=transpose(result);

function y = expand(lambda)
% A simple helper function to provide the coefficients of a polynomial
% produced by raising the function (x+y) to the nth power.
% The argument, lambda is a vector of constants that are needed to
% derive the coefficients.
% Determines the number of iterations, as well as the degree
% of the polynomial in question
m = length(lambda);
coeff = zeros(m+1,1);
% A temporary variable is necessary b/c the iterations that follow
% require information from the previous iteration...
temp = coeff;
% Must run the index to m+1 b/c MATLAB uses a 1-based index
for i= 1:m+1
    if i ~= 1
        coeff(1) = temp(1) * lambda(i-1);
        for j = 2:m+1
            coeff(j) = temp(j)*lambda(i-1) + temp(j-1);
            if j == i-1
                coeff(j) = temp(j-1) + lambda(i-1);
            end
        end
    end
    coeff(i) = 1;
    temp = coeff;
end
y = coeff;

function y = f(case_number, time)
% This is the source function for our solution. It works just like rho.m
% in that the case_number determines the function to use.
% case_number=1 : f = 0
if case_number == 1
    result = 0;
end
if case_number == 2
    result = 7.5e-2;
end
y=result;

function y = inhour(lambda, L, beta, rho , init_root)
% This function begins by taking the arguments and converting them into
% the correct m-degree polynomial inorder to take advantage of the given
% method of finding the roots of said polynomial.
m = length(lambda);
sum = zeros(m,1);
coeff = expand(lambda);
coeff_2 = zeros(m+2,1);
for i = 2:m+1
    coeff_2(i) = rho*coeff(i) - L*coeff(i-1);
end
coeff_2(1) = rho*coeff(1);
coeff_2(m+2) = -L * coeff(m+1);
for i = 1:m
    temp_lambda = trunc(lambda, i);
    temp = beta(i)* expand(temp_lambda);
    sum = temp + sum;
end
sum = -1*sum;
res = zeros(m+2,1);
for i =1:m
    res(i+1) = coeff_2(i+1) + sum(i);
end
res(1) = coeff_2(1);
res(m+2) = coeff_2(m+2);
e_vals = rootfinder(res, init_root, .00001);
y = e_vals;

function y = myDeriv(coeff)
% A simple function that calculates the derivitive
% coefficient vector for a given polynomial.
deg = length(coeff);
if deg ~= 1
    result = zeros(1,deg - 1);
    for i= 1: (deg-1)
        result(i) = coeff(i+1) * i;
    end
end
if deg == 1
    result = 0;
end
y = result;

function y = myEval(coeff, x)
% Evaluates the polynomial expressed as coeff at the value x.
deg = length(coeff);
sum = coeff(1);
if deg ~= 1
    for i = 2:deg
        sum = sum + coeff(i) * x^(i-1);
    end
end
y = sum;

function y=myHorner(a,z,n)
%Applies a functional implementation of the Horner mehtod
%The user supplies a(The poly), z(The root), and n(The degree)
%This program uses Horner's method to write p(x) = (x-z)q(x)+c
%Where p and q are polynomials of degree n and n-1 respectively
b = zeros(1,n+1);
b(n) = a(n+1);
if n>0
    for i = 1:n
        b(n+1-i) = a(n-i+2) + b(n+2-i)*z;
    end
    c= a(1) + b(1)*z;
end
a = b;
ret = a;
ret(n) = ret(n) + c; %add the constant
y=ret;

function y = newton(val, poly, tol)
% A simple implementation of Newton's Method
eps = 1;
x = val;
deriv = myDeriv(poly);
while eps > tol
    temp = x - (myEval(poly,x) / myEval(deriv, x));
    eps = abs(x - temp);
    x = temp;
end
y = x;

% function y = rho(beta_sum, t, t_insert, v, rho_ex)
% % This function represents the time-dependent reactivity function
% % for the point kinetics equation.
% c = [-0.48158745700046   0.90405552235442 ...
%       0.07713328005310  -1.00502376863748];
% R = @(t) c(1) + c(2)*sin(c(3)*(height-v*t/60) + c(4));
% result = R(t-t_insert);
% result = result*beta_sum;
% y=result;

function y = rootfinder(coeff, init, tol)
% This is a simple wrapper function that takes an coefficent vector
% and uses Newton's method to find all of the real roots of said poly.
% The function takes advantage of Horner's method to deflate the poly
% at each step to expedite computation. The argument init is a vector
% of initial values that are used in Newton's method.
deg = length(coeff) - 1;
result = zeros(deg,1);
counter=1;
while deg > 1
    result(counter) = newton(init(counter),coeff,tol);
    coeff = myHorner(coeff, result(counter), deg) ;
    deg = deg - 1;
    counter = counter + 1;
end
result(counter) = -coeff(1)/coeff(2);
y = result;

% function y = swap(arg)
% % This simple function takes the 1st element of arg and puts it in the
% % mth place of the resultant vector, and puts the 2nd in the m-1st...and
% % so on...
% m = length(arg);
% res = zeros(1,m);
% for i = 1:m
%     res(m-i+1) =arg(i);
% end
% y=res;

function y = trunc(var, t)
% This is a simple helper method that removes the ith element
% from the vector var.
m = length(var); flag = 0;
temp = zeros(1,m-1);
for i = 1:m
    if i ~= t && flag == 0
        temp(i) = var(i);
    end
    if i ~= t && flag == 1
        temp(i-1) = var(i);
    end
    if i == t
        flag = 1;
    end
end
y = transpose(temp);