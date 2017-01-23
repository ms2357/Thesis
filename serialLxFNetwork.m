function [abserr]=test(a,b,N,L)

%initialize parameters
t0 = 0;
tf = 1;

%define space mesh
dx = 1 / L;
x = a: dx : b;

% define time mesh
dt = (tf - t0) / N;
t = t0 : dt : tf;

%wave speed
a = .8;
%CFL number
mu = a * dt / dx;

%set inital value fucntion  
f = sin(pi * x);
%set boundary value function 
g = sin(-pi * a * t);

%preallocate u, set initial u value to f, boundary value to g
u = zeros(N+1,L+1);
uexact = zeros(N+1,L+1);
u(1,:) = f;
u(:,1) = g;


%loop through time
for n=1:N
    %loop through space
    for j=2:L
        u(n + 1 , j) = (u(n , j - 1) + u(n , j + 1)) / 2 ...
            - (mu / 2) * (u(n , j + 1) - u(n , j - 1));
     end
    
    %ghost cell values, outflow boundary 
    u(n+1 , j + 1) = 2 * u(n , L) - u(n , L - 1);
   
    X = x - a * t(n + 1);
    uexact(n + 1 , :) = sin(pi * X);
    
    %calculate the abs error and infinity norm of the error
    abserr = abs(uexact-u); 
    
end
u(N , L + 1) = 2 * u(n , L) - u(n , L - 1);

