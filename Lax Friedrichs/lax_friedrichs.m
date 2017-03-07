%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Description: This function implements the Lax Friedrich method to 
% numerically solve the scalar advection equation u_t + au_x = 0. 
% It stores all the parameters that define the desired 
% method.  Stable provided mu=abs(ak/h)<=1 
% Lax Friedrichs scheme
% u_new(j) = (u_old(j-1)+u_old(j+1))/2 - (mu/2)*(u_old(j+1)-u_old(j-1))
%
%
% parameters:
%   t0:intitial time
%   tf:final time
%   L:number of spatial grid points 
%   N: number of time grid points
%   f: function that gives inital data
%   g:function that computes ghost cell values 
%   a:wave speed, a>0 for right moving wave
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%initialize parameters
t0 = 0;
tf = 1;
L = 100;
N = 100;

%define space mesh
dx = 1 / L;
x = 0: dx : 1;

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
uexact(1,:) = u(1,:);

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


clf
for i=1:L
    plot(x , u( i , :) , 'b' , x , uexact(i , :) , 'r')
	axis([0 1 -1 1])
    pause(.05)
    drawnow
end

%calculate L1 error norm in space
row_err=dx*trapz(abserr,2);
%calculate L1 error norm in time 
err=dt*trapz(row_err);
disp(err)