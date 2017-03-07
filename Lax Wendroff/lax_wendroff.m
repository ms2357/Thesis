
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Description: This script implements the Lax Friedrich method to 
% numerically solve the scalar advection equation u_t + au_x = 0. 
% It stores all the parameters that define the desired 
% method.  Stable provided mu=abs(ak/h)<=1 
% Lax Friedrichs scheme
% u_new(j) = (u_old(j-1)+u_old(j+1))/2 - (mu/2)*(u_old(j+1)-u_old(j-1))

% parameters:
%   t0:intitial time
%   tf:final time
%   L:number of spatial grid points 
%   N: number of time grid points
%   f: function that gives inital data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%initialize parameters
t0=0;
tf=1;
L=100;
N=100;


%define space mesh
dx = 1/L;
x = 0:dx:1;
x = x';

%set inital funciton  
f = exp(-(x-.5).^2);

% define time mesh
dt = (tf-t0)/N;
t = t0:dt:tf;

%wave speed
a = .8;
%CFL number
mu = a*dt/dx;

%preallocate u, set initial u value to f
u = zeros(N+1,L+1);
u(:,1) = f;



%loop through time
for j=1:L
    %loop through space
    for k=2:N
        u(k,j+1) = u(k,j) - (mu/2)*(u(k+1,j)-u(k-1,j)) +...
            (mu^2/2)*(u(k+1,j)-2*u(k,j)+u(k-1,j));
    end
    % BC's step function
    u(1,j+1)=1;
    u(N+1,j+1)=0;
     
 
end


clf
for i=1:L
   
    plot(x,u(:,i))
	axis([0 1 .7 1])
    pause(.05)
    drawnow
end
