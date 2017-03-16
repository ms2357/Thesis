
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
f = sin(pi * x);

% define time mesh
dt = (tf-t0)/N;
t = t0:dt:tf;

%wave speed
a = .8;

%CFL number
mu = a*dt/dx;

%preallocate u, set initial u value to f
u = zeros(N+1,L+1);
u(1,:) = f;


%loop through time
for n=1:N
    %loop through space
    for j=2:L
        
        u(n + 1 , j) = (u(n , j - 1) + u(n , j + 1)) / 2 ...
            - (mu / 2) * (u(n , j + 1) - u(n , j - 1));
        
     end
    
    
    %Set values at boundaries using periodic BC's
    j=1; 
    
    % u(0) == u(N)
    %u(n+1,j)=(u(n,j+1)+u(n,L))/2-(mu/2)*(u(n,j+1)-u(n,L));
    u(n+1,j)=(u(n,j+1)+u(n,L-1))/2-(mu/2)*(u(n,j+1)-u(n,L-1));
    
    j = L+1;
    
    % u(N+2) == u(2)
    u(n+1,j)=(u(n,2)+u(n,j-1))/2-(mu/2)*(u(n,2)-u(n,j-1));
  
end


clf
%figure(glf)
% hold on
for i=1:L
   
    plot(x,u(i,:))
	axis([0 1 0 1])
    pause(.05)
    drawnow
end
hold off