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

%set inital funcitons  
InitialPressure = sin(pi * x)';
InitialVelocity = cos(pi * x)';

% define time mesh
dt = (tf-t0)/N;
t = t0:dt:tf;

%wave speed
A=[0 1; .5 0];

%CFL number
mu = A*dt/dx;

%preallocate u, set initial u value to f
U =struct('pressure', zeros(N+1,L+1), 'velocity',zeros(N+1,L+1));
U(:).pressure(1,:)=InitialPressure;
U(:).velocity(1,:)=InitialVelocity;


%loop through time
for n=1:N
    %loop through space
    for j=2:L
%u(n + 1,j)= (u(n,j - 1) + u(n,j + 1))/2 -(mu/2)*(u(n,j + 1) - u(n,j - 1))
      
      prev=(cell2mat({U.pressure(n,j-1),U.velocity(n,j-1)}))';
      %disp(first)  
      next=(cell2mat({U.pressure(n,j+1),U.velocity(n,j+1)}))';
      %disp(second)
      term1=(prev+next)/2;
      term2=-(mu/2)*(next-prev);
      new=term1+term2;

      U.pressure(n + 1 , j)=new(1);
      U.velocity(n+1,j)=new(2);
    end

    
    %Set values at boundaries using periodic BC's
    j=1; 
    
    % u(0) == u(N)
    prev=cell2mat({U.pressure(n , L ),U.velocity(n , L )})';
    next=cell2mat({U.pressure(n , j+1),U.velocity(n , j+1)})';
    new= (prev + next)/ 2 - (mu / 2) * (next - prev);
    U.pressure(n + 1 , j)=new(1);
    U.velocity(n + 1,j) = new(2); 
    %U.pressure(n + 1 , j)=U.pressure(n,L+1);
    %U.velocity(n + 1,j) = U.velocity(n,L+1);
    
        
      
    j = L;
    
    % u(N+1) == u(2)
    prev=cell2mat({U.pressure(n,j),U.velocity(n,j)})';
    next=cell2mat({U.pressure(n,2),U.velocity(n,2)})';
    new= (prev + next)./ 2 - (mu/ 2) * (next - prev);
    U.pressure(n + 1 , j + 1)=new(1);
    U.velocity(n + 1 , j + 1) = new(2); 
    %U.pressure(n + 1 , j+1)=U.pressure(n,1);
    %U.velocity(n + 1,j+1) = U.velocity(n,1);
    
    
end

clf
%figure(glf)
% hold on
for i=1:L
   
    plot(x,U.pressure(i,:),x,U.velocity(i,:))
	axis([0 1 -1 2])
    pause(.05)
    drawnow
end
hold off