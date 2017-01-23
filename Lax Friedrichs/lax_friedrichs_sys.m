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
dt = dx;
t = t0:dt:tf;

%wave speed details
rho0 = 1;                 %density
K0 = 1;                   %bulk modulus
Z0 = rho0 * K0;             %impedence
A = [0 K0; 1 / rho0 0];       %wave speed
R = [-Z0 Z0 ; 1 1];         

%impose boundary functions for chars
bd0=R*[sin(pi*t);cos(pi*t)];
bdL=R*[-sin(pi*(L-t));-cos(pi*(L-t))];
f=bd0(1,:)';
g=bdL(2,:)';
 
%CFL number
mu = A*dt/dx;

%preallocate u, set initial u value to f
U =struct('pressure', zeros(N+1,L+1), 'velocity',zeros(N+1,L+1));
U(:).pressure(1,:) = InitialPressure;
U(:).velocity(1,:) = InitialVelocity;


%struct to store imposed and interpolated values at boundaries
%col 1 stores values at left bdd col 2 stores values at right bdd
%w left moving 
%z right moving
C=struct('w', zeros(N+1,2), 'z',zeros(N+1,2));

%set col 2 of z char to imposed bdd function values 
C.w(:,2)=g;
%set col 1 of z char to imposed bdd function values 
C.z(:,1)=f;



%loop through time
for n=1:N
    %loop through space
    for j=2:L
    %u(n + 1,j)=(u(n,j - 1) + u(n,j + 1))/2
    %           -(mu/2)*(u(n,j + 1) - u(n,j - 1))
      
      prev=(cell2mat({U.pressure(n,j - 1),U.velocity(n,j - 1)}))';
      next=(cell2mat({U.pressure(n,j + 1),U.velocity(n,j + 1)}))';
      newU= (prev + next)/ 2 - (mu/ 2) * (next - prev);

      U.pressure(n + 1 , j)=newU(1);
      U.velocity(n + 1 , j)=newU(2);
    end
    
    %Boundary conditions
    j=1;
    
    prev=R*(cell2mat({U.pressure(n,j),U.velocity(n,j)}))';
    next=R*(cell2mat({U.pressure(n,j+1),U.velocity(n,j+1)}))';
    
    %z(x^n)=mu*(z_(L)-z_(L-1))+z_L...linear interpolation for unk char
    curr = mu*(prev-next) + prev ;
    C.w(n,j)=curr(1);
    
    %calc press/vel using chars
    newC=(cell2mat({C.w(n,j),C.z(n,j)}))';
    newU=R\newC;
    U.pressure(n+1,j)=newU(1);
    U.velocity(n+1,j)=newU(2);
    
       
    j = L+1;
    prev=R*(cell2mat({U.pressure(n,j-1),U.velocity(n,j-1)}))';
    next=R*(cell2mat({U.pressure(n,j),U.velocity(n,j)}))';
    
    %w(x^n)=mu*(z_(L)-z_(L-1))+z_L...linear interpolation for unk char
    curr = mu*(prev-next) +prev ;
    C.z(n,2)=curr(2);
    
    %calc press/vel using chars
    newC=(cell2mat({C.w(n,2),C.z(n,2)}))';
    newU=R\newC;
    U.pressure(n+1,j)=newU(1);
    U.velocity(n+1,j)=newU(2);
    
end

clf
for i=1:L
   
    plot(x,U.pressure(i,:),x,U.velocity(i,:))
	axis([0 1 -1.5 1.5])
    pause(.05)
    drawnow
end
hold off

