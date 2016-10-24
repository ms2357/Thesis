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

%wave speed details
rho0 = 1;                 %density
K0 = 1;                   %bulk modulus
Z0 = rho0 * K0;             %impedence
A = [0 K0; 1 / rho0 0];       %wave speed
R = [-Z0 Z0 ; 1 1];         

%impose boundary functions w(0,t)= first value of R*u(0,t)=f(t) 
% and z(1,t)= second value of R*p(1,t)=g(t)
%essentially p+u=f(t)@0 p+u=g(t) @1
bdd=R*[sin(pi*t);cos(pi+pi*t)];
f=bdd(1,:)';
g=bdd(2,:)';

%CFL number
mu = A*dt/dx;

%preallocate u, set initial u value to f
U =struct('pressure', zeros(N+1,L+1), 'velocity',zeros(N+1,L+1));
U(:).pressure(1,:) = InitialPressure;
U(:).velocity(1,:) = InitialVelocity;



C=struct('w', zeros(N+1,2), 'z',zeros(N+1,2));
C.w(:,1)=f;
C.z(:,2)=g;



%loop through time
for n=1:N
    %loop through space
    for j=2:L
%u(n + 1,j)= (u(n,j - 1) + u(n,j + 1))/2 -(mu/2)*(u(n,j + 1) - u(n,j - 1))
      
      prev=(cell2mat({U.pressure(n,j-1),U.velocity(n,j-1)}))';
      next=(cell2mat({U.pressure(n,j+1),U.velocity(n,j+1)}))';
      newU= (prev + next)/ 2 - (mu/ 2) * (next - prev);

      U.pressure(n + 1 , j)=newU(1);
      U.velocity(n+1,j)=newU(2);
    end
    
    %Set values at boundaries using periodic BC's
    j=1;
    
    %need to do linear interp for z and then back solve for w then get u
    %and p from those here and on other boundary
    prev=R*(cell2mat({U.pressure(n,j),U.velocity(n,j)}))';
    next=R*(cell2mat({U.pressure(n,j+1),U.velocity(n,j+1)}))';
    curr = 2 * prev - next;
    C.z(n,j)=curr(2);
    newC=(cell2mat({C.w(n,j),C.z(n,j)}))';
    newU=R\newC;
    U.pressure(n+1,j)=newU(1);
    U.velocity(n+1,j)=newU(2);
    
       
    j = L;
    % u(L+1) == u(2)
    prev=R*(cell2mat({U.pressure(n,j),U.velocity(n,j)}))';
    next=R*(cell2mat({U.pressure(n,j-1),U.velocity(n,j-1)}))';
    curr = 2 * prev - next;
    C.w(n,2)=curr(1);
    newC=(cell2mat({C.w(n,2),C.z(n,2)}))';
    newU=R\newC;
    U.pressure(n+1,j)=newU(1);
    U.velocity(n+1,j)=newU(2);
    
end

clf
%figure(glf)
% hold on
for i=1:L
   
    plot(x,U.pressure(i,:),x,U.velocity(i,:))
	axis([0 1 -1 1])
    pause(.05)
    drawnow
end
hold off

