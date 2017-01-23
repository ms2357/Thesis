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
rho0 = 2;                 %density
K0 = .25;                   %bulk modulus
rho1=3;
K1=.33;
%impedence, does this vary????
Z0 = rho0 * K0;             
Z1=rho1*K1;
%Zmatrix=[Z1 Z0;1 1];

A = [0 K0; 1 / rho0 0];       %wave speed
R = [-Z0 Z0 ; 1 1];         

%CFL number
mu = A*dt/dx;

%preallocate u, set initial value using intital condition
U =struct('PressureEdge1', zeros(N+1,L+1), 'VelocityEdge1',zeros(N+1,L+1),...
    'PressureEdge2', zeros(N+1,L+1), 'VelocityEdge2',zeros(N+1,L+1));
U(:).PressureEdge1(1,:) = InitialPressure(1:N+1);
U(:).VelocityEdge1(1,:) = InitialVelocity(1:N+1);
U(:).PressureEdge2(1,:) = InitialPressure(N+1:end);
U(:).VelocityEdge2(1,:) = InitialVelocity(N+1:end);



%struct to store w,z at junctions
%w left moving 
%z right moving
C=struct('w0', zeros(N+1,2), 'z0',zeros(N+1,2),...
    'w1', zeros(N+1,2), 'z1',zeros(N+1,2));

%loop through time
for n=1:N
    %loop through space
    for j=2:L
%u(n + 1,j)= (u(n,j - 1) + u(n,j + 1))/2 -(mu/2)*(u(n,j + 1) - u(n,j - 1))

      prev=(cell2mat({U.PressureEdge1(n,j-1),U.VelocityEdge1(n,j-1)}))';
      next=(cell2mat({U.PressureEdge1(n,j+1),U.VelocityEdge1(n,j+1)}))';
      newU= (prev + next)/ 2 - (mu/ 2) * (next - prev);

      U.PressureEdge1(n + 1 , j)=newU(1);
      U.VelocityEdge1(n + 1 , j)=newU(2);
      
      prev=(cell2mat({U.PressureEdge2(n,j-1),U.VelocityEdge2(n,j-1)}))';
      next=(cell2mat({U.PressureEdge2(n,j+1),U.VelocityEdge2(n,j+1)}))';
      newU= (prev + next)/ 2 - (mu/ 2) * (next - prev);

      U.PressureEdge2(n + 1 , j)=newU(1);
      U.VelocityEdge2(n + 1 , j)=newU(2);
    end

    %Jucntions
    j=1;

    %calc chars using pres/vel
    prev=(cell2mat({U.PressureEdge1(n,j),U.VelocityEdge1(n,j)}))';
    prev= R*prev;
    next=(cell2mat({U.PressureEdge1(n,j+1),U.VelocityEdge1(n,j+1)}))';
    next= R*next;
    
    %z(x^n)=mu*(z_(L)-z_(L-1))+z_L...linear interpolation for known char
    curr = prev - mu*(prev - next);
    %set w char value at junction 
    C.w0(n+1,j)=curr(1);

    %calc chars using pres/vel
    j = L+1;
    prev=(cell2mat({U.PressureEdge1(n,j-1),U.VelocityEdge1(n,j-1)}))';
    prev=R*prev;
    next=(cell2mat({U.PressureEdge1(n,j),U.VelocityEdge1(n,j)}))';
    next=R*next;
    
    %w(x^n)=mu*(z_(L)-z_(L-1))+z_L...linear interpolation for known char
    curr = mu*(prev - next) + prev;
    
    %set w char value at junction 
    C.z1(n+1,2)=curr(2);
    
    %calc unk chars at juction using known chars.  Solve system Ax=b using
    %interpolated chars 
    knownchars = (cell2mat({C.w0(n+1,1),C.z1(n+1,2)}))';
    knownchars = [-Z0 -Z1;-1 1]*knownchars;
    unkchars = [-Z1 -Z0; -1 1] \ knownchars ;
    C.w1(n+1,2) = unkchars(1);
    C.z0(n+1,1) = unkchars(2);
  
    %calc press/vel using chars
    j=1;
    newC = (cell2mat({C.w0(n+1,1),C.z0(n+1,1)}))';
    newU = R \ newC;
    U.PressureEdge1(n+1,j) = newU(1);
    U.VelocityEdge1(n+1,j) = newU(2);
    
    %calc press/vel using chars
    j=L+1;
    newC = (cell2mat({C.w1(n+1,2),C.z1(n+1,2)}))';
    newU = R \ newC;
    U.PressureEdge1(n+1,j) = newU(1);
    U.VelocityEdge1(n+1,j) = newU(2);

end

clf
%figure(glf)
% hold on
for i=1:L
   
    plot(x(1:N+1),U.PressureEdge1(i,:),x(1:N+1),U.VelocityEdge1(i,:))
	axis([0 1 -1.5 1.5])
    pause(.05)
    drawnow
    
end
hold off
