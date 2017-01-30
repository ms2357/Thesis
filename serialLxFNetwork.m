a=1;
b=2;

%initialize parameters
t0=0;
tf=1;
L=100;
N=100;

%define space mesh
dx = 1/L;
x = 0:dx/2:1;
x = x';

%set inital funcitons  
InitialPressure = sin(pi * x)';
%InitialPressure = x;
%InitialVelocity = ones(1,length(x));
InitialVelocity = cos(pi*x)';


% define time mesh
dt = dx;
t = t0:dt:tf;

%wave speed details
rho1 = 2;                 %density
K1 = .25;                   %bulk modulus
rho2=2;
K2=.25;
%impedence
Z1 = rho1 * K1;             
Z2=rho2*K2;
%Zmatrix=[Z1 Z0;1 1];

A = [0 K1; 1 / rho1 0];       %wave speed
R = [-Z1 Z1 ; 1 1];         

%CFL number
mu = A*dt/dx;

%preallocate u, set initial value using intital condition
U =struct('PressureEdge1', zeros(N+1,L+1), 'VelocityEdge1',zeros(N+1,L+1),...
    'PressureEdge2', zeros(N+1,L+1), 'VelocityEdge2',zeros(N+1,L+1));
U(:).PressureEdge1(1,:) = InitialPressure(1:N+1);
U(:).VelocityEdge1(1,:) = InitialVelocity(1:N+1);
U(:).PressureEdge2(1,:) = InitialPressure(N+1:end);
U(:).VelocityEdge2(1,:) = InitialVelocity(N+1:end);



%struct to store w,z at junctions, index gives egde, column 1 of struct is
%the char at the left endpoint 'a' and column 2 gives the char at the 
% right endpoint 'b'
%w left moving 
%z right moving
C=struct('w1', zeros(N+1,2), 'z1',zeros(N+1,2),...
    'w2', zeros(N+1,2), 'z2',zeros(N+1,2));

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

    %calc chars using pres/vel for egde 1
    prev=(cell2mat({U.PressureEdge1(n,j),U.VelocityEdge1(n,j)}))';
    prev= R*prev;
    next=(cell2mat({U.PressureEdge1(n,j+1),U.VelocityEdge1(n,j+1)}))';
    next= R*next;
    %w(x^n)=mu*(z_(L)-z_(L-1))+z_L...linear interpolation for known char
    curr = prev - mu*(prev - next);
    %set w char value at junction 
    C.w1(n+1,a)=curr(1);
    
    %calc chars using pres/vel for edge 2
    prev=(cell2mat({U.PressureEdge2(n,j),U.VelocityEdge2(n,j)}))';
    prev= R*prev;
    next=(cell2mat({U.PressureEdge2(n,j+1),U.VelocityEdge2(n,j+1)}))';
    next= R*next;
    %w(x^n)=mu*(z_(L)-z_(L-1))+z_L...linear interpolation for known char
    curr = prev - mu*(prev - next);
    %set w char value at junction 
    C.w2(n+1,b)=curr(1);
    
    
    
    
    
    
    %right junction 
    j = L+1;
    %calc chars using pres/vel for edge 1
    prev=(cell2mat({U.PressureEdge1(n,j-1),U.VelocityEdge1(n,j-1)}))';
    prev=R*prev;
    next=(cell2mat({U.PressureEdge1(n,j),U.VelocityEdge1(n,j)}))';
    next=R*next;
    %z(x^n)=mu*(z_(L)-z_(L-1))+z_L...linear interpolation for known char
    curr = mu*(prev - next) + prev;
    %set w char value at junction 
    C.z1(n+1,b)=curr(2);
    
    %calc chars using pres/vel for edge 2
    prev=(cell2mat({U.PressureEdge2(n,j-1),U.VelocityEdge2(n,j-1)}))';
    prev=R*prev;
    next=(cell2mat({U.PressureEdge2(n,j),U.VelocityEdge2(n,j)}))';
    next=R*next;
    %z(x^n)=mu*(z_(L)-z_(L-1))+z_L...linear interpolation for known char
    curr = mu*(prev - next) + prev;
    %set w char value at junction 
    C.z2(n+1,a)=curr(2);
    
    
    
    
    %calc unk chars at juction using known chars.  Solve system Ax=b using
    %interpolated chars 
    knownchars = (cell2mat({C.w1(n+1,a),C.z2(n+1,a)}))';
    knownchars = [Z1 Z2;-1 1]*knownchars;
    unkchars = [Z2 Z1; -1 1] \ knownchars ;
    C.w2(n+1,a) = unkchars(1);
    C.z1(n+1,a) = unkchars(2);
  
    %calc press/vel using chars
    j=1;
    newC = (cell2mat({C.w1(n+1,a),C.z1(n+1,a)}))';
    newU = R \ newC;
    U.PressureEdge1(n+1,j) = newU(1);
    U.VelocityEdge1(n+1,j) = newU(2);
    %U.PressureEdge2(n+1,L+1)=U.PressureEdge1(n+1,j);
    
    newC = (cell2mat({C.w2(n+1,a),C.z2(n+1,a)}))';
    newU = R \ newC;
    U.PressureEdge2(n+1,L+1) = newU(1);
    U.VelocityEdge2(n+1,L+1) = newU(2);
    %U.PressureEdge1(n+1,L+1)=U.PressureEdge2(n+1,j);
    
    %calc press/vel using chars
    j=L+1;
    
    %calc unk chars at juction using known chars.  Solve system Ax=b using
    %interpolated chars 
    knownchars = (cell2mat({C.w2(n+1,b),C.z1(n+1,b)}))';
    knownchars = [-Z2 -Z1;1 -1]*knownchars;
    unkchars = [-Z1 -Z2; 1 -1] \ knownchars ;
    C.w1(n+1,b) = unkchars(1);
    C.z2(n+1,b) = unkchars(2);
    
    
    %set chars on edge 2 at L+1 to calc'd press/vel  
    newC = (cell2mat({C.w1(n+1,2),C.z1(n+1,2)}))';
    newU = R \ newC;
    U.PressureEdge1(n+1,j) = newU(1);
    U.VelocityEdge1(n+1,j) = newU(2);
    %U.PressureEdge2(n+1,1) = newU(1);
    
    newC = (cell2mat({C.w2(n+1,2),C.z2(n+1,2)}))';
    newU = R \ newC;
    U.PressureEdge2(n+1,1) = newU(1);
    U.VelocityEdge2(n+1,1) = newU(2);

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
