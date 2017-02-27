%initialize parameters
t0 = 0;
tf = 1;
L=4;
N=4;

%define space mesh
a=0;
b=1;
dx = (b-a) / L;
x = a: dx : b;



% define time mesh
%dt = (tf - t0) / N;
dt=dx;
t = t0 : dt : tf;

%wave speed
c = .8;
%CFL number
mu = c * dt / dx;

%set inital value fucntion  
f = sin(pi * x);
%set boundary value function 
g = sin(-pi * c * t);

%preallocate u, set initial u value to f, boundary value to g
u = zeros(N+1,L+1);
uexact = zeros(N+1,L+1);
u(1,:) = f;
u(:,1) = g;
uexact(1,:)=u(1,:);
row_err=zeros(N+1,1);
%Error=zeros(N+1,1);
abserr=zeros(N+1,L+1);
   
for mm=1:N+1
    
        
        %loop through time
    for n=1:N
        %loop through space
        for j=2:L
            u(n + 1 , j) = (u(n , j - 1) + u(n , j + 1)) / 2 ...
                - (mu / 2) * (u(n , j + 1) - u(n , j - 1));
         end

        %ghost cell values, outflow boundary 
        u(n+1 , j + 1) = 2 * u(n , L) - u(n , L - 1);

        X = x - c * t(n + 1);
        uexact(n + 1 , :) = sin(pi * X);
    
    end
    
    %calculate the abs error and infinity norm of the error
    abserr = abs(uexact-u);
    %calculate L1 error norm in space
    row_err=dx*trapz(abserr,2);
    %calculate L1 error norm in time 
    col_err = dx*trapz(row_err);
    Error(mm)=col_err;
    
    %reassign parameters
    %space mesh
    dx=dx/2;
    x = a: dx : b;
    N=2*N;
    %time mesh
    dt=dx;
    t = t0 : dt : tf;
    L=2*L;
    %CFL number
    mu = c * dt / dx;
    
    %set inital value fucntion  
    f = sin(pi * x);
    %set boundary value function 
    g = sin(-pi * c * t);

    %preallocate u, set initial u value to f, boundary value to g
    u = zeros(N+1,L+1);
    uexact = zeros(N+1,L+1);
    u(1,:) = f;
    u(:,1) = g;
    uexact(1,:)=u(1,:);
    row_err=zeros(N+1,1);
    %Error=zeros(N,1);
    abserr=zeros(N+1,L+1);

end       



%calculate the order ratio for each value of N, except the first
N=4;
Ri = zeros(N+1,1);
Ri(1) = Error(1);
for k=2:N+1
    Ri(k) = (1/log(2))*(log(Error(k - 1) / Error(k)));
end