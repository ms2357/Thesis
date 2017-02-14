t = [ 0  .25 0.5000  .75 1.0000];
c0w=-1;
c0z=1;
N=length(t);
L=length(x);
w=zeros(N,2);
z=zeros(N,2);

for n=1:N;

    w(n,1) =-sin(pi*(0-c0w*t(n)))+1;
    z(n,1) = sin(pi*(0-c0z*t(n)))+1;
    
    w(n,2) =-sin(pi*(1-c0w*t(n)))+1;
    z(n,2) = sin(pi*(1-c0z*t(n)))+1;
end