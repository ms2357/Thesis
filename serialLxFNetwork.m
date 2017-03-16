a = 1;
b = 2;

%initialize parameters
t0 = 0;
tf = 1;
L = 8;
 

%define space mesh
h = 1 / L;
x = 0 : h  : 1;
x = x';

%set inital funcitons  
InitialPressure = sin( pi * x );
InitialVelocity = ones( 1 , length( x ) );

% define time mesh
k = h  ;
t = t0 : k : tf;
N = length( t ) - 1;

%wave speed details
rho1 = 1;                 %density
K1 = 1;                   %bulk modulus
rho2 = 1;
K2 = 1;
%impedence
Z1 = rho1 * K1;             
Z2 = rho2 * K2;

A = [ 0 K1 ; 1 / rho1 0 ];       %wave speed
R = [ -Z1 Z1 ; 1 1 ];
evaluesA = eig( A );
c0w = evaluesA( 1 );
c0z = evaluesA( 2 );

%CFL number
mu = A * k / h;

%preallocate u, set initial value using intital condition
U =struct('PressureEdge1',zeros( N + 1 , L + 1 ),...
          'VelocityEdge1',zeros( N + 1 , L + 1 ),...
          'PressureEdge2',zeros( N + 1 , L + 1 ),...
          'VelocityEdge2',zeros( N + 1 , L + 1 ));
U( : ).PressureEdge1( 1 , : ) = InitialPressure( 1 : L + 1 );
U( : ).VelocityEdge1( 1 , : ) = InitialVelocity( 1 : L + 1 );
U( : ).PressureEdge2( 1 , : ) = InitialPressure( L + 1 : end );
U( : ).VelocityEdge2( 1 , : ) = InitialVelocity( L + 1 : end );



%struct to store w,z at junctions, index gives egde, column 1 of struct is
%the char at the left endpoint 'a' and column 2 gives the char at the 
% right endpoint 'b'
%w left moving 
%z right moving
C=struct( 'w1',zeros( N + 1 , 2 ),...
          'z1',zeros( N + 1 , 2 ),...
          'w2',zeros( N + 1 , 2 ),...
          'z2',zeros( N + 1 , 2 ) );

%loop through time
for n=1:N
    %loop through space
    for j=2:L
%u(n + 1,j)= (u(n,j - 1) + u(n,j + 1))/2 -(mu/2)*(u(n,j + 1) - u(n,j - 1))
      
      %calculate interior nodes for edge 1
      prev = ( cell2mat( { U.PressureEdge1( n , j - 1 ),...
                           U.VelocityEdge1( n , j - 1 ) } ) )';
                  
      next = ( cell2mat( { U.PressureEdge1( n , j + 1 ),...
                           U.VelocityEdge1( n , j + 1 ) } ) )';
                  
      newU = ( prev + next ) / 2 - ( mu / 2 ) * ( next - prev );

      U.PressureEdge1( n + 1 , j ) = newU( 1 );
      U.VelocityEdge1( n + 1 , j ) = newU( 2 );
      
      %calculate interior nodes for edge 2
      prev = ( cell2mat( { U.PressureEdge2( n , j - 1),...
                           U.VelocityEdge2( n , j - 1 ) } ) )';
                    
      next = ( cell2mat( { U.PressureEdge2( n , j + 1 ),...
                           U.VelocityEdge2( n , j + 1 ) } ) )';
                    
      newU = ( prev + next ) / 2 - ( mu / 2 ) * ( next - prev );

      U.PressureEdge2( n + 1 , j ) = newU( 1 );
      U.VelocityEdge2( n + 1 , j ) = newU( 2 );
      

    end

    %JUNCTIONS
    
    j = 1;

    %calc  w char using pres/vel for egde 1 at vertex a
    prev = R * ( cell2mat( { U.PressureEdge1( n , j ),...
                             U.VelocityEdge1( n , j ) } ) )';
    
    %check to ensure press/vel for w/z satisfy IC
    C.w1( 1 , a ) = prev( 1 );
    C.z1( 1 , a ) = prev( 2 );
    
    next = R *( cell2mat( { U.PressureEdge1( n , j + 1 ),...
                             U.VelocityEdge1( n , j + 1 ) } ) )';
   
    %w(x^n)=c0(dt/dx)*(w_(L)-w_(L-1))+w_L...linear interpolation for known char
    %where c0=-1, an eigen value of A and the slope of the characteristic
    curr = c0w * ( k / h ) * ( prev - next ) + prev;
    
    %set w char value at vertex  
    C.w1( n + 1 , a ) = curr( 1 );
        
    %calc w char using pres/vel for edge 2 at vertex b
    prev = R * ( cell2mat( { U.PressureEdge2( n , j ),...
                             U.VelocityEdge2( n , j ) } ) )';
      
    next = R * ( cell2mat( { U.PressureEdge2( n , j + 1 ),...
                             U.VelocityEdge2( n , j + 1 ) } ) )';
   
    %check to ensure press/vel for w/z satisfy IC
    C.w2( 1 , b ) = prev( 1 );
    C.z2( 1 , b ) = prev( 2 );
      
    %w(x^n)=c0(dt/dx)*(z_(L)-z_(L-1))+z_L...linear interpolation for known char
    %where c0=-1, an eigen value of A and the slope of the characteristic
    curr = c0w * ( k / h ) * ( prev - next ) + prev;
    
    %set w char value at vertex 
    C.w2( n + 1 , b ) = curr( 1 );
    
    
    %JUNCTIONS
    
    j = L + 1;
    
    %calc z char using pres/vel for edge 1 at vertex b 
    prev = R * ( cell2mat( { U.PressureEdge1( n , j - 1 ),...
                             U.VelocityEdge1( n , j - 1 ) } ) )';
   
    next = R * (cell2mat( { U.PressureEdge1( n , j ),...
                            U.VelocityEdge1( n , j ) } ) )';
                      
    %check to ensure press/vel for w/z satisfy IC
    C.w1( 1 , b ) = next( 1 );
    C.z1( 1 , b ) = next( 2 );
    
    %z(x^n)=mu*(z_(L)-z_(L-1))+z_L...linear interpolation for known char
    curr = c0z * ( k / h ) * ( prev - next ) + next;
    
    %set w char value at vertex 
    C.z1( n + 1 , b ) = curr( 2 );
    
    %calc chars using pres/vel for edge 2 at vertex a
    prev = R * ( cell2mat( { U.PressureEdge2( n , j - 1 ),...
                             U.VelocityEdge2( n , j - 1 ) } ) )';
                         
    next = R * ( cell2mat( { U.PressureEdge2( n , j ),...
                             U.VelocityEdge2( n , j ) } ) )';
                         
    %check to ensure press/vel for w/z satisfy IC
    C.w2( 1 , a ) = next( 1 );
    C.z2( 1 , a ) = next( 2 );
   
    %z(x^n)=mu*(z_(L)-z_(L-1))+z_L...linear interpolation for known char
    curr = c0z * ( k / h ) * ( next - prev ) + next;
    
    %set w char value at vertex 
    C.z2( n + 1 , a ) = curr( 2 );
    
    
    
    %Vertex a
    %calc unk chars at vertex using juction conditions. Solve system 
    %Ax=b using interpolated chars 
    knownchars = ( cell2mat( { C.w1( n + 1 , a ),...
                               C.z2( n + 1 , a ) } ) )';
    %multiply by coeff matrix                       
    knownchars = [ Z1  Z2 ; -1 1 ] * knownchars;
    %solve
    unkchars = [ Z2 Z1 ; -1 1 ] \ knownchars ;
    %set new chars at vertex a 
    C.w2( n + 1 , a ) = unkchars( 1 );
    C.z1( n + 1 , a ) = unkchars( 2 );
  
    %calc press/vel using chars on edge 1 at vertex a
    newC = ( cell2mat( { C.w1( n + 1 , a ),...
                         C.z1( n + 1 , a ) } ) )';
                     
    newU = R \ newC;
    U.PressureEdge1( n + 1 , 1 ) = newU( 1 );
    U.VelocityEdge1( n + 1 , 1 ) = newU( 2 );
    
    %calc press/vel using chars on edge 2 at vertex a
    newC = ( cell2mat( { C.w2( n + 1 , a ),...
                         C.z2( n + 1 , a ) } ) )';
                     
    newU = R \ newC;
    U.PressureEdge2( n + 1 , L + 1 ) = newU( 1 );
    U.VelocityEdge2( n + 1 , L + 1 ) = newU( 2 );
    
    
    %Vertex b
    %calc unk chars at vertex using juction conditions. Solve system 
    %Ax=b using interpolated chars 
    knownchars = ( cell2mat( { C.w2( n + 1 , b ),...
                               C.z1( n + 1 , b ) } ) )';
    
    %multiply by coefficient matrix
    knownchars = [ Z2 Z1 ; 1 -1 ] * knownchars;
    %solve
    unkchars = [ Z1 Z2 ; 1 -1 ] \ knownchars;
    C.w1( n + 1 , b ) = unkchars( 1 );
    C.z2( n + 1 , b ) = unkchars( 2 );
    
    
    %calc press/vel using chars on edge 1 at vertex b  
    newC = ( cell2mat( { C.w1( n + 1 , b ),...
                         C.z1( n + 1 , b ) } ) )';
                     
    newU = R \ newC;
    U.PressureEdge1( n + 1 , L + 1 ) = newU( 1 );
    U.VelocityEdge1( n + 1 , L + 1 ) = newU( 2 );
   
    %calc press/vel using chars on edge 2 at vertex b
    newC = ( cell2mat( { C.w2( n + 1 , b ),...
                         C.z2( n + 1 , b ) } ) )';
                     
    newU = R \ newC;
    U.PressureEdge2( n + 1 , 1 ) = newU( 1 );
    U.VelocityEdge2( n + 1 , 1 ) = newU( 2 );

end

clf
%figure(glf)
% hold on
for i=1:N+1
    
    plot( x ( 1 : L + 1 ),U.PressureEdge1( i , : ),'b',...
          x ( L + 1 : end ),U.PressureEdge2( i , : ),'g',...
          x ( 1 : L + 1 ),U.VelocityEdge1( i , : ),'r',...
          x ( L + 1 : end ),U.VelocityEdge2( i , : ),'m')
      
	axis( [ 0  1  0 1.5 ])
    pause( .05 )
    drawnow

end
hold off
