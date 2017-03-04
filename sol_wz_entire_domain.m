a = 1;
b = 2;

%initialize parameters
t0 = 0;
tf = 1;
L = 4;
 

%define space mesh
h = 1 / L;
x = 0 : h / 2 : 1;
x = x';

%set inital funcitons  
InitialPressure = sin( pi * x );
InitialVelocity = ones( 1 , length( x ) );

% define time mesh
k = h / 4 ;
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
RI = -1/det(R)*R;
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
C=struct( 'w1',zeros( N + 1 , L + 1 ),...
          'z1',zeros( N + 1 , L + 1 ),...
          'w2',zeros( N + 1 , L + 1 ),...
          'z2',zeros( N + 1 , L + 1 ));
W= (-sin(pi*(x)) + 1)';
Z =  (sin(pi*(x)) + 1)';
C.w1( 1 , : ) = W( 1 , 1 : L + 1 );
C.z1( 1 , : ) = Z( 1 , 1 : L + 1 );
C.w2( 1 , : ) = W( 1 , L + 1 : end);
C.z2( 1 , : ) = Z( 1 , L + 1 : end );
 
 
%struct to hold the exaxt solution
E =struct('ExactPressureEdge1',zeros( N + 1 , L + 1 ),...
          'ExactVelocityEdge1',zeros( N + 1 , L + 1 ),...
          'ExactPressureEdge2',zeros( N + 1 , L + 1 ),...
          'ExactVelocityEdge2',zeros( N + 1 , L + 1 ));
      
%loop through time
for n=1:N
    
   %calculate interior nodes for edge 1
   W( n + 1 , : ) = -sin( pi * ( x - c0w * t( n + 1 )))+1;
   Z( n + 1 , : ) = sin( pi * ( x - c0z * t( n + 1 )))+1;
   
   C.w1( n+1 , : ) = W( n+1 , 1 : L + 1 );
   C.z1( n+1 , : ) = Z( n+1 , 1 : L + 1 );
   C.w2( n+1 , : ) = W( n+1 , L + 1 : end);
   C.z2( n+1 , : ) = Z( n+1 , L + 1 : end );
   
   newEdge1 = RI * [ C.w1(n,:) ; C.z1(n,:) ];
   newEdge2 = RI * [ C.w2(n,:) ; C.z2(n,:) ];
   
   
   U.PressureEdge1(n,:) = newEdge1( 1, : );
   U.VelocityEdge1(n,:) = newEdge1( 2, : );
   U.PressureEdge2(n,:) = newEdge2( 1, : );
   U.VelocityEdge2(n,:) = newEdge2( 2, : );
   
   
   %calculate exact solution using known solutions for w,z
   wX = x - c0w * t( n );
   zX = x - c0z * t( n );
   wexact = ( -sin( pi * wX) + 1 )';
   zexact = ( sin ( pi * zX ) + 1 )';
   uexact= RI* [ wexact ; zexact ];
    
   E.ExactPressureEdge1( n , : ) =  uexact( 1 , 1 : L + 1 );
   E.ExactVelocityEdge1( n , : ) =  uexact( 2 , 1 : L + 1 );
   E.ExactPressureEdge2( n , : ) =  uexact( 1 , L + 1 : end );
   E.ExactVelocityEdge2( n , : ) =  uexact( 2 , L + 1 : end ); 
   
   

    
end

errorPE1 = max( max ( E.ExactPressureEdge1 - U.PressureEdge1 ));
errorVE1 = max( max ( E.ExactVelocityEdge1 - U.VelocityEdge1 ));
errorPE2 = max( max ( E.ExactPressureEdge2 - U.PressureEdge2 ));
errorVE2 = max( max ( E.ExactVelocityEdge2 - U.VelocityEdge2 ));



clf
%figure(glf)
% hold on
for i=1:N+1
    
    plot( x ( 1 : L + 1 ),U.PressureEdge1( i , : ),'b',...
          x ( L + 1 : end ),U.PressureEdge2( i , : ),'g',...
          x ( 1 : L + 1 ),U.VelocityEdge1( i , : ),'r',...
          x ( L + 1 : end ),U.VelocityEdge2( i , : ),'m',...
          x ( 1 : L + 1 ) , E.ExactPressureEdge1( i , : ) , 'b',...
          x ( L + 1 : end ) , E.ExactPressureEdge2( i , : ) , 'g',...
          x ( 1 : L + 1 ) , E.ExactVelocityEdge1( i , : ) , 'r',...
          x ( L + 1 : end ) , E.ExactVelocityEdge2( i , : ) , 'm')
      
	axis( [ 0  1  -2 2 ])
    pause( .05 )
    drawnow

end
hold off
    
    

