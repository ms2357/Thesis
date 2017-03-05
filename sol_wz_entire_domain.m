a = 1;
b = 2;

%initialize parameters
t0 = 0;
tf = 1;
L = 2;
 

%define space mesh
h = 1 / L;
x = 0 : h / 2 : 1;
x = x';

%set inital funcitons  
InitialPressure = ( sin( pi * x ) )';
InitialVelocity = ones( 1 , length( x ) );

% define time mesh
k = h / 2 ;
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

A = [ 0        K1 ; 
      1 / rho1 0 ];       %wave speed
  
R = [ -Z1 Z1 ;
       1   1 ];
RI = (1/(-Z1*1-1*Z1))*[ 1   -Z1;
                       -1  -Z1];
evaluesA = eig( A );
c0w = evaluesA( 1 );
c0z = evaluesA( 2 );

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
      
IC = RI* [ InitialPressure ; InitialVelocity ];
W= (-.5*sin(pi*(x)) + .5)';
Z =  (.5*sin(pi*(x)) + .5)';
C.w1( 1 , : ) = IC( 1 , 1 : L + 1 );
C.z1( 1 , : ) = IC( 2 , 1 : L + 1 );
C.w2( 1 , : ) = IC( 1 , L + 1 : end);
C.z2( 1 , : ) = IC( 2 , L + 1 : end );
 
 
%struct to hold the exaxt solution
E =struct('ExactPressureEdge1',zeros( N + 1 , L + 1 ),...
          'ExactVelocityEdge1',zeros( N + 1 , L + 1 ),...
          'ExactPressureEdge2',zeros( N + 1 , L + 1 ),...
          'ExactVelocityEdge2',zeros( N + 1 , L + 1 ));
      
%loop through time
for n=1:N+1
    
   %calculate interior nodes for edge 1
   W( n  , : ) = -.5*sin( pi * ( x - c0w * t( n ))) + .5;
   Z( n  , : ) =  .5*sin( pi * ( x - c0z * t( n ))) + .5;
   
   C.w1( n , : ) = W( n , 1 : L + 1 );
   C.z1( n , : ) = Z( n , 1 : L + 1 );
   C.w2( n , : ) = W( n , L + 1 : end);
   C.z2( n , : ) = Z( n , L + 1 : end );
   
   newEdge1 = R * [ C.w1(n,:) ; C.z1(n,:) ];
   newEdge2 = R * [ C.w2(n,:) ; C.z2(n,:) ];
   
   
   U.PressureEdge1(n,:) = newEdge1( 1, : );
   U.VelocityEdge1(n,:) = newEdge1( 2, : );
   U.PressureEdge2(n,:) = newEdge2( 1, : );
   U.VelocityEdge2(n,:) = newEdge2( 2, : );
end
for n=1:N+1   
   %calculate exact solution using known solutions for w,z
   wX = x - c0w * t( n );
   zX = x - c0z * t( n );
   wexact = ( -.5*sin( pi * wX) + .5 )';
   zexact = ( .5*sin ( pi * zX ) + .5 )';
   uexact= R* [ wexact ; zexact ];
    
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
          x ( 1 : L + 1 ) , E.ExactPressureEdge1( i , : ) , '.',...
          x ( L + 1 : end ) , E.ExactPressureEdge2( i , : ) , '.',...
          x ( 1 : L + 1 ) , E.ExactVelocityEdge1( i , : ) , '.',...
          x ( L + 1 : end ) , E.ExactVelocityEdge2( i , : ) , '.')
      
	axis( [ 0  1  -2 2 ])
    pause( .05 )
    drawnow

end
hold off
    
    

