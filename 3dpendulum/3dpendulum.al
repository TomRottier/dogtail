% 3D pendulum
% Starting point for a model of a dog's tail which uses a 3D 
% double pendulum with specified motion of the base
%
%   Tom Rottier, June 2021
%--------------------------------------------------------------------
%   Settings

%--------------------------------------------------------------------
%   Physical declarations
newtonian n
bodies a
points o, p1 ,p2
mass a=m
inertia a,ix,iy,iz
%--------------------------------------------------------------------
%   Mathematical declarations
variables q{6}'
variables u{6}'
variables rx,ry,rz
constants g,lao,la
%--------------------------------------------------------------------
%   Geometry relating unit vectors
dircos(n,a,space123,q4,q5,q6)
%--------------------------------------------------------------------
%   Position vectors
p_o_p1> = q1*n1> + q2*n2> + q3*n3>  % From origin to base of pendulum
p_p1_ao> = lao*a1>                  % From base to CoM
p_p1_p2> = la*a1>                   % From base to tip

%--------------------------------------------------------------------
%    Angular velocity of pendulum
w_a_n> = u4*a1> + u5*a2> + u6*a3>
%--------------------------------------------------------------------
%   Kinematical differential equations
q1' = u1
q2' = u2
q3' = u3

kindiffs(n,a,space123,q4,q5,q6)   
%--------------------------------------------------------------------
%   Velocities
v_o_n> = 0>
v_p1_n> = dt(p_o_p1>, n)
v_ao_n> = dt(p_o_ao>, n)
%--------------------------------------------------------------------
%   Motion constraints
auxiliary[1] = u1       % Constrains velocity of p1 to zero
auxiliary[2] = u2
auxiliary[3] = u3

constrain( auxiliary[u1,u2,u3] )
%--------------------------------------------------------------------
%   Angular accelerations
alf_a_n> = dt(w_a_n>, n)
%--------------------------------------------------------------------
%   Accelerations
a_o_n> = 0>
a_p1_n> = dt(v_p1_n>, n)
a_ao_n> = dt(v_ao_n>, n)
%--------------------------------------------------------------------
%   Forces
gravity(g*n3>)
force(p1, rx*n1> + ry*n2> + rz*n3>)
%--------------------------------------------------------------------
%   Equations of motion
zero = fr() + frstar()
kane(rx,ry,rz)
%--------------------------------------------------------------------
%   Energy and momentum
ke = ke()                   % Kinetic energy
pe = -mass(a)*g * dot(p_o_ao>, n3>)   % Potential energy
te = ke + pe
angmom> = momentum(angular, p1)
%--------------------------------------------------------------------
%   Desired quantities
p1x = dot(p_o_p1>, n1>)
p1y = dot(p_o_p1>, n2>)
p1z = dot(p_O_p1>, n3>)
p2x = dot(p_o_p2>, n1>)
p2y = dot(p_o_p2>, n2>)
p2z = dot(p_O_p2>, n3>)
%--------------------------------------------------------------------
%   outputs
output t,p1x,p1y,p1z,p2x,p2y,p2z
output t,rx,ry,rz,ke,pe,te
%--------------------------------------------------------------------
%   Input constants
input abserr=1.0e-08, relerr=1.0e-07
%--------------------------------------------------------------------
%   generate code
code dynamics() 3dpendulum.f, nosubs


