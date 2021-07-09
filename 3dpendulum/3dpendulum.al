% 3D pendulum
% Starting point for a model of a dog's tail which uses a 3D 
% double pendulum with specified motion of the base
%
%   Tom Rottier, June 2021
%--------------------------------------------------------------------
%   Settings
overwrite on
degrees off
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
q4 = 0
dircos(n,a,body123,q4,q5,q6)
%--------------------------------------------------------------------
%   Position vectors
p_o_p1> = q1*n1> + q2*n2> + q3*n3>  % From origin to base of pendulum
p_p1_ao> = lao*a1>                  % From base to CoM
p_p1_p2> = la*a1>                   % From base to tip
%--------------------------------------------------------------------
%   Motion constraints
auxiliary[1] = u1       % Constrains velocity of p1 to zero
auxiliary[2] = u2
auxiliary[3] = u3
dependent[1] = u4

constrain( auxiliary[u1,u2,u3], dependent[u4] )
%--------------------------------------------------------------------
%    Angular velocity of pendulum
w_a_n> = u4*a1> + u5*a2> + u6*a3>
%--------------------------------------------------------------------
%   Kinematical differential equations
q1' = u1
q2' = u2
q3' = u3
kindiffs(n,a,body123,q4,q5,q6)   
%--------------------------------------------------------------------
%   Velocities
v_o_n> = 0>
v_p1_n> = dt(p_o_p1>, n)
v_ao_n> = dt(p_o_ao>, n)
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
amomx = dot(angmom>, n1>)
amomy = dot(angmom>, n2>)
amomz = dot(angmom>, n3>)
%--------------------------------------------------------------------
%   Desired quantities to output
p1x = dot(p_o_p1>, n1>)
p1y = dot(p_o_p1>, n2>)
p1z = dot(p_O_p1>, n3>)
p2x = dot(p_o_p2>, n1>)
p2y = dot(p_o_p2>, n2>)
p2z = dot(p_O_p2>, n3>)
%--------------------------------------------------------------------
%   Set up system of non-linear equations to solve for Q's and U's
vp2x = dot( dt(p_o_p2>, n), n1>)
vp2y = dot( dt(p_o_p2>, n), n2>)
vp2z = dot( dt(p_o_p2>, n), n3>)

out[1] = vp2x - rhs(vp2x)
out[2] = vp2y - rhs(vp2y)
out2[1] = vp2y - rhs(vp2y)
out2[2] = vp2z - rhs(vp2z)

%evaluate(u5,u1=0,u2=0,u3=0,q5=1,q6=1,la=1,vp2x=1,vp2y=2,vp2z=0)
%evaluate(u6,u1=0,u2=0,u3=0,q5=1,q6=1,la=1,vp2x=1,vp2y=2,vp2z=0)
%--------------------------------------------------------------------
%   Outputs
output t,p1x,p1y,p1z,p2x,p2y,p2z
output t,rx,ry,rz,ke,pe,te,amomx,amomy,amomz
%--------------------------------------------------------------------
%   Input constants
input abserr=1.0e-08, relerr=1.0e-07
%--------------------------------------------------------------------
%   generate code
code dynamics() 3dpendulum.f, nosubs


