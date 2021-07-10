% Simple model of dog tail
% A double pendulum with the motion of the base specified.
% Spring-dampers connect the rigid segments.
%
%   Tom Rottier, June 2021
%--------------------------------------------------------------------
%   Settings
overwrite on
degrees off
autorhs off
%--------------------------------------------------------------------
%   Physical declarations
newtonian n
bodies a,b
points o, p1, p2, p3, cm
mass a=ma, b=mb
inertia a,ixa,iya,iza
inertia b,ixb,iyb,izb
%--------------------------------------------------------------------
%   Mathematical declarations
variables q{6}'
variables u{6}'
variables rx,ry,rz
constants lao,la,lbo,lb
constants ka,ba,kb,bb       % Stiffness and damping for each spring
constants g
specified px'', py'', pz''  % Position and derivatives of base
%--------------------------------------------------------------------
%   Geometry relating unit vectors
dircos(n, a, body123, q1, q2, q3)
dircos(a, b, body123, q4, q5, q6)
%--------------------------------------------------------------------
%   Position vectors
p_o_p1> = px*n1> + py*n2> + pz*n3>  % From origin to base of first seg
p_p1_ao> = lao*a1>                  % From base to CoM of first seg 
p_p1_p2> = la*a1>                   % From base of first to base of second seg
p_p2_bo> = lbo*b1>  
p_p2_p3> = lb*b1>   
%--------------------------------------------------------------------
%   Motion constraints
dependent[1] = u4
constrain(dependent[u4])
%--------------------------------------------------------------------
%    Angular velocity of pendulum
w_a_n> = u1*a1> + u2*a2> + u3*a3>
w_b_a> = u4*b1> + u5*b2> + u6*b3>
%--------------------------------------------------------------------
%   Kinematical differential equations
kindiffs(n, a, body123, q1, q2, q3)
kindiffs(a, b, body123, q4, q5, q6)
%--------------------------------------------------------------------
%   Velocities
v_o_n> = 0>
v_p1_n> = dt(p_o_p1>, n)
v_ao_n> = dt(p_o_ao>, n)
v_p2_n> = dt(p_o_p2>, n)
v_bo_n> = dt(p_o_bo>, n)
%--------------------------------------------------------------------
%   Angular accelerations
alf_a_n> = dt(w_a_n>, n)
alf_b_n> = dt(w_b_n>, n)
%--------------------------------------------------------------------
%   Accelerations
a_o_n> = 0>
a_p1_n> = dt(v_p1_n>, n)
a_ao_n> = dt(v_ao_n>, n)
a_p2_n> = dt(v_p2_n>, n)
a_bo_n> = dt(v_bo_n>, n)
%--------------------------------------------------------------------
%   Forces and torques
gravity(g*n3>)

ator> = -ka*(q1*a1> + q2*a2> + q3*a3>) - ba*w_a_n>
btor> = -kb*(q4*b1> + q5*b2> + q6*b3>) - bb*w_b_a>
torque(a, ator>)
torque(a/b, btor>)
%--------------------------------------------------------------------
%   Equations of motion
zero = fr() + frstar()
kane()
%--------------------------------------------------------------------
%   Energy and momentum
ke = ke()                   % Kinetic energy
pe = -mass()*g * (dot(cm(o), n3>) + (lao+lbo))  % Potential energy
te = ke + pe
angmom> = momentum(angular, p1)
amomx = dot(angmom>, n1>)
amomy = dot(angmom>, n2>)
amomz = dot(angmom>, n3>)
%--------------------------------------------------------------------
%   Desired quantities to output
p1x = dot(p_o_p1>, n1>)
p1y = dot(p_o_p1>, n2>)
p1z = dot(p_o_p1>, n3>)
p2x = dot(p_o_p2>, n1>)
p2y = dot(p_o_p2>, n2>)
p2z = dot(p_o_p2>, n3>)
p3x = dot(p_o_p3>, n1>)
p3y = dot(p_o_p3>, n2>)
p3z = dot(p_o_p3>, n3>)
%--------------------------------------------------------------------
%   Set up system of non-linear equations to solve for Q's and U's
vp2x = dot( dt(p_o_p2>, n), n1>)
vp2y = dot( dt(p_o_p2>, n), n2>)
vp2z = dot( dt(p_o_p2>, n), n3>)
vp3x = dot( dt(p_o_p3>, n), n1>)
vp3y = dot( dt(p_o_p3>, n), n2>)
vp3z = dot( dt(p_o_p3>, n), n3>)

out[1] = vp2x - rhs(vp2x)
out[2] = vp2y - rhs(vp2y)
out[3] = vp2z - rhs(vp2z)
out[4] = vp3x - rhs(vp3x)
out[5] = vp3y - rhs(vp3y)
out[6] = vp3z - rhs(vp3z)
%--------------------------------------------------------------------
%   Outputs
output t,p1x,p1y,p1z,p2x,p2y,p2z,p3x,p3y,p3z
output t,ke,pe,te,amomx,amomy,amomz
%--------------------------------------------------------------------
%   Input constants
input abserr=1.0e-08, relerr=1.0e-07
%--------------------------------------------------------------------
%   generate code
%code nonlinear(out,u1,u2,u3,u4,u5,u6) invkin.f
code dynamics() dogtail.f, nosubs
