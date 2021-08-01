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
autoz on
%--------------------------------------------------------------------
%   Physical declarations
newtonian n
bodies a,b
frames f
points o, p1, p2, p3, cm
mass a=ma, b=mb
inertia a,ixa,iya,iza
inertia b,ixb,iyb,izb
%--------------------------------------------------------------------
%   Mathematical declarations
variables q{6}'
variables u{6}'
variables rx,ry,rz,atorx,atory,atorz,btorx,btory,btorz
constants lao,la,lbo,lb
constants ka,ba,kb,bb       % Stiffness and damping for each spring
constants eqx,eqy,eqz
constants g
specified ox'', oy'', oz''  % Orientation of frame f about n
specified px'', py'', pz''  % Position and derivatives of base
%--------------------------------------------------------------------
%   Geometry relating unit vectors
dircos(n, f, body123, ox, oy, oz)
dircos(f, a, body123, q1, q2, q3)
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
constrain( dependent[u4] )
%--------------------------------------------------------------------
%    Angular velocity of bodies
angvel(n, f, body123, ox, oy, oz)
w_a_f> = u1*a1> + u2*a2> + u3*a3>
w_b_a> = u4*b1> + u5*b2> + u6*b3>
%--------------------------------------------------------------------
%   Kinematical differential equations
kindiffs(f, a, body123, q1, q2, q3)
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
alf_f_n> = dt(w_f_n>, n)
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
atorx = eqx*eqy*eqz*ka*ba*t
atory = ka*ba*t
atorz = ka*ba*t
btorx = kb*bb*t
btory = kb*bb*t
btorz = kb*bb*t
torque(a, atorx*a1> + atory*a2> + atorz*a3>)
torque(a/b, btorx*b1> + btory*b2> + btorz*b3>)
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
%   Outputs
output t,p1x,p1y,p1z,p2x,p2y,p2z,p3x,p3y,p3z
output t,ke,pe,te,amomx,amomy,amomz
%--------------------------------------------------------------------
%   Input constants
input abserr=1.0e-08, relerr=1.0e-07
%--------------------------------------------------------------------
%   generate code
code dynamics() model/dogtail.f, nosubs