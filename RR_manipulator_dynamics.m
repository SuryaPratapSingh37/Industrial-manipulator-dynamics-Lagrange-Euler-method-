% code for dynamic model of RR manipulator
function [torque1,torque2] = RR_manipulator_dynamics(L1,L2,deg,m1,m2,vel,acc);

%syms c1 s1 c2 s2 g m1 m2 v1 v2 L1 L2 alph1 alph2 TAU t1 t2 I1  I2 

%Input values
% t1 is theta 1 in degrees
% t2 is theta 2 in degrees
% m1 is mass of link 1
% m2 is mass of link 2
% v1 is angular velocity of 1st link
% v2 is angular velocity of 2nd link
% alph1 is angular acceleration of 1st link
% alph2 is angular acceleration of 2nd link
% L1 is length of 1st link
% L2 is length of 2nd link

v1=vel;
v2=vel;
alph1=acc;
alph2=acc;
t1=deg;
t2=deg;

g=9.81;%m/s^2 acceleration of gravity
c1=cosd(t1);
c2=cosd(t2);
s1=sind(t1);
s2=sind(t2);

%Transformation matrices
T00=eye(4)
T01=sym([c1,-s1,0,L1*c1;s1,c1,0,L1*s1;0,0,1,0;0,0,0,1])
T12=sym([c2,-s2,0,L2*c2;s2,c2,0,L2*s2;0,0,1,0;0,0,0,1])
T02=T01*T12

Q1=[0,-1,0,0;1,0,0,0;0,0,0,0;0,0,0,0] % Rotary joint
Q2=Q1 % Because 2nd joint also rotary

% I1=[Ix1 0 0 0;0 Iy1  0 0;0 0 Iz1 0;0 0 0 1]
% I2=[Ix2 0 0 0;0 Iy2  0 0;0 0 Iz2 0;0 0 0 1]

% I1=[0 0 0 0;0 0 0 0;0 0 0 0;0 0 0 m]
% I2=[0 0 0 0;0 0 0 0;0 0 0 0;0 0 0 m]

% I1=[(1/3)*m1*L^2 0 0 -m1*L/2;0 0  0 0;0 0 0 0;-m1*L/2 0 0 m1]
% I2=[(1/3)*m2*L^2 0 0 -m2*L/2;0 0  0 0;0 0 0 0;-m2*L/2 0 0 m2]

I1=[(1/3)*m1*L1^2 0 0 -m1*L1/2;0 0  0 0;0 0 0 0;-m1*L1/2 0 0 m1] %Inertia matrix of 1st link
I2=[(1/3)*m2*L2^2 0 0 -m2*L2/2;0 0  0 0;0 0 0 0;-m2*L2/2 0 0 m2] %Inertia matrix of 2nd link

d11=T00*Q1*T01
d12=0
d21=T00*Q1*T02
d22=T01*Q2*T12

%Inertia matrix coefficients
M11=trace(d11*I1*transpose(d11))+trace(d21*I2*transpose(d21))
M22=trace(d22*I2*transpose(d22))
M12=trace(d22*I2*transpose(d21))
M21=trace(d21*I2*transpose(d22))
M=[M11 M12;M21 M22]

simplify(M)

 % Velocity coupling coefficients
 h111=0;
 h112=trace(T00*Q1*T01*Q2*T12*I2*transpose(d21))
 h121=trace(Q1*T01*Q2*T12*I2*transpose(d21))
 h122=trace(T01*Q2*Q2*T12*I2*transpose(d21))
 h211=trace(T00*Q1*Q1*T02*I2*transpose(d22))
 h212=trace(T00*Q1*T01*Q2*T12*I2*transpose(d22))
 h221=trace(Q1*T01*Q2*T12*I2*transpose(d22))
 h222=0;
 
 H1=h111*v1^2+h112*v1*v2+h121*v1*v2+h122*v2^2
 H2=h211*v1^2+h212*v1*v2+h221*v1*v2+h222*v2^2
 
 H=[H1; H2]
 
 r1=[-L1/2; 0; 0; 1]
 r2=[-L2/2; 0; 0; 1]
 
 gr=[0 -g 0 0]
 
 %Gravity loading
 G1=-(m1*gr*d11*r1+m2*gr*d21*r2)
 G2=-m2*gr*d22*r2
 
 G=[G1; G2]
 
 
 TAU=M*[alph1 ;alph2]+H+G;
 torque1=TAU(1)
 torque2=TAU(2)
end

 