%SIR Please help me with the meanings of the variables used like
%w_max,t1,t2,t2,a1,a2,a3
tic
clear all
clc
close all

% VT Graph parameters
speed_rpm=5
w_max=speed_rpm*6*pi/180

t1=0.5
t2=11.8
t3=12

a1=w_max/t1
a2=0
a3=w_max/(t3-t2)
torque=zeros(2);

t=0:0.1:t3;

for i=1:size(t,2)
    if t(i)<=t1;
        vel(i)=a1*t(i);
        theta(i)=a1*0.5*t(i)^2;
        acc(i)=a1;
    elseif (t(i)<=t2)&(t(i)>t1);
        vel(i)=w_max;
        theta(i)=w_max*(t(i)-t1)+(a1*0.5*t1*t1);
        acc(i)=0;
    else
        vel(i)=w_max-a3*(t(i)-t2);
         %theta(i)=(w_max*(t(i)-t2))-(a3*0.5*(t(i)-t2)^2)+(a3*t2*(t(i)-t2));
        %theta(i)=(2*pi)-(0.5*a3*(t3-t(i))^2);
        s=t(i);
        theta(i)=w_max*(s-t2)+(a3/6)*(-s*s*s+t2^3-3*s*t2^2+3*t2*s^2)+(w_max*(t2-t1)+(a1*0.5*t1*t1));
        acc(i)=-a3;
    end
end


L1=0.27;
L2=0.27;
m1=0.2;
m2=0.2;
deg=linspace(0,360,size(t,2));
for i=1:size(t,2)
 [torque(1,i),torque(2,i)]=RR_manipulator_dynamics(L1,L2,deg(i),m1,m2,vel(i),acc(i));
end
figure (1)
plot(t,theta)
hold on
figure (2)
plot(t,vel)
hold on
figure (3)
plot(t,acc)
hold on; 
figure (4)
plot(t,torque(1,:),'-m')
figure (5)
plot(t,torque(2,:),'-m')

toc
