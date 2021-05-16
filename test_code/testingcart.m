clc; close all; clear all;
%%
qMatrix = zeros(1,6);
L1 = Link('d',0.5,'a',0,'alpha',pi/2,'offset',0,'qlim',deg2rad([-360 360]));
L2 = Link('d',0,'a',0.84,'alpha',0,'offset',pi/2,'qlim',deg2rad([-40 120]));
L3 = Link('d',0,'a',0.3,'alpha',-pi/2,'offset',0,'qlim',deg2rad([-140 90]));
L4 = Link('d',0.84,'a',0,'alpha',pi/2,'offset',0,'qlim',deg2rad([-180 180]));
L5 = Link('d',0,'a',0,'alpha',-pi/2,'offset',0,'qlim',deg2rad([-120 120]));
L6 = Link('d',0.28,'a',0,'alpha',0,'offset',0,'qlim',deg2rad([-360 360]));
robotDenso = SerialLink([L1 L2 L3 L4 L5 L6],'name','VP6242');
q0 = deg2rad([0 0 0 0 0 0]);
scale=0.6;
workspace=[-10 10 -10 10 -0.1 10];
robotDenso.base = transl(0,0,3.2);
view(-15,21);
robotDenso.plot(q0,'workspace',workspace,'scale',scale);
x(1,1) =  1.5;% Points in x input
x(2,1) =  1.5;% Points in y input
x(3,1) =  1.5;   % Points in z input
theta(1,1) = -176*pi/180                % Roll angle
theta(2,1) = 1*pi/180             % Pitch angle
theta(3,1) = -0.01*pi/18 %Yaw
T = [rpy2r(theta(1,1),theta(2,1),theta(3,1)) x(:,1);zeros(1,3) 1];          % Create transformation of first point and angle
q0 = zeros(1,6);                                                          % Initial guess for joint angles
qMatrix = robotDenso.ikcon(T,q0)