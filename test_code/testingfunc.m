clc;clear all; close all;
%%
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
%%
% Create image target (points in the image plane) 
pStar = [662 362 362 662; 362 362 662 662];
%Create 3D points
P=[1.8,1.8,1.8,1.8;
-0.25,0.25,0.25,-0.25;
 1.25,1.25,0.75,0.75];
%initalise the camera
cam = CentralCamera('focal', 0.08, 'pixel', 10e-5, ...
'resolution', [1024 1024], 'centre', [512 512],'name', 'camera');
fps = 30; %frame rate
lambda = 0.6; %gain of the controler
depth = mean (P(1,:)); %depth of the IBVS
Tc0= robotDenso.fkine(q0);
cam.plot_camera('Tcam',Tc0, 'label','scale',0.5);
%%
% for i = 1:100-1 %RMRC
%     qMatrix = 0; %(i,:) whatever the column is initalising as i dont have initial qMatrix 
%     qAdd = 0.1 %4x4 matrix adding either x y or z
%     while qMatrix < 5
%         qMatrix = qMatrix + qAdd
%     end
% end