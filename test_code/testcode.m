clc
clear all
close all

L1 = Link('d',0.5,'a',0,'alpha',pi/2,'offset',0,'qlim',deg2rad([-160 160]));
L2 = Link('d',0,'a',0.84,'alpha',0,'offset',pi/2,'qlim',deg2rad([-40 120]));
L3 = Link('d',0,'a',0.3,'alpha',-pi/2,'offset',-pi/2,'qlim',deg2rad([30 160]));
L4 = Link('d',0.84,'a',0,'alpha',pi/2,'offset',0,'qlim',deg2rad([-160 160]));
L5 = Link('d',0,'a',0,'alpha',-pi/2,'offset',0,'qlim',deg2rad([-120 120]));
L6 = Link('d',0.28,'a',0,'alpha',0,'offset',0,'qlim',deg2rad([-360 360]));

VP6242 = SerialLink([L1 L2 L3 L4 L5 L6],'name','VP6242');

workspace = [-5 5 -5 5 -5 5];   
scale = 0.5; 
VP6242.teach()
% q = (deg2rad([30,30,35,-30,0,0]));  
% VP6242.plot(q,'workspace',workspace,'scale',scale);
%EEpose = robot.fkine(deg2rad(q))
% Brick1T=transl(0.1,0.2,0.1);
% initialq = [0 0 0 0 0 0];
% Startq = VP6242.getpos;
% brickPose = troty(pi);
% Brickq= VP6242.ikcon(Brick1T, initialq); % Pick up the brick
% steps = 100;
% s = lspb(0,1,steps);
% JointLink = nan(steps,VP6242.n);
% initialq = VP6242.getpos
%      % from start position to brick
% for j = 1:steps
%     JointLink(j,:) = (1-s(j))*initialq + s(j)*Brickq;
%     VP6242.animate(JointLink(j,:))
%     pause(0.01)
%    if j == steps
%         disp(Brickq) %brick location
%         disp(UR3.model.fkine(UR3.model.getpos))%end pos
%         %disp(self.brick1Transform(1:3,4))%location
%         pause
%    end
% %    if j == steps
% %         disp(Brickq) %brick location
% %         disp(UR3.model.fkine(UR3.model.getpos))%end pos
% %         %disp(self.brick1Transform(1:3,4))%location
% %         pause
% %    end
% end
% BrickPos=VP6242.getpos
% Brickq2= VP6242.ikcon(brickPose,BrickPos); % Pick up the brick
% s = lspb(0,1,steps);
% JointLink = nan(steps,VP6242.n);
% currentq = VP6242.getpos;
%      % from start position to brick
% for j = 1:steps
%     JointLink(j,:) = (1-s(j))*currentq + s(j)*Brickq;
%     VP6242.animate(JointLink(j,:))
%     pause(0.01)
% end

