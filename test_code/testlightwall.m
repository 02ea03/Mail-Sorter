clc; clear all; close all;
camlight;
view(3)
axis equal;
hold on;
workspace=[-10 10 -10 10 -0.1 10];
%% option 1
%setting up the light curtain
n=10;
posLight = zeros(n,9);
for i = 1:n
    x1(i,:) = 1;
    x2(i,:) = 1;
    x3(i,:) = -1;
    x4(i,:) = -1;    
    y1(i,:) = 1;
    y2(i,:) = -1;
    y3(i,:) = -1;
    y4(i,:) = 1;    
    z(i,:) = i;
    
    posLight(i,1) = x1(i,:);
    posLight(i,2) = x2(i,:);
    posLight(i,3) = x3(i,:);
    posLight(i,4) = x4(i,:);
    posLight(i,5) = y1(i,:);
    posLight(i,6) = y2(i,:);
    posLight(i,7) = y3(i,:);
    posLight(i,8) = y4(i,:);
    posLight(i,9) = z(i,:);
end
%plotting

for i = 1:n
    xT = [posLight(i,1), posLight(i,2), posLight(i,3), posLight(i,4),posLight(i,1)];
    yT = [posLight(i,5),posLight(i,6),posLight(i,7),posLight(i,8),posLight(i,5)];
    zT = [posLight(i,9),posLight(i,9),posLight(i,9),posLight(i,9),posLight(i,9)];
    plot3([xT],[yT],[zT]);
end

%% Collision
boxLight = collisionBox(2,2,10); %(length x,width y,hight z)
ballLight = collisionSphere(0.5); %(radius)
transBall = trvec2tform([10 1 0]); %translating (transl)
transBox = trvec2tform([0 0 5]); %translating (transl)
ballBox.Pose = transBox;
ballLight.Pose = transBall;
% [~,patchObj] = show(boxLight); % if i never show this it still exists but you cannot see it
hold on;
axis([-5 15 -5 5 -1 10]);
tempBall = show(ballLight);
n = 10;
for i = 0:0.5:n
    newX=n-i;
    newY=1-i/n;
    newZ=2*sin(i/4);
    transBall = trvec2tform([newX newY newZ]);
    ballLight.Pose = transBall;
    tempBall = show(ballLight);
    pause(0.1);
    [areIntersecting,dist,witnessPoints] = checkCollision(boxLight,ballLight); %collision check
    if areIntersecting == 1
        display('intersection robot stopping')
        pause(); %robot has to stop
        for i = 0:0.5:n %bouncing back / ball no longer intersecting with the light curtain
            newXup=newX+i;
            newYup=-i/n;
            newZup=newZ+2*sin(i/4);
            transBall = trvec2tform([newXup newYup newZup]);
            ballLight.Pose = transBall;
            tempBall = show(ballLight);
            pause(0.1);
        end
        break;
    end
end
pause();
delete(tempBall);




