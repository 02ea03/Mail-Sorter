clf;clc;clear all;close all;
%doCameraSpin = false;
%requires optimization tool box add on
% Turn on a light (only turn on 1, don't keep turning them on), and make axis equal
camlight;
view(3)
axis equal;
hold on;

%% e-stop
location  = [-2,-3.55,3.1];        % location of proposed safety infrastructure % mounting position for fences,desk
[f,v,data] = plyread('button01.ply','tri');
vertexColours = [data.vertex.red, data.vertex.green, data.vertex.blue] / 255;

for zOffset = location(1,3)
    for yOffset = location(1,2)
        for xOffset = location(1,1)
         trisurf(f,v(:,1) + xOffset,v(:,2) + yOffset, v(:,3) + zOffset ...
        ,'FaceVertexCData',vertexColours,'EdgeColor','interp','EdgeLighting','flat');
        end
    end
end
hold on;
%% fire_extinguisher
    [f,v,data] = plyread('fire_extinguisher01.ply','tri');
    vertexColours = [data.vertex.red, data.vertex.green, data.vertex.blue] / 255;

    for xOffset = [-5]
        for yOffset = [4]
            trisurf(f,v(:,1)+ xOffset,v(:,2) + yOffset, v(:,3) ...
                ,'FaceVertexCData',vertexColours,'EdgeColor','interp','EdgeLighting','flat');
        end
    end
  
  
%% box pink

    [f,v,data] = plyread('boxpink.ply','tri');
    vertexColours = [data.vertex.red, data.vertex.green, data.vertex.blue] / 255;

    for xOffset = [-1]
        for yOffset = [-1.3]
                 for zOffset = [3.1]
           trisurf(f,v(:,1) + xOffset,v(:,2) + yOffset, v(:,3) + zOffset ...
                ,'FaceVertexCData',vertexColours,'EdgeColor','interp','EdgeLighting','flat');
                 end
        end
    end
%% box Blue
    [f,v,data] = plyread('boxblue.ply','tri');
    vertexColours = [data.vertex.red, data.vertex.green, data.vertex.blue] / 255;

    for xOffset = [-1]
        for yOffset = [1.3]
                 for zOffset = [3.1]
           trisurf(f,v(:,1) + xOffset,v(:,2) + yOffset, v(:,3) + zOffset ...
                ,'FaceVertexCData',vertexColours,'EdgeColor','interp','EdgeLighting','flat');
                 end
        end
    end
    
    
    %% Belt Conveyor
    [f,v,data] = plyread('belt01.ply','tri');
    vertexColours = [data.vertex.red, data.vertex.green, data.vertex.blue] / 255;

    for xOffset = [1.5]
        for yOffset = [0]
                 for zOffset = [3.1]
           trisurf(f,v(:,1) + xOffset,v(:,2) + yOffset, v(:,3) + zOffset ...
                ,'FaceVertexCData',vertexColours,'EdgeColor','interp','EdgeLighting','flat');
                 end
        end
    end
    %% wall alarm
    [f,v,data] = plyread('wallalarm.ply','tri');
    vertexColours = [data.vertex.red, data.vertex.green, data.vertex.blue] / 255;
    for xOffset = [7.8]
        for yOffset = [0]
                 for zOffset = [3]
           trisurf(f,v(:,1) + xOffset,v(:,2) + yOffset, v(:,3) + zOffset ...
                ,'FaceVertexCData',vertexColours,'EdgeColor','interp','EdgeLighting','flat');
                 end
        end
    end
     %% first aid
    [f,v,data] = plyread('firstaidbox.ply','tri');
    vertexColours = [data.vertex.red, data.vertex.green, data.vertex.blue] / 255;
    for xOffset = [0]
        for yOffset = [7.8]
                 for zOffset = [3]
           trisurf(f,v(:,1) + xOffset,v(:,2) + yOffset, v(:,3) + zOffset ...
                ,'FaceVertexCData',vertexColours,'EdgeColor','interp','EdgeLighting','flat');
                 end
        end
    end
    %% red fence
    [f,v,data] = plyread('redfence.ply','tri');
    vertexColours = [data.vertex.red, data.vertex.green, data.vertex.blue] / 255;

    for xOffset = [-6.5]
        for yOffset = [3]
            trisurf(f,v(:,1)+ xOffset,v(:,2) + yOffset, v(:,3) ...
                ,'FaceVertexCData',vertexColours,'EdgeColor','interp','EdgeLighting','flat');
        end
    end
    
         %% wall 2 

    [fWall,vWall,data] = plyread('Wall.ply','tri');
    vertexColours = [data.vertex.red, data.vertex.green, data.vertex.blue] / 255;

    for xOffset = [0.3]
        for yOffset = [2.1]
                 for zOffset = [3.1]
           wallV([1]) = trisurf(fWall,vWall(:,1) + xOffset,vWall(:,2) + yOffset, vWall(:,3) + zOffset ...
                ,'FaceVertexCData',vertexColours,'EdgeColor','interp','EdgeLighting','flat');
                 end
        end
    end
        faceNormalsWall = zeros(size(fWall,1),3);
    for faceIndexWall = 1:size(fWall,1)
        vWall1 = vWall(fWall(faceIndexWall,1)',:);
        vWall2 = vWall(fWall(faceIndexWall,2)',:);
        vWall3 = vWall(fWall(faceIndexWall,3)',:);
        faceNormalsWall(faceIndexWall,:) = unit(cross(vWall2-vWall1,vWall3-vWall1));
    end
    
        %% red frence 02
    [f,v,data] = plyread('redfence.ply','tri');
    vertexColours = [data.vertex.red, data.vertex.green, data.vertex.blue] / 255;

    for xOffset = [-6.5]
        for yOffset = [-2]
            trisurf(f,v(:,1)+ xOffset,v(:,2) + yOffset, v(:,3) ...
                ,'FaceVertexCData',vertexColours,'EdgeColor','interp','EdgeLighting','flat');
        end
    end
    
    
        %% robotDenso human
    [f,v,data] = plyread('human01.ply','tri');
    vertexColours = [data.vertex.red, data.vertex.green, data.vertex.blue] / 255;

    for xOffset = [-6.3]
        for yOffset = [-5.3]
            trisurf(f,v(:,1)+ xOffset,v(:,2) + yOffset, v(:,3) ...
                ,'FaceVertexCData',vertexColours,'EdgeColor','interp','EdgeLighting','flat');
        end
    end
    
%% Table
    [f,v,data] = plyread('table01.ply','tri');
    vertexColours = [data.vertex.red, data.vertex.green, data.vertex.blue] / 255;

    for xOffset = [0]
        for yOffset = [0]
            trisurf(f,v(:,1)+ xOffset,v(:,2) + yOffset, v(:,3) ...
                ,'FaceVertexCData',vertexColours,'EdgeColor','interp','EdgeLighting','flat');
        end
    end

%% Wall Stick 
surf([-8 8; -8 8],... %X
[8 8; 8 8],... %Y
[8 8; 0 0],... %Z (height in Z;plot in x)
'CData',imread('robot_arm.jpg'),'FaceColor','texturemap');
hold on;

surf([8 8; 8 8],... %X
[-7 8; -7 8],... %Y
[8 8; 0 0],... %Z (height in Z;plot in x)
'CData',imread('bg06.jpg'),'FaceColor','texturemap');
hold on;

surf([-8,-8;8,8],... %X
[-7,8;-7,8],... %Y
[0,0;0,0],... %Z (height in Z;plot in x)
'CData',imread('green.jpg'),'FaceColor','texturemap');
hold on;
%% Papers
[f,v,data] = plyread('paperCircle.ply','tri');
vertexColours = [data.vertex.red, data.vertex.green, data.vertex.blue] / 255;
for zOffset = [3.5]
    for yOffset = [0]
        for xOffset = [1.6]            
            paperCircle([1])= trisurf(f,v(:,1) + xOffset,v(:,2) + yOffset, v(:,3) + zOffset ...
        ,'FaceVertexCData',vertexColours,'EdgeColor','interp','EdgeLighting','flat');
        end
    end
end

[f,v,data] = plyread('paperCircle.ply','tri');
vertexColours = [data.vertex.red, data.vertex.green, data.vertex.blue] / 255;
for zOffset = [3.53]
    for yOffset = [0]
        for xOffset = [1.6]            
            paperCircle([2])= trisurf(f,v(:,1) + xOffset,v(:,2) + yOffset, v(:,3) + zOffset ...
        ,'FaceVertexCData',vertexColours,'EdgeColor','interp','EdgeLighting','flat');
        end
    end
end

[f,v,data] = plyread('paperTri.ply','tri');
vertexColours = [data.vertex.red, data.vertex.green, data.vertex.blue] / 255;
for zOffset = [3.51]
    for yOffset = [0]
        for xOffset = [1.6]          
            paperTri([3]) = trisurf(f,v(:,1) + xOffset,v(:,2) + yOffset, v(:,3) + zOffset ... %delete(paperTri([1]) to delete the
        ,'FaceVertexCData',vertexColours,'EdgeColor','interp','EdgeLighting','flat');
        end
    end
end

[f,v,data] = plyread('paperTri.ply','tri');
vertexColours = [data.vertex.red, data.vertex.green, data.vertex.blue] / 255;
for zOffset = [3.52]
    for yOffset = [0]
        for xOffset = [1.6]          
            paperTri([4]) = trisurf(f,v(:,1) + xOffset,v(:,2) + yOffset, v(:,3) + zOffset ... %delete(paperTri([1]) to delete the
        ,'FaceVertexCData',vertexColours,'EdgeColor','interp','EdgeLighting','flat');
        end
    end
end
display('Finished Loading Background')
%% Robot
%robotDenso = VP6242(false);
%q = deg2rad([0 0 0 0 0 0]);
%robotDenso.model.base = transl(0,0,3.2);
view(-15,21);
%robotDenso.model.animate(q);
L1 = Link('d',0.5,'a',0,'alpha',pi/2,'offset',0,'qlim',deg2rad([-360 360]));
L2 = Link('d',0,'a',0.84,'alpha',0,'offset',pi/2,'qlim',deg2rad([-40 120]));
L3 = Link('d',0,'a',0.3,'alpha',-pi/2,'offset',0,'qlim',deg2rad([-140 90]));
L4 = Link('d',0.84,'a',0,'alpha',pi/2,'offset',0,'qlim',deg2rad([-180 180]));
L5 = Link('d',0,'a',0,'alpha',-pi/2,'offset',0,'qlim',deg2rad([-120 120]));
L6 = Link('d',0.28,'a',0,'alpha',0,'offset',0,'qlim',deg2rad([-360 360]));
robotDenso = SerialLink([L1 L2 L3 L4 L5 L6],'name','VP6242');
q = deg2rad([0 0 0 0 0 0]);
scale=0.6;
workspace=[-10 10 -10 10 -0.1 10];
robotDenso.base = transl(0,0,3.2);
view(-15,21);
robotDenso.plot(q,'workspace',workspace,'scale',scale);
hold on;
pause(0.01);
view(-15,21);
%% Standing Position to Mail
% 1.1) Set parameters for the simulation
t = 10;             % Total time (s)
deltaT = 0.1;      % Control frequency
steps = t/deltaT;   % No. of steps for simulation
delta = 2*pi/steps; % Small angle change
epsilon = 0.1;      % Threshold value for manipulability/Damped Least Squares
W = diag([1 1 1 0.1 0.1 0.1]);    % Weighting matrix for the velocity vector

% 1.2) Allocate array data
m = zeros(steps,1);             % Array for Measure of Manipulability
qMatrix = zeros(steps,6);       % Array for joint angles
qdot = zeros(steps,6);          % Array for joint velocities
theta = zeros(3,steps);         % Array for roll-pitch-yaw angles
x = zeros(3,steps);             % Array for x-y-z trajectory
positionError = zeros(3,steps); % For plotting trajectory error
angleError = zeros(3,steps);    % For plotting trajectory error

% 1.3) Set up trajectory, initial pose
s = lspb(0,1,steps);                % Trapezoidal trajectory scalar

for i=1:steps
    x(1,i) = (1-s(i))*-1.12 + s(i)*1.6  ; % Points in x
    x(2,i) = (1-s(i))*0 + 1.5*sin(i*delta/2); % Points in y
    x(3,i) = (1-s(i))*4.84 +sin(i*delta/2)+ s(i)*4; % Points in z
    theta(1,i) = (1-s(i))*0 + s(i)*-176*pi/180;                 % Roll angle
    theta(2,i) = (1-s(i))*-pi/2 + s(i)*1*pi/180;            % Pitch angle
    theta(3,i) = (1-s(i))*0 + s(i)*-0.01*pi/18;            % Yaw angle
end

T = [rpy2r(theta(1,1),theta(2,1),theta(3,1)) x(:,1);zeros(1,3) 1];          % Create transformation of first point and angle
q0 = zeros(1,6);                                                            % Initial guess for joint angles
qMatrix(1,:) = robotDenso.ikcon(T,q0);                                      % Solve joint angles to achieve first waypoint

% 1.4) Track the trajectory with RMRC
for i = 1:steps-1
    T = robotDenso.fkine(qMatrix(i,:));                                     % Get forward transformation at current joint state
    deltaX = x(:,i+1) - T(1:3,4);                                         	% Get position error from next waypoint
    Rd = rpy2r(theta(1,i+1),theta(2,i+1),theta(3,i+1));                     % Get next RPY angles, convert to rotation matrix
    Ra = T(1:3,1:3);                                                        % Current end-effector rotation matrix
    Rdot = (1/deltaT)*(Rd - Ra);                                                % Calculate rotation matrix error
    S = Rdot*Ra';                                                           % Skew symmetric!
    linear_velocity = (1/deltaT)*deltaX;
    angular_velocity = [S(3,2);S(1,3);S(2,1)];                              % Check the structure of Skew Symmetric matrix!!
    deltaTheta = tr2rpy(Rd*Ra');                                            % Convert rotation matrix to RPY angles
    xdot = W*[linear_velocity;angular_velocity];                          	% Calculate end-effector velocity to reach next waypoint.
    J = robotDenso.jacob0(qMatrix(i,:));                 % Get Jacobian at current joint state
    m(i) = sqrt(round(det(J*J')));
    if m(i) < epsilon  % If manipulability is less than given threshold
        lambda = (1 - m(i)/epsilon)*5E-2;
    else
        lambda = 0;
    end
    invJ = inv(J'*J + lambda *eye(6))*J';                                   % DLS Inverse
    qdot(i,:) = (invJ*xdot)';                                                % Solve the RMRC equation (you may need to transpose the         vector)
    for j = 1:6                                                             % Loop through joints 1 to 6
        if qMatrix(i,j) + deltaT*qdot(i,j) < robotDenso.qlim(j,1)                     % If next joint angle is lower than joint limit...
            qdot(i,j) = 0; % Stop the motor
        elseif qMatrix(i,j) + deltaT*qdot(i,j) > robotDenso.qlim(j,2)                 % If next joint angle is greater than joint limit ...
            qdot(i,j) = 0; % Stop the motor
        end
    end
    qMatrix(i+1,:) = qMatrix(i,:) + deltaT*qdot(i,:);                         	% Update next joint state based on joint velocities
    positionError(:,i) = x(:,i+1) - T(1:3,4);                               % For plotting
    angleError(:,i) = deltaTheta;                                           % For plotting
end

% Collsion Checking and Re-mapping if needed
while(1)
    faces = fWall;
    vertex = vWall;
    faceNormals = faceNormalsWall;
    goUp = [0,0,0,0;0,0,0,0;0,0,0,0.1;0,0,0,0];

    if IsCollision(robotDenso,qMatrix(i,:),faces,vertex,faceNormals,false)
    disp('Collision detected!!');
    newRoute = robotDenso.fkine(qMatrix)+ goUp;
    newRoute(:,:,1) = robotDenso.fkine(qMatrix)- goUp;
    newRoute(:,:,100) = robotDenso.fkine(qMatrix)- goUp;
    qMatrix = robotDenso.ikcon(newRoute);
  
    else
    words=['No collision found in the trajectory!'];
    disp(words); 
    break;
    end
    
end

% Plot the results
for i=1:steps-1
    robotDenso.animate(qMatrix(i,:)); %plot the robot
end
pause(1);
%% Mail to Pink Box
%setting the ply file (for some reason it reads the last ply file read)
[f,v,data] = plyread('paperCircle.ply','tri');
vertexColours = [data.vertex.red, data.vertex.green, data.vertex.blue] / 255;
% 1.1) Set parameters for the simulation
t = 10;             % Total time (s)
deltaT = 0.1;      % Control frequency
steps = t/deltaT;   % No. of steps for simulation
delta = 2*pi/steps; % Small angle change
epsilon = 0.1;      % Threshold value for manipulability/Damped Least Squares
W = diag([1 1 1 0.1 0.1 0.1]);    % Weighting matrix for the velocity vector

% 1.2) Allocate array data
m = zeros(steps,1);             % Array for Measure of Manipulability
qMatrix = zeros(steps,6);       % Array for joint angles
qdot = zeros(steps,6);          % Array for joint velocities
theta = zeros(3,steps);         % Array for roll-pitch-yaw angles
x = zeros(3,steps);             % Array for x-y-z trajectory
positionError = zeros(3,steps); % For plotting trajectory error
angleError = zeros(3,steps);    % For plotting trajectory error

% 1.3) Set up trajectory, initial pose
s = lspb(0,1,steps);                % Trapezoidal trajectory scalar

for i=1:steps
    x(1,i) = (1-s(i))*1.6 + s(i)*-0.9; % Points in x
    x(2,i) = (1-s(i))*0.01 - s(i)*1.1; % Points in y
    x(3,i) = (1-s(i))*4 + s(i)*4.2;     % Points in z
    theta(1,i) = (1-s(i))*-176*pi/180 + s(i)*-176*pi/180;                 % Roll angle 
    theta(2,i) = (1-s(i))*1*pi/180 + s(i)*0.29*pi/180;            % Pitch angle
    theta(3,i) = (1-s(i))*-0.01*pi/18 + s(i)*80.48*pi/180; %yaw angle
end
T = [rpy2r(theta(1,1),theta(2,1),theta(3,1)) x(:,1);zeros(1,3) 1];          % Create transformation of first point and angle
q0 = zeros(1,6);                                                            % Initial guess for joint angles
qMatrix(1,:) = robotDenso.ikcon(T,q0);                                            % Solve joint angles to achieve first waypoint

% 1.4) Track the trajectory with RMRC
for i = 1:steps-1
    T = robotDenso.fkine(qMatrix(i,:));                                           % Get forward transformation at current joint state
    deltaX = x(:,i+1) - T(1:3,4);                                         	% Get position error from next waypoint
    Rd = rpy2r(theta(1,i+1),theta(2,i+1),theta(3,i+1));                     % Get next RPY angles, convert to rotation matrix
    Ra = T(1:3,1:3);                                                        % Current end-effector rotation matrix
    Rdot = (1/deltaT)*(Rd - Ra);                                                % Calculate rotation matrix error
    S = Rdot*Ra';                                                           % Skew symmetric!
    linear_velocity = (1/deltaT)*deltaX;
    angular_velocity = [S(3,2);S(1,3);S(2,1)];                              % Check the structure of Skew Symmetric matrix!!
    deltaTheta = tr2rpy(Rd*Ra');                                            % Convert rotation matrix to RPY angles
    xdot = W*[linear_velocity;angular_velocity];                          	% Calculate end-effector velocity to reach next waypoint.
    J = robotDenso.jacob0(qMatrix(i,:));                 % Get Jacobian at current joint state
    m(i) = sqrt(round(det(J*J')));
    if m(i) < epsilon  % If manipulability is less than given threshold
        lambda = (1 - m(i)/epsilon)*5E-2;
    else
        lambda = 0;
    end
    invJ = inv(J'*J + lambda *eye(6))*J';                                   % DLS Inverse
    qdot(i,:) = (invJ*xdot)';                                                % Solve the RMRC equation (you may need to transpose the         vector)
    for j = 1:6                                                             % Loop through joints 1 to 6
        if qMatrix(i,j) + deltaT*qdot(i,j) < robotDenso.qlim(j,1)           % If next joint angle is lower than joint limit...
            qdot(i,j) = 0; % Stop the motor
        elseif qMatrix(i,j) + deltaT*qdot(i,j) > robotDenso.qlim(j,2)       % If next joint angle is greater than joint limit ...
            qdot(i,j) = 0; % Stop the motor
        end
    end
    qMatrix(i+1,:) = qMatrix(i,:) + deltaT*qdot(i,:);                         	% Update next joint state based on joint velocities
    positionError(:,i) = x(:,i+1) - T(1:3,4);                               % For plotting
    angleError(:,i) = deltaTheta;                                           % For plotting
end

% Collsion Checking and Re-mapping if needed
while(1)
    faces = fWall;
    vertex = vWall;
    faceNormals = faceNormalsWall;
    goUp = [0,0,0,0;0,0,0,0;0,0,0,0.1;0,0,0,0];
    
    if IsCollision(robotDenso,qMatrix(i,:),faces,vertex,faceNormals,false)
    disp('Collision detected!!');
    newRoute = robotDenso.fkine(qMatrix)+ goUp;
    newRoute(:,:,1) = robotDenso.fkine(qMatrix)- goUp;
    newRoute(:,:,100) = robotDenso.fkine(qMatrix)- goUp;
    qMatrix = robotDenso.ikcon(newRoute);
  
    else
    words=['No collision found in the trajectory!'];
    disp(words); 
    break;
    end
    
end

% Plot the results
for i=1:steps-1
    delete(paperCircle([2]));   %Deleting paper
    posEE=robotDenso.fkine(qMatrix(i,:));   %gets the end-effector in terms of 4x4 matrix
    xOffset = posEE(1,4);   %x of end-effector
    yOffset = posEE(2,4);   %y of end-effector
    zOffset = posEE(3,4)-0.175; %z of end-effector and positioning paper to be under the gripper
    paperCircle([2])= trisurf(f,v(:,1) + xOffset,v(:,2) + yOffset, v(:,3) + zOffset ... 
        ,'FaceVertexCData',vertexColours,'EdgeColor','interp','EdgeLighting','flat');   %Replot the paper
    robotDenso.animate(qMatrix(i,:)); %plot the robot
end
pause(1);
delete(paperCircle([2]));
zOffset = [3.15];    %drop off (z coodinate change)
paperCircle([2])= trisurf(f,v(:,1) + xOffset,v(:,2) + yOffset, v(:,3) + zOffset ... 
        ,'FaceVertexCData',vertexColours,'EdgeColor','interp','EdgeLighting','flat'); %Replot the paper
pause(1);
%% Pink Box to Mail
% 1.1) Set parameters for the simulation
t = 10;             % Total time (s)
deltaT = 0.1;      % Control frequency
steps = t/deltaT;   % No. of steps for simulation
delta = 2*pi/steps; % Small angle change
epsilon = 0.1;      % Threshold value for manipulability/Damped Least Squares
W = diag([1 1 1 0.1 0.1 0.1]);    % Weighting matrix for the velocity vector

% 1.2) Allocate array data
m = zeros(steps,1);             % Array for Measure of Manipulability
qMatrix = zeros(steps,6);       % Array for joint angles
qdot = zeros(steps,6);          % Array for joint velocities
theta = zeros(3,steps);         % Array for roll-pitch-yaw angles
x = zeros(3,steps);             % Array for x-y-z trajectory
positionError = zeros(3,steps); % For plotting trajectory error
angleError = zeros(3,steps);    % For plotting trajectory error

% 1.3) Set up trajectory, initial pose
s = lspb(0,1,steps);                % Trapezoidal trajectory scalar

for i=1:steps
    x(1,i) = (1-s(i))*-0.9 + s(i)*1.6; % Points in x
    x(2,i) = (1-s(i))*-1.16 +s(i)*0; % Points in y
    x(3,i) = (1-s(i))*4.2 + s(i)*4;     % Points in z
    theta(1,i) = (1-s(i))*-176*pi/180 + s(i)*-176*pi/180; % Roll angle 
    theta(2,i) = (1-s(i))*0.29*pi/180 + s(i)*1*pi/180;  % Pitch angle
    theta(3,i) = (1-s(i))*80.48*pi/180 + s(i)*-0.01*pi/18; %yaw angle
end
T = [rpy2r(theta(1,1),theta(2,1),theta(3,1)) x(:,1);zeros(1,3) 1];          % Create transformation of first point and angle
q0 = zeros(1,6);                                                            % Initial guess for joint angles
qMatrix(1,:) = robotDenso.ikcon(T,q0);                                            % Solve joint angles to achieve first waypoint

% 1.4) Track the trajectory with RMRC
for i = 1:steps-1
    T = robotDenso.fkine(qMatrix(i,:));                                           % Get forward transformation at current joint state
    deltaX = x(:,i+1) - T(1:3,4);                                         	% Get position error from next waypoint
    Rd = rpy2r(theta(1,i+1),theta(2,i+1),theta(3,i+1));                     % Get next RPY angles, convert to rotation matrix
    Ra = T(1:3,1:3);                                                        % Current end-effector rotation matrix
    Rdot = (1/deltaT)*(Rd - Ra);                                                % Calculate rotation matrix error
    S = Rdot*Ra';                                                           % Skew symmetric!
    linear_velocity = (1/deltaT)*deltaX;
    angular_velocity = [S(3,2);S(1,3);S(2,1)];                              % Check the structure of Skew Symmetric matrix!!
    deltaTheta = tr2rpy(Rd*Ra');                                            % Convert rotation matrix to RPY angles
    xdot = W*[linear_velocity;angular_velocity];                          	% Calculate end-effector velocity to reach next waypoint.
    J = robotDenso.jacob0(qMatrix(i,:));                 % Get Jacobian at current joint state
    m(i) = sqrt(round(det(J*J')));
    if m(i) < epsilon  % If manipulability is less than given threshold
        lambda = (1 - m(i)/epsilon)*5E-2;
    else
        lambda = 0;
    end
    invJ = inv(J'*J + lambda *eye(6))*J';                                   % DLS Inverse
    qdot(i,:) = (invJ*xdot)';                                                % Solve the RMRC equation (you may need to transpose the         vector)
    for j = 1:6                                                             % Loop through joints 1 to 6
        if qMatrix(i,j) + deltaT*qdot(i,j) < robotDenso.qlim(j,1)                     % If next joint angle is lower than joint limit...
            qdot(i,j) = 0; % Stop the motor
        elseif qMatrix(i,j) + deltaT*qdot(i,j) > robotDenso.qlim(j,2)                 % If next joint angle is greater than joint limit ...
            qdot(i,j) = 0; % Stop the motor
        end
    end
    qMatrix(i+1,:) = qMatrix(i,:) + deltaT*qdot(i,:);                         	% Update next joint state based on joint velocities
    positionError(:,i) = x(:,i+1) - T(1:3,4);                               % For plotting
    angleError(:,i) = deltaTheta;                                           % For plotting
end

% Collsion Checking and Re-mapping if needed
while(1)
    faces = fWall;
    vertex = vWall;
    faceNormals = faceNormalsWall;
    goUp = [0,0,0,0;0,0,0,0;0,0,0,0.1;0,0,0,0];

    if IsCollision(robotDenso,qMatrix(i,:),faces,vertex,faceNormals,false)
    disp('Collision detected!!');
    newRoute = robotDenso.fkine(qMatrix)+ goUp;
    newRoute(:,:,1) = robotDenso.fkine(qMatrix)- goUp;
    newRoute(:,:,100) = robotDenso.fkine(qMatrix)- goUp;
    qMatrix = robotDenso.ikcon(newRoute);
  
    else
    words=['No collision found in the trajectory!'];
    disp(words); 
    break;
    end
    
end

% Plot the results
for i=1:steps-1
    robotDenso.animate(qMatrix(i,:)); % plotting the robot movements
end
pause(1);
%% Mail to Blue Box
%setting the ply file (for some reason it reads the last ply file read)
[f,v,data] = plyread('paperTri.ply','tri');
vertexColours = [data.vertex.red, data.vertex.green, data.vertex.blue] / 255;
% 1.1) Set parameters for the simulation
t = 10;             % Total time (s)
deltaT = 0.1;      % Control frequency
steps = t/deltaT;   % No. of steps for simulation
delta = 2*pi/steps; % Small angle change
epsilon = 0.1;      % Threshold value for manipulability/Damped Least Squares
W = diag([1 1 1 0.1 0.1 0.1]);    % Weighting matrix for the velocity vector

% 1.2) Allocate array data
m = zeros(steps,1);             % Array for Measure of Manipulability
qMatrix = zeros(steps,6);       % Array for joint angles
qdot = zeros(steps,6);          % Array for joint velocities
theta = zeros(3,steps);         % Array for roll-pitch-yaw angles
x = zeros(3,steps);             % Array for x-y-z trajectory
positionError = zeros(3,steps); % For plotting trajectory error
angleError = zeros(3,steps);    % For plotting trajectory error

% 1.3) Set up trajectory, initial pose
s = lspb(0,1,steps);                % Trapezoidal trajectory scalar

for i=1:steps
    x(1,i) = (1-s(i))*1.6 + s(i)*-0.9; % Points in x
    x(2,i) = (1-s(i))*0.01 + s(i)*1.16; % Points in y
    x(3,i) = (1-s(i))*4 + 0.5*sin(i*delta/2) + s(i)*4.2;     % Points in z
    theta(1,i) = (1-s(i))*-176*pi/180 + s(i)*-176*pi/180;   % Roll angle 
    theta(2,i) = (1-s(i))*1*pi/180 + s(i)*0.29*pi/180;      % Pitch angle
    theta(3,i) = (1-s(i))*-0.01*pi/18 + s(i)*80.48*pi/180;  %yaw angle
end
T = [rpy2r(theta(1,1),theta(2,1),theta(3,1)) x(:,1);zeros(1,3) 1];          % Create transformation of first point and angle
q0 = zeros(1,6);                                                            % Initial guess for joint angles
qMatrix(1,:) = robotDenso.ikcon(T,q0);                                            % Solve joint angles to achieve first waypoint

% 1.4) Track the trajectory with RMRC
for i = 1:steps-1
    T = robotDenso.fkine(qMatrix(i,:));                                           % Get forward transformation at current joint state
    deltaX = x(:,i+1) - T(1:3,4);                                         	% Get position error from next waypoint
    Rd = rpy2r(theta(1,i+1),theta(2,i+1),theta(3,i+1));                     % Get next Roll Pitch Yaw angles, convert to rotation matrix
    Ra = T(1:3,1:3);                                                        % Current end-effector rotation matrix
    Rdot = (1/deltaT)*(Rd - Ra);                                                % Calculate rotation matrix error
    S = Rdot*Ra';                                                           % Skew symmetric!
    linear_velocity = (1/deltaT)*deltaX;
    angular_velocity = [S(3,2);S(1,3);S(2,1)];                              % Check the structure of Skew Symmetric matrix!!
    deltaTheta = tr2rpy(Rd*Ra');                                            % Convert rotation matrix to RPY angles
    xdot = W*[linear_velocity;angular_velocity];                          	% Calculate end-effector velocity to reach next waypoint.
    J = robotDenso.jacob0(qMatrix(i,:));                 % Get Jacobian at current joint state
    m(i) = sqrt(round(det(J*J')));
    if m(i) < epsilon  % If manipulability is less than given threshold
        lambda = (1 - m(i)/epsilon)*5E-2;
    else
        lambda = 0;
    end
    invJ = inv(J'*J + lambda *eye(6))*J';                                   % DLS Inverse
    qdot(i,:) = (invJ*xdot)';                                                % Solve the RMRC equation (you may need to transpose the         vector)
    for j = 1:6                                                             % Loop through joints 1 to 6
        if qMatrix(i,j) + deltaT*qdot(i,j) < robotDenso.qlim(j,1)                     % If next joint angle is lower than joint limit...
            qdot(i,j) = 0; % Stop the motor
        elseif qMatrix(i,j) + deltaT*qdot(i,j) > robotDenso.qlim(j,2)                 % If next joint angle is greater than joint limit ...
            qdot(i,j) = 0; % Stop the motor
        end
    end
    qMatrix(i+1,:) = qMatrix(i,:) + deltaT*qdot(i,:);                         	% Update next joint state based on joint velocities
    positionError(:,i) = x(:,i+1) - T(1:3,4);                               % For plotting
    angleError(:,i) = deltaTheta;                                           % For plotting
end

% Collsion Checking and Re-mapping if needed
while(1)
    faces = fWall;
    vertex = vWall;
    faceNormals = faceNormalsWall;
    goUp = [0,0,0,0;0,0,0,0;0,0,0,0.1;0,0,0,0];

    if IsCollision(robotDenso,qMatrix(i,:),faces,vertex,faceNormals,false)
    disp('Collision detected!!');
    newRoute = robotDenso.fkine(qMatrix)+ goUp;
    newRoute(:,:,1) = robotDenso.fkine(qMatrix)- goUp;
    newRoute(:,:,100) = robotDenso.fkine(qMatrix)- goUp;

    qMatrix = robotDenso.ikcon(newRoute);
   
    else
    words=['No collision found in the trajectory!'];
    disp(words); 
    break;
    end
    
end
% Plot the results
for i=1:steps-1
    delete(paperTri([4]));
    posEE=robotDenso.fkine(qMatrix(i,:)); %gets the end-effector in terms of 4x4 matrix
    xOffset = posEE(1,4);   %x of end-effector
    yOffset = posEE(2,4);   %y of end-effector
    zOffset = posEE(3,4)-0.175; %z of end-effector and positioning paper to be under the gripper
    paperTri([4])= trisurf(f,v(:,1) + xOffset,v(:,2) + yOffset, v(:,3) + zOffset ... 
        ,'FaceVertexCData',vertexColours,'EdgeColor','interp','EdgeLighting','flat'); %replot paperTri
    robotDenso.animate(qMatrix(i,:)); %replot robot
end
pause(1);
delete(paperTri([4]));
zOffset = [3.15]; %setting the z to drop the paper off
paperTri([4])= trisurf(f,v(:,1) + xOffset,v(:,2) + yOffset, v(:,3) + zOffset ... 
        ,'FaceVertexCData',vertexColours,'EdgeColor','interp','EdgeLighting','flat'); %replor paper
pause(1);
%% Blue to Mail
% 1.1) Set parameters for the simulation
t = 10;             % Total time (s)
deltaT = 0.1;      % Control frequency
steps = t/deltaT;   % No. of steps for simulation
delta = 2*pi/steps; % Small angle change
epsilon = 0.1;      % Threshold value for manipulability/Damped Least Squares
W = diag([1 1 1 0.1 0.1 0.1]);    % Weighting matrix for the velocity vector

% 1.2) Allocate array data
m = zeros(steps,1);             % Array for Measure of Manipulability
qMatrix = zeros(steps,6);       % Array for joint angles
qdot = zeros(steps,6);          % Array for joint velocities
theta = zeros(3,steps);         % Array for roll-pitch-yaw angles
x = zeros(3,steps);             % Array for x-y-z trajectory
positionError = zeros(3,steps); % For plotting trajectory error
angleError = zeros(3,steps);    % For plotting trajectory error

% 1.3) Set up trajectory, initial pose
s = lspb(0,1,steps);                % Trapezoidal trajectory scalar

for i=1:steps
    x(1,i) = (1-s(i))*-0.9 + s(i)*1.6; % Points in x
    x(2,i) = (1-s(i))*1.16 +s(i)*0; % Points in y
    x(3,i) = (1-s(i))*4.2  +0.5*sin(i*delta/2)+ s(i)*4;     % Points in z
    theta(1,i) = (1-s(i))*-176*pi/180 + s(i)*-176*pi/180; % Roll angle 
    theta(2,i) = (1-s(i))*0.29*pi/180 + s(i)*1*pi/180;  % Pitch angle
    theta(3,i) = (1-s(i))*80.48*pi/180 + s(i)*-0.01*pi/18; %yaw angle
end
T = [rpy2r(theta(1,1),theta(2,1),theta(3,1)) x(:,1);zeros(1,3) 1];          % Create transformation of first point and angle
q0 = zeros(1,6);                                                            % Initial guess for joint angles
qMatrix(1,:) = robotDenso.ikcon(T,q0);                                            % Solve joint angles to achieve first waypoint

% 1.4) Track the trajectory with RMRC
for i = 1:steps-1
    T = robotDenso.fkine(qMatrix(i,:));                                           % Get forward transformation at current joint state
    deltaX = x(:,i+1) - T(1:3,4);                                         	% Get position error from next waypoint
    Rd = rpy2r(theta(1,i+1),theta(2,i+1),theta(3,i+1));                     % Get next RPY angles, convert to rotation matrix
    Ra = T(1:3,1:3);                                                        % Current end-effector rotation matrix
    Rdot = (1/deltaT)*(Rd - Ra);                                                % Calculate rotation matrix error
    S = Rdot*Ra';                                                           % Skew symmetric!
    linear_velocity = (1/deltaT)*deltaX;
    angular_velocity = [S(3,2);S(1,3);S(2,1)];                              % Check the structure of Skew Symmetric matrix!!
    deltaTheta = tr2rpy(Rd*Ra');                                            % Convert rotation matrix to RPY angles
    xdot = W*[linear_velocity;angular_velocity];                          	% Calculate end-effector velocity to reach next waypoint.
    J = robotDenso.jacob0(qMatrix(i,:));                 % Get Jacobian at current joint state
    m(i) = sqrt(round(det(J*J')));
    if m(i) < epsilon  % If manipulability is less than given threshold
        lambda = (1 - m(i)/epsilon)*5E-2;
    else
        lambda = 0;
    end
    invJ = inv(J'*J + lambda *eye(6))*J';                                   % DLS Inverse
    qdot(i,:) = (invJ*xdot)';                                                % Solve the RMRC equation (you may need to transpose the         vector)
    for j = 1:6                                                             % Loop through joints 1 to 6
        if qMatrix(i,j) + deltaT*qdot(i,j) < robotDenso.qlim(j,1)                     % If next joint angle is lower than joint limit...
            qdot(i,j) = 0; % Stop the motor
        elseif qMatrix(i,j) + deltaT*qdot(i,j) > robotDenso.qlim(j,2)                 % If next joint angle is greater than joint limit ...
            qdot(i,j) = 0; % Stop the motor
        end
    end
    qMatrix(i+1,:) = qMatrix(i,:) + deltaT*qdot(i,:);                         	% Update next joint state based on joint velocities
    positionError(:,i) = x(:,i+1) - T(1:3,4);                               % For plotting
    angleError(:,i) = deltaTheta;                                           % For plotting
end

% Collsion Checking and Re-mapping if needed
while(1)
    faces = fWall;
    vertex = vWall;
    faceNormals = faceNormalsWall;
    goUp = [0,0,0,0;0,0,0,0;0,0,0,0.1;0,0,0,0];

    if IsCollision(robotDenso,qMatrix(i,:),faces,vertex,faceNormals,false)
    disp('Collision detected!!');
    newRoute = robotDenso.fkine(qMatrix)+ goUp;
    newRoute(:,:,1) = robotDenso.fkine(qMatrix)- goUp;
    newRoute(:,:,100) = robotDenso.fkine(qMatrix)- goUp;

    qMatrix = robotDenso.ikcon(newRoute);
   
    else
    words=['No collision found in the trajectory!'];
    disp(words); 
    break;
    end
    
end

% Plot the results
for i=1:steps-1
    robotDenso.animate(qMatrix(i,:)); %(i,:) need for loop
end
pause(1);
%% Mail to Blue Box (second time)
%setting the ply file (for some reason it reads the last ply file read)
[f,v,data] = plyread('paperTri.ply','tri');
vertexColours = [data.vertex.red, data.vertex.green, data.vertex.blue] / 255;
% 1.1) Set parameters for the simulation
t = 10;             % Total time (s)
deltaT = 0.1;      % Control frequency
steps = t/deltaT;   % No. of steps for simulation
delta = 2*pi/steps; % Small angle change
epsilon = 0.1;      % Threshold value for manipulability/Damped Least Squares
W = diag([1 1 1 0.1 0.1 0.1]);    % Weighting matrix for the velocity vector

% 1.2) Allocate array data
m = zeros(steps,1);             % Array for Measure of Manipulability
qMatrix = zeros(steps,6);       % Array for joint angles
qdot = zeros(steps,6);          % Array for joint velocities
theta = zeros(3,steps);         % Array for roll-pitch-yaw angles
x = zeros(3,steps);             % Array for x-y-z trajectory
positionError = zeros(3,steps); % For plotting trajectory error
angleError = zeros(3,steps);    % For plotting trajectory error

% 1.3) Set up trajectory, initial pose
s = lspb(0,1,steps);                % Trapezoidal trajectory scalar

for i=1:steps
    x(1,i) = (1-s(i))*1.6 + s(i)*-0.9; % Points in x
    x(2,i) = (1-s(i))*0.01 + s(i)*1.16; % Points in y
    x(3,i) = (1-s(i))*4 +0.5*sin(i*delta/2)+ s(i)*4.2;     % Points in z
    theta(1,i) = (1-s(i))*-176*pi/180 + s(i)*-176*pi/180;   % Roll angle 
    theta(2,i) = (1-s(i))*1*pi/180 + s(i)*0.29*pi/180;      % Pitch angle
    theta(3,i) = (1-s(i))*-0.01*pi/18 + s(i)*80.48*pi/180;  %yaw angle
end
T = [rpy2r(theta(1,1),theta(2,1),theta(3,1)) x(:,1);zeros(1,3) 1];          % Create transformation of first point and angle
q0 = zeros(1,6);                                                            % Initial guess for joint angles
qMatrix(1,:) = robotDenso.ikcon(T,q0);                                            % Solve joint angles to achieve first waypoint

% 1.4) Track the trajectory with RMRC
for i = 1:steps-1
    T = robotDenso.fkine(qMatrix(i,:));                                           % Get forward transformation at current joint state
    deltaX = x(:,i+1) - T(1:3,4);                                         	% Get position error from next waypoint
    Rd = rpy2r(theta(1,i+1),theta(2,i+1),theta(3,i+1));                     % Get next RPY angles, convert to rotation matrix
    Ra = T(1:3,1:3);                                                        % Current end-effector rotation matrix
    Rdot = (1/deltaT)*(Rd - Ra);                                                % Calculate rotation matrix error
    S = Rdot*Ra';                                                           % Skew symmetric!
    linear_velocity = (1/deltaT)*deltaX;
    angular_velocity = [S(3,2);S(1,3);S(2,1)];                              % Check the structure of Skew Symmetric matrix!!
    deltaTheta = tr2rpy(Rd*Ra');                                            % Convert rotation matrix to RPY angles
    xdot = W*[linear_velocity;angular_velocity];                          	% Calculate end-effector velocity to reach next waypoint.
    J = robotDenso.jacob0(qMatrix(i,:));                 % Get Jacobian at current joint state
    m(i) = sqrt(round(det(J*J')));
    if m(i) < epsilon  % If manipulability is less than given threshold
        lambda = (1 - m(i)/epsilon)*5E-2;
    else
        lambda = 0;
    end
    invJ = inv(J'*J + lambda *eye(6))*J';                                   % DLS Inverse
    qdot(i,:) = (invJ*xdot)';                                                % Solve the RMRC equation (you may need to transpose the         vector)
    for j = 1:6                                                             % Loop through joints 1 to 6
        if qMatrix(i,j) + deltaT*qdot(i,j) < robotDenso.qlim(j,1)                     % If next joint angle is lower than joint limit...
            qdot(i,j) = 0; % Stop the motor
        elseif qMatrix(i,j) + deltaT*qdot(i,j) > robotDenso.qlim(j,2)                 % If next joint angle is greater than joint limit ...
            qdot(i,j) = 0; % Stop the motor
        end
    end
    qMatrix(i+1,:) = qMatrix(i,:) + deltaT*qdot(i,:);                         	% Update next joint state based on joint velocities
    positionError(:,i) = x(:,i+1) - T(1:3,4);                               % For plotting
    angleError(:,i) = deltaTheta;                                           % For plotting
end

% Collsion Checking and Re-mapping if needed
while(1)
    faces = fWall;
    vertex = vWall;
    faceNormals = faceNormalsWall;
    goUp = [0,0,0,0;0,0,0,0;0,0,0,0.1;0,0,0,0];

    if IsCollision(robotDenso,qMatrix(i,:),faces,vertex,faceNormals,false)
    disp('Collision detected!!');
    newRoute = robotDenso.fkine(qMatrix)+ goUp;
    newRoute(:,:,1) = robotDenso.fkine(qMatrix)- goUp;
    newRoute(:,:,100) = robotDenso.fkine(qMatrix)- goUp;

    qMatrix = robotDenso.ikcon(newRoute);
   
    else
    words=['No collision found in the trajectory!'];
    disp(words); 
    break;
    end
    
end

%  Plot the results
for i=1:steps-1
    delete(paperTri([3]));
    posEE=robotDenso.fkine(qMatrix(i,:)); %gets the end-effector in terms of 4x4 matrix
    xOffset = posEE(1,4);   %x of end-effector
    yOffset = posEE(2,4);   %y of end-effector
    zOffset = posEE(3,4)-0.175; %z of end-effector and positioning paper to be under the gripper
    paperTri([3])= trisurf(f,v(:,1) + xOffset,v(:,2) + yOffset, v(:,3) + zOffset ... 
        ,'FaceVertexCData',vertexColours,'EdgeColor','interp','EdgeLighting','flat'); %Replot paperTri
    robotDenso.animate(qMatrix(i,:)); %replot the robot
end
pause(1);
delete(paperTri([3])); 
zOffset = [3.16]; %setting the z to drop the paper off
paperTri([3])= trisurf(f,v(:,1) + xOffset,v(:,2) + yOffset, v(:,3) + zOffset ... 
        ,'FaceVertexCData',vertexColours,'EdgeColor','interp','EdgeLighting','flat'); %Replot paperTri
pause(1);
%% Blue to Mail (second time)
% 1.1) Set parameters for the simulation
t = 10;             % Total time (s)
deltaT = 0.1;      % Control frequency
steps = t/deltaT;   % No. of steps for simulation
delta = 2*pi/steps; % Small angle change
epsilon = 0.1;      % Threshold value for manipulability/Damped Least Squares
W = diag([1 1 1 0.1 0.1 0.1]);    % Weighting matrix for the velocity vector

% 1.2) Allocate array data
m = zeros(steps,1);             % Array for Measure of Manipulability
qMatrix = zeros(steps,6);       % Array for joint angles
qdot = zeros(steps,6);          % Array for joint velocities
theta = zeros(3,steps);         % Array for roll-pitch-yaw angles
x = zeros(3,steps);             % Array for x-y-z trajectory
positionError = zeros(3,steps); % For plotting trajectory error
angleError = zeros(3,steps);    % For plotting trajectory error

% 1.3) Set up trajectory, initial pose
s = lspb(0,1,steps);                % Trapezoidal trajectory scalar

for i=1:steps
    x(1,i) = (1-s(i))*-0.9 + s(i)*1.6; % Points in x
    x(2,i) = (1-s(i))*1.16 +s(i)*0; % Points in y
    x(3,i) = (1-s(i))*4.2 +0.5*sin(i*delta/2)+ s(i)*4;     % Points in z
    theta(1,i) = (1-s(i))*-176*pi/180 + s(i)*-176*pi/180; % Roll angle 
    theta(2,i) = (1-s(i))*0.29*pi/180 + s(i)*1*pi/180;  % Pitch angle
    theta(3,i) = (1-s(i))*80.48*pi/180 + s(i)*-0.01*pi/18; %yaw angle
end
T = [rpy2r(theta(1,1),theta(2,1),theta(3,1)) x(:,1);zeros(1,3) 1];          % Create transformation of first point and angle
q0 = zeros(1,6);                                                            % Initial guess for joint angles
qMatrix(1,:) = robotDenso.ikcon(T,q0);                                            % Solve joint angles to achieve first waypoint

% 1.4) Track the trajectory with RMRC
for i = 1:steps-1
    T = robotDenso.fkine(qMatrix(i,:));                                           % Get forward transformation at current joint state
    deltaX = x(:,i+1) - T(1:3,4);                                         	% Get position error from next waypoint
    Rd = rpy2r(theta(1,i+1),theta(2,i+1),theta(3,i+1));                     % Get next RPY angles, convert to rotation matrix
    Ra = T(1:3,1:3);                                                        % Current end-effector rotation matrix
    Rdot = (1/deltaT)*(Rd - Ra);                                                % Calculate rotation matrix error
    S = Rdot*Ra';                                                           % Skew symmetric!
    linear_velocity = (1/deltaT)*deltaX;
    angular_velocity = [S(3,2);S(1,3);S(2,1)];                              % Check the structure of Skew Symmetric matrix!!
    deltaTheta = tr2rpy(Rd*Ra');                                            % Convert rotation matrix to RPY angles
    xdot = W*[linear_velocity;angular_velocity];                          	% Calculate end-effector velocity to reach next waypoint.
    J = robotDenso.jacob0(qMatrix(i,:));                 % Get Jacobian at current joint state
    m(i) = sqrt(round(det(J*J')));
    if m(i) < epsilon  % If manipulability is less than given threshold
        lambda = (1 - m(i)/epsilon)*5E-2;
    else
        lambda = 0;
    end
    invJ = inv(J'*J + lambda *eye(6))*J';                                   % DLS Inverse
    qdot(i,:) = (invJ*xdot)';                                                % Solve the RMRC equation (you may need to transpose the         vector)
    for j = 1:6                                                             % Loop through joints 1 to 6
        if qMatrix(i,j) + deltaT*qdot(i,j) < robotDenso.qlim(j,1)                     % If next joint angle is lower than joint limit...
            qdot(i,j) = 0; % Stop the motor
        elseif qMatrix(i,j) + deltaT*qdot(i,j) > robotDenso.qlim(j,2)                 % If next joint angle is greater than joint limit ...
            qdot(i,j) = 0; % Stop the motor
        end
    end
    qMatrix(i+1,:) = qMatrix(i,:) + deltaT*qdot(i,:);                         	% Update next joint state based on joint velocities
    positionError(:,i) = x(:,i+1) - T(1:3,4);                               % For plotting
    angleError(:,i) = deltaTheta;                                           % For plotting
end

% Collsion Checking and Re-mapping if needed
while(1)
    faces = fWall;
    vertex = vWall;
    faceNormals = faceNormalsWall;
    goUp = [0,0,0,0;0,0,0,0;0,0,0,0.1;0,0,0,0];

    if IsCollision(robotDenso,qMatrix(i,:),faces,vertex,faceNormals,false)
    disp('Collision detected!!');
    newRoute = robotDenso.fkine(qMatrix)+ goUp;
    newRoute(:,:,1) = robotDenso.fkine(qMatrix)- goUp;
    newRoute(:,:,100) = robotDenso.fkine(qMatrix)- goUp;

    qMatrix = robotDenso.ikcon(newRoute);
   
    else
    words=['No collision found in the trajectory!'];
    disp(words); 
    break;
    end
    
end
% Plot the results
for i=1:steps-1
    robotDenso.animate(qMatrix(i,:)); %(i,:) need for loop
end
pause(1);
%% Mail to Pink Box (second time)
%setting the ply file (for some reason it reads the last ply file read)
[f,v,data] = plyread('paperCircle.ply','tri');
vertexColours = [data.vertex.red, data.vertex.green, data.vertex.blue] / 255;
% 1.1) Set parameters for the simulation
t = 10;             % Total time (s)
deltaT = 0.1;      % Control frequency
steps = t/deltaT;   % No. of steps for simulation
delta = 2*pi/steps; % Small angle change
epsilon = 0.1;      % Threshold value for manipulability/Damped Least Squares
W = diag([1 1 1 0.1 0.1 0.1]);    % Weighting matrix for the velocity vector

% 1.2) Allocate array data
m = zeros(steps,1);             % Array for Measure of Manipulability
qMatrix = zeros(steps,6);       % Array for joint angles
qdot = zeros(steps,6);          % Array for joint velocities
theta = zeros(3,steps);         % Array for roll-pitch-yaw angles
x = zeros(3,steps);             % Array for x-y-z trajectory
positionError = zeros(3,steps); % For plotting trajectory error
angleError = zeros(3,steps);    % For plotting trajectory error

% 1.3) Set up trajectory, initial pose
s = lspb(0,1,steps);                % Trapezoidal trajectory scalar

for i=1:steps
    x(1,i) = (1-s(i))*1.6 + s(i)*-0.9; % Points in x
    x(2,i) = (1-s(i))*0.01 -s(i)*1.1; % Points in y
    x(3,i) = (1-s(i))*4 + s(i)*4.2;     % Points in z
    theta(1,i) = (1-s(i))*-176*pi/180 + s(i)*-176*pi/180;                 % Roll angle 
    theta(2,i) = (1-s(i))*1*pi/180 + s(i)*0.29*pi/180;            % Pitch angle
    theta(3,i) = (1-s(i))*-0.01*pi/18 + s(i)*80.48*pi/180; %yaw angle
end
T = [rpy2r(theta(1,1),theta(2,1),theta(3,1)) x(:,1);zeros(1,3) 1];          % Create transformation of first point and angle
q0 = zeros(1,6);                                                            % Initial guess for joint angles
qMatrix(1,:) = robotDenso.ikcon(T,q0);                                            % Solve joint angles to achieve first waypoint

% 1.4) Track the trajectory with RMRC
for i = 1:steps-1
    T = robotDenso.fkine(qMatrix(i,:));                                           % Get forward transformation at current joint state
    deltaX = x(:,i+1) - T(1:3,4);                                         	% Get position error from next waypoint
    Rd = rpy2r(theta(1,i+1),theta(2,i+1),theta(3,i+1));                     % Get next RPY angles, convert to rotation matrix
    Ra = T(1:3,1:3);                                                        % Current end-effector rotation matrix
    Rdot = (1/deltaT)*(Rd - Ra);                                                % Calculate rotation matrix error
    S = Rdot*Ra';                                                           % Skew symmetric!
    linear_velocity = (1/deltaT)*deltaX;
    angular_velocity = [S(3,2);S(1,3);S(2,1)];                              % Check the structure of Skew Symmetric matrix!!
    deltaTheta = tr2rpy(Rd*Ra');                                            % Convert rotation matrix to RPY angles
    xdot = W*[linear_velocity;angular_velocity];                          	% Calculate end-effector velocity to reach next waypoint.
    J = robotDenso.jacob0(qMatrix(i,:));                 % Get Jacobian at current joint state
    m(i) = sqrt(round(det(J*J')));
    if m(i) < epsilon  % If manipulability is less than given threshold
        lambda = (1 - m(i)/epsilon)*5E-2;
    else
        lambda = 0;
    end
    invJ = inv(J'*J + lambda *eye(6))*J';                                   % DLS Inverse
    qdot(i,:) = (invJ*xdot)';                                                % Solve the RMRC equation (you may need to transpose the         vector)
    for j = 1:6                                                             % Loop through joints 1 to 6
        if qMatrix(i,j) + deltaT*qdot(i,j) < robotDenso.qlim(j,1)                     % If next joint angle is lower than joint limit...
            qdot(i,j) = 0; % Stop the motor
        elseif qMatrix(i,j) + deltaT*qdot(i,j) > robotDenso.qlim(j,2)                 % If next joint angle is greater than joint limit ...
            qdot(i,j) = 0; % Stop the motor
        end
    end
    qMatrix(i+1,:) = qMatrix(i,:) + deltaT*qdot(i,:);                         	% Update next joint state based on joint velocities
    positionError(:,i) = x(:,i+1) - T(1:3,4);                               % For plotting
    angleError(:,i) = deltaTheta;                                           % For plotting
end
% Collsion Checking and Re-mapping if needed
while(1)
    faces = fWall;
    vertex = vWall;
    faceNormals = faceNormalsWall;
    goUp = [0,0,0,0;0,0,0,0;0,0,0,0.1;0,0,0,0];

    if IsCollision(robotDenso,qMatrix(i,:),faces,vertex,faceNormals,false)
    disp('Collision detected!!');
    newRoute = robotDenso.fkine(qMatrix)+ goUp;
    newRoute(:,:,1) = robotDenso.fkine(qMatrix)- goUp;
    newRoute(:,:,100) = robotDenso.fkine(qMatrix)- goUp;

    qMatrix = robotDenso.ikcon(newRoute);
   
    else
    words=['No collision found in the trajectory!'];
    disp(words); 
    break;
    end
    
end

% Plot the results
for i=1:steps-1
    delete(paperCircle([1]));
    posEE=robotDenso.fkine(qMatrix(i,:)); %gets the end-effector in terms of 4x4 matrix
    xOffset = posEE(1,4);   %x of end-effector
    yOffset = posEE(2,4);   %y of end-effector
    zOffset = posEE(3,4)-0.175; %z of end-effector and positioning paper to be under the gripper
    paperCircle([1])= trisurf(f,v(:,1) + xOffset,v(:,2) + yOffset, v(:,3) + zOffset ... 
        ,'FaceVertexCData',vertexColours,'EdgeColor','interp','EdgeLighting','flat'); %replot the paperCircle
    robotDenso.animate(qMatrix(i,:)); %replot the robot
end
pause(1);
delete(paperCircle([1]));
zOffset = [3.16]; %setting the z to drop the paper off
paperCircle([1])= trisurf(f,v(:,1) + xOffset,v(:,2) + yOffset, v(:,3) + zOffset ... 
        ,'FaceVertexCData',vertexColours,'EdgeColor','interp','EdgeLighting','flat'); %replot the paperCircle
pause(1);
%% Pink Box to Mail (second time)
% 1.1) Set parameters for the simulation
t = 10;             % Total time (s)
deltaT = 0.1;      % Control frequency
steps = t/deltaT;   % No. of steps for simulation
delta = 2*pi/steps; % Small angle change
epsilon = 0.1;      % Threshold value for manipulability/Damped Least Squares
W = diag([1 1 1 0.1 0.1 0.1]);    % Weighting matrix for the velocity vector

% 1.2) Allocate array data
m = zeros(steps,1);             % Array for Measure of Manipulability
qMatrix = zeros(steps,6);       % Array for joint angles
qdot = zeros(steps,6);          % Array for joint velocities
theta = zeros(3,steps);         % Array for roll-pitch-yaw angles
x = zeros(3,steps);             % Array for x-y-z trajectory
positionError = zeros(3,steps); % For plotting trajectory error
angleError = zeros(3,steps);    % For plotting trajectory error

% 1.3) Set up trajectory, initial pose
s = lspb(0,1,steps);                % Trapezoidal trajectory scalar

for i=1:steps
    x(1,i) = (1-s(i))*-0.9 + s(i)*1.6; % Points in x
    x(2,i) = (1-s(i))*-1.16 +s(i)*0; % Points in y
    x(3,i) = (1-s(i))*4.2 + s(i)*4;     % Points in z
    theta(1,i) = (1-s(i))*-176*pi/180 + s(i)*-176*pi/180; % Roll angle 
    theta(2,i) = (1-s(i))*0.29*pi/180 + s(i)*1*pi/180;  % Pitch angle
    theta(3,i) = (1-s(i))*80.48*pi/180 + s(i)*-0.01*pi/18; %yaw angle
end
T = [rpy2r(theta(1,1),theta(2,1),theta(3,1)) x(:,1);zeros(1,3) 1];          % Create transformation of first point and angle
q0 = zeros(1,6);                                                            % Initial guess for joint angles
qMatrix(1,:) = robotDenso.ikcon(T,q0);                                            % Solve joint angles to achieve first waypoint

% 1.4) Track the trajectory with RMRC
for i = 1:steps-1
    T = robotDenso.fkine(qMatrix(i,:));                                           % Get forward transformation at current joint state
    deltaX = x(:,i+1) - T(1:3,4);                                         	% Get position error from next waypoint
    Rd = rpy2r(theta(1,i+1),theta(2,i+1),theta(3,i+1));                     % Get next RPY angles, convert to rotation matrix
    Ra = T(1:3,1:3);                                                        % Current end-effector rotation matrix
    Rdot = (1/deltaT)*(Rd - Ra);                                                % Calculate rotation matrix error
    S = Rdot*Ra';                                                           % Skew symmetric!
    linear_velocity = (1/deltaT)*deltaX;
    angular_velocity = [S(3,2);S(1,3);S(2,1)];                              % Check the structure of Skew Symmetric matrix!!
    deltaTheta = tr2rpy(Rd*Ra');                                            % Convert rotation matrix to RPY angles
    xdot = W*[linear_velocity;angular_velocity];                          	% Calculate end-effector velocity to reach next waypoint.
    J = robotDenso.jacob0(qMatrix(i,:));                 % Get Jacobian at current joint state
    m(i) = sqrt(round(det(J*J')));
    if m(i) < epsilon  % If manipulability is less than given threshold
        lambda = (1 - m(i)/epsilon)*5E-2;
    else
        lambda = 0;
    end
    invJ = inv(J'*J + lambda *eye(6))*J';                                   % DLS Inverse
    qdot(i,:) = (invJ*xdot)';                                                % Solve the RMRC equation (you may need to transpose the         vector)
    for j = 1:6                                                             % Loop through joints 1 to 6
        if qMatrix(i,j) + deltaT*qdot(i,j) < robotDenso.qlim(j,1)                     % If next joint angle is lower than joint limit...
            qdot(i,j) = 0; % Stop the motor
        elseif qMatrix(i,j) + deltaT*qdot(i,j) > robotDenso.qlim(j,2)                 % If next joint angle is greater than joint limit ...
            qdot(i,j) = 0; % Stop the motor
        end
    end
    qMatrix(i+1,:) = qMatrix(i,:) + deltaT*qdot(i,:);                         	% Update next joint state based on joint velocities
    positionError(:,i) = x(:,i+1) - T(1:3,4);                               % For plotting
    angleError(:,i) = deltaTheta;                                           % For plotting
end
% Collsion Checking and Re-mapping if needed
while(1)
    faces = fWall;
    vertex = vWall;
    faceNormals = faceNormalsWall;
    goUp = [0,0,0,0;0,0,0,0;0,0,0,0.1;0,0,0,0];

    if IsCollision(robotDenso,qMatrix(i,:),faces,vertex,faceNormals,false)
    disp('Collision detected!!');
    newRoute = robotDenso.fkine(qMatrix)+ goUp;
    newRoute(:,:,1) = robotDenso.fkine(qMatrix)- goUp;
    newRoute(:,:,100) = robotDenso.fkine(qMatrix)- goUp;

    qMatrix = robotDenso.ikcon(newRoute);
   
    else
    words=['No collision found in the trajectory!'];
    disp(words); 
    break;
    end
    
end

% Plot the results
for i=1:steps-1
    robotDenso.animate(qMatrix(i,:)); %(i,:) need for loop
end
pause(1);





%% IsCollision
% This is based upon the output of questions 2.5 and 2.6
% Given a robot model (robot), and trajectory (i.e. joint state vector) (qMatrix)
% and triangle obstacles in the environment (faces,vertex,faceNormals)
function result = IsCollision(robot,qMatrix,faces,vertex,faceNormals,returnOnceFound)
if nargin < 6
    returnOnceFound = true;
end
result = false;

for qIndex = 1:size(qMatrix,1)
    % Get the transform of every joint (i.e. start and end of every link)
     tr = GetLinkPoses(qMatrix(qIndex,:), robot);

    % Go through each link and also each triangle face
    for i = 1 : size(tr,3)-1    
        for faceIndex = 1:size(faces,1)
            vertOnPlane = vertex(faces(faceIndex,1)',:);
            [intersectP,check] = LinePlaneIntersection(faceNormals(faceIndex,:),vertOnPlane,tr(1:3,4,i)',tr(1:3,4,i+1)'); 
            if check == 1 && IsIntersectionPointInsideTriangle(intersectP,vertex(faces(faceIndex,:)',:))
                plot3(intersectP(1),intersectP(2),intersectP(3),'g*');
                display('Intersection');
                result = true;
                if returnOnceFound
                    
                    return
                end
            end
        end    
    end
end
end
%% GetLinkPoses
% q - robot joint angles
% robot -  seriallink robot model
% transforms - list of transforms
function [ transforms ] = GetLinkPoses( q, robotDenso)

links = robotDenso.links;
transforms = zeros(4, 4, length(links) + 1);
transforms(:,:,1) = robotDenso.base;

for i = 1:length(links)
    L = links(1,i);
    
    current_transform = transforms(:,:, i);
    
    current_transform = current_transform * trotz(q(1,i) + L.offset) * ...
    transl(0,0, L.d) * transl(L.a,0,0) * trotx(L.alpha);
    transforms(:,:,i + 1) = current_transform;
end
end
%% IsIntersectionPointInsideTriangle
% Given a point which is known to be on the same plane as the triangle
% determine if the point is 
% inside (result == 1) or 
% outside a triangle (result ==0 )
function result = IsIntersectionPointInsideTriangle(intersectP,triangleVerts)

u = triangleVerts(2,:) - triangleVerts(1,:);
v = triangleVerts(3,:) - triangleVerts(1,:);

uu = dot(u,u);
uv = dot(u,v);
vv = dot(v,v);

w = intersectP - triangleVerts(1,:);
wu = dot(w,u);
wv = dot(w,v);

D = uv * uv - uu * vv;

% Get and test parametric coords (s and t)
s = (uv * wv - vv * wu) / D;
if (s < 0.0 || s > 1.0)        % intersectP is outside Triangle
    result = 0;
    return;
end

t = (uv * wu - uu * wv) / D;
if (t < 0.0 || (s + t) > 1.0)  % intersectP is outside Triangle
    result = 0;
    return;
end

result = 1;                      % intersectP is in Triangle
end