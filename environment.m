clf;clc;clear all;close all;
%doCameraSpin = false;

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
        for yOffset = [-2]
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
        for yOffset = [2]
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
    %% red frence
    [f,v,data] = plyread('redfence.ply','tri');
    vertexColours = [data.vertex.red, data.vertex.green, data.vertex.blue] / 255;

    for xOffset = [-6.5]
        for yOffset = [3]
            trisurf(f,v(:,1)+ xOffset,v(:,2) + yOffset, v(:,3) ...
                ,'FaceVertexCData',vertexColours,'EdgeColor','interp','EdgeLighting','flat');
        end
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
    
    
        %% robot human
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
        for xOffset = [1.4]            
            trisurf(f,v(:,1) + xOffset,v(:,2) + yOffset, v(:,3) + zOffset ...
        ,'FaceVertexCData',vertexColours,'EdgeColor','interp','EdgeLighting','flat');
        end
    end
end

[f,v,data] = plyread('paperTri.ply','tri');
vertexColours = [data.vertex.red, data.vertex.green, data.vertex.blue] / 255;
for zOffset = [3.51]
    for yOffset = [0]
        for xOffset = [1.4]          
            paperTri([1]) = trisurf(f,v(:,1) + xOffset,v(:,2) + yOffset, v(:,3) + zOffset ... %delete(paperTri([1]) to delete the
        ,'FaceVertexCData',vertexColours,'EdgeColor','interp','EdgeLighting','flat');
        end
    end
end
%% Robot
robot = VP6242(false);
q = zeros(1,6);
robot.model.base = transl(0,0,3.2);
view(-15,21);
robot.model.animate(q);
hold on;
view(-15,21);
%% Standing Position to Mail
%function PosInit2Mail()
% 1.1) Set parameters for the simulation
t = 10;             % Total time (s)
deltaT = 0.1;      % Control frequency
steps = t/deltaT;   % No. of steps for simulation
delta = 2*pi/steps; % Small angle change
epsilon = 0.1;      % Threshold value for manipulability/Damped Least Squares
W = diag([1 1 1 0.1 0.1 0.1]);    % Weighting matrix for the velocity vector

%POS1 = transl(2,0,5);
%POS2 = transl(2,0,3.5);
% 1.2) Allocate array data
m = zeros(steps,1);             % Array for Measure of Manipulability
qMatrix = zeros(steps,6);       % Array for joint anglesR
qdot = zeros(steps,6);          % Array for joint velocities
theta = zeros(3,steps);         % Array for roll-pitch-yaw angles
x = zeros(3,steps);             % Array for x-y-z trajectory
positionError = zeros(3,steps); % For plotting trajectory error
angleError = zeros(3,steps);    % For plotting trajectory error

% 1.3) Set up trajectory, initial pose
s = lspb(0,1,steps);                % Trapezoidal trajectory scalar

for i = 1 : size(steps,3)-1            % Collision (steps was tr)
    for faceIndex = 1:size(faces,1)
        vertOnPlane = vertex(faces(faceIndex,1)',:);
        [intersectP,check] = LinePlaneIntersection(faceNormals(faceIndex,:),vertOnPlane,steps(1:3,4,i)',steps(1:3,4,i+1)'); %(steps was tr)
        if check == 1 && IsIntersectionPointInsideTriangle(intersectP,vertex(faces(faceIndex,:)',:))
            plot3(intersectP(1),intersectP(2),intersectP(3),'g*');
            display('Intersection');
        end
    end    
end

for i=1:steps
    x(1,i) = (1-s(i))*2 + s(i)*2; % Points in x
    x(2,i) = (1-s(i))*0 + s(i)*0; % Points in y
    x(3,i) = (1-s(i))*5 + s(i)*3.6; % Points in z
    theta(1,i) = 0;                 % Roll angle 
    theta(2,i) = 5*pi/9;            % Pitch angle
    theta(3,i) = 0;                 % Yaw angle
    %x(i,:) = (1-s(i))*POS1 + s(i)*POS2;
end
 
T = [rpy2r(theta(1,1),theta(2,1),theta(3,1)) x(:,1);zeros(1,3) 1];          % Create transformation of first point and angle
q0 = zeros(1,6);                                                            % Initial guess for joint angles
qMatrix(1,:) = robot.model.ikcon(T,q0);                                            % Solve joint angles to achieve first waypoint

% 1.4) Track the trajectory with RMRC
for i = 1:steps-1
    T = robot.model.fkine(qMatrix(i,:));                                           % Get forward transformation at current joint state
    deltaX = x(:,i+1) - T(1:3,4);                                         	% Get position error from next waypoint
    Rd = rpy2r(theta(1,i+1),theta(2,i+1),theta(3,i+1));                     % Get next RPY angles, convert to rotation matrix
    Ra = T(1:3,1:3);                                                        % Current end-effector rotation matrix
    Rdot = (1/deltaT)*(Rd - Ra);                                                % Calculate rotation matrix error
    S = Rdot*Ra';                                                           % Skew symmetric!
    linear_velocity = (1/deltaT)*deltaX;
    angular_velocity = [S(3,2);S(1,3);S(2,1)];                              % Check the structure of Skew Symmetric matrix!!
    deltaTheta = tr2rpy(Rd*Ra');                                            % Convert rotation matrix to RPY angles
    xdot = W*[linear_velocity;angular_velocity];                          	% Calculate end-effector velocity to reach next waypoint.
    J = robot.model.jacob0(qMatrix(i,:));                 % Get Jacobian at current joint state
    m(i) = sqrt(round(det(J*J')));
    if m(i) < epsilon  % If manipulability is less than given threshold
        lambda = (1 - m(i)/epsilon)*5E-2;
    else
        lambda = 0;
    end
    invJ = inv(J'*J + lambda *eye(6))*J';                                   % DLS Inverse
    qdot(i,:) = (invJ*xdot)';                                                % Solve the RMRC equation (you may need to transpose the         vector)
    for j = 1:6                                                             % Loop through joints 1 to 6
        if qMatrix(i,j) + deltaT*qdot(i,j) < robot.model.qlim(j,1)                     % If next joint angle is lower than joint limit...
            qdot(i,j) = 0; % Stop the motor
        elseif qMatrix(i,j) + deltaT*qdot(i,j) > robot.model.qlim(j,2)                 % If next joint angle is greater than joint limit ...
            qdot(i,j) = 0; % Stop the motor
        end
    end
    qMatrix(i+1,:) = qMatrix(i,:) + deltaT*qdot(i,:);                         	% Update next joint state based on joint velocities
    positionError(:,i) = x(:,i+1) - T(1:3,4);                               % For plotting
    angleError(:,i) = deltaTheta;                                           % For plotting
end
% 1.5) Plot the results
for 0:1:steps
robot.model.animate(qMatrix(i,:))
paperCircle %look video in resources
end
view(-15,21);
%end

%% Mail to Right Box
%how to pick up the mail and move it in motion

%% Right Box to Mail


%% Mail to Left Box


%% Left Box to Mail