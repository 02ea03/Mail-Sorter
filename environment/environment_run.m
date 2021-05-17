%% Environment building 
% clf;clc;clear all;close all;
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
disp('Finished Loading Background')

%% option 1
%setting up the light curtain
n=10;
posLight = zeros(n,9);
for i = 3:6
    x1(i,:) = 2.4;
    x2(i,:) = 2.4;
    x3(i,:) = -2.4;
    x4(i,:) = -2.4;    
    y1(i,:) = 3.6;
    y2(i,:) = -3.6;
    y3(i,:) = -3.6;
    y4(i,:) = 3.6;    
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
view(-15,21);