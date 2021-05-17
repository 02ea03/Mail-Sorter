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