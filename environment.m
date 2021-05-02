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
%% Robot
robot = VP6242(false);
q = zeros(1,6);
robot.model.base = transl(0,0,0);
robot.model.animate(q);
hold on;
