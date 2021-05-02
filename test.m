clc;clear all;close all;
camlight;
axis equal;
view(3);
hold on;
robot = VP6242(false);
q = zeros(1,6);
robot.model.base = transl(0,0,0);
robot.model.animate(q);


% %% parts check
%     [f,v,data] = plyread('VP6242Link0.ply','tri');
%     vertexColours = [data.vertex.red, data.vertex.green, data.vertex.blue] / 255;
% 
%     for xOffset = [0]
%         for yOffset = [0]
%             trisurf(f,v(:,1)+ xOffset,v(:,2) + yOffset, v(:,3) ...
%                 ,'FaceVertexCData',vertexColours,'EdgeColor','interp','EdgeLighting','flat');
%         end
%     end
%         [f,v,data] = plyread('VP6242Link1.ply','tri');
%     vertexColours = [data.vertex.red, data.vertex.green, data.vertex.blue] / 255;
% 
%     for xOffset = [0]
%         for yOffset = [0]
%             trisurf(f,v(:,1)+ xOffset,v(:,2) + yOffset, v(:,3) ...
%                 ,'FaceVertexCData',vertexColours,'EdgeColor','interp','EdgeLighting','flat');
%         end
%     end
%         [f,v,data] = plyread('VP6242Link2.ply','tri');
%     vertexColours = [data.vertex.red, data.vertex.green, data.vertex.blue] / 255;
% 
%     for xOffset = [0]
%         for yOffset = [0]
%             trisurf(f,v(:,1)+ xOffset,v(:,2) + yOffset, v(:,3) ...
%                 ,'FaceVertexCData',vertexColours,'EdgeColor','interp','EdgeLighting','flat');
%         end
%     end
%         [f,v,data] = plyread('VP6242Link3.ply','tri');
%     vertexColours = [data.vertex.red, data.vertex.green, data.vertex.blue] / 255;
% 
%     for xOffset = [0]
%         for yOffset = [0]
%             trisurf(f,v(:,1)+ xOffset,v(:,2) + yOffset, v(:,3) ...
%                 ,'FaceVertexCData',vertexColours,'EdgeColor','interp','EdgeLighting','flat');
%         end
%     end
%         [f,v,data] = plyread('VP6242Link4.ply','tri');
%     vertexColours = [data.vertex.red, data.vertex.green, data.vertex.blue] / 255;
% 
%     for xOffset = [0]
%         for yOffset = [0]
%             trisurf(f,v(:,1)+ xOffset,v(:,2) + yOffset, v(:,3) ...
%                 ,'FaceVertexCData',vertexColours,'EdgeColor','interp','EdgeLighting','flat');
%         end
%     end
%         [f,v,data] = plyread('VP6242Link5.ply','tri');
%     vertexColours = [data.vertex.red, data.vertex.green, data.vertex.blue] / 255;
% 
%     for xOffset = [0]
%         for yOffset = [0]
%             trisurf(f,v(:,1)+ xOffset,v(:,2) + yOffset, v(:,3) ...
%                 ,'FaceVertexCData',vertexColours,'EdgeColor','interp','EdgeLighting','flat');
%         end
%     end
