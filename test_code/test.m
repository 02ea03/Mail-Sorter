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


a = transl(0,0,0);
b = transl(1,1,1);
qInitial = zero(1,6);
steps = 100;
s = lspb(0,1,steps);
for i = 1 : size(tr,3)-1    
    for faceIndex = 1:size(faces,1)
        vertOnPlane = vertex(faces(faceIndex,1)',:);
        [intersectP,check] = LinePlaneIntersection(faceNormals(faceIndex,:),vertOnPlane,tr(1:3,4,i)',tr(1:3,4,i+1)'); 
        if check == 1 && IsIntersectionPointInsideTriangle(intersectP,vertex(faces(faceIndex,:)',:))
            plot3(intersectP(1),intersectP(2),intersectP(3),'g*');
            display('Intersection');
        end
    end    
end
q0=robot.model.ikcon(a,qInitial);
q1=robot.model.ikcon(b,qInitial);

for i=1:steps
    x(i,:) = (1-s(i))*a + s(i)*b;
end



