classdef VP6242 < handle
    properties
        %> Robot model
        model;
        
        %>
        workspace = [-0.5 0.5 -0.5 0.5 -0.5 0.5];   
        
        %> Flag to indicate if gripper is used
        useGripper = false;        
    end
    
    methods%% Class for VP6242 robot simulation
function self = VP6242(useGripper)
    if nargin < 1
        useGripper = false;
    end
    self.useGripper = useGripper;
%> Define the boundaries of the workspace

        
% robot = 
self.GetVP6242Robot();
% robot = 
self.PlotAndColourRobot();%robot,workspace);
end

%% GetVP6242Robot
% Given a name (optional), create and return a VP6242 robot model
function GetVP6242Robot(self)
%     if nargin < 1
        % Create a unique name (ms timestamp after 1ms pause)
        pause(0.001);
        name = ['VP_6242_',datestr(now,'yyyymmddTHHMMSSFFF')];
%     end

    % DnH

    L1 = Link('d',0.125,'a',0,'alpha',0,'offset',0,'qlim',deg2rad([-160 160]));
    L2 = Link('d',0,'a',0.21,'alpha',0,'offset',0,'qlim',deg2rad([-120 120]));
    L3 = Link('d',0,'a',-0.075,'alpha',-pi/2,'offset',0,'qlim',deg2rad([19 160]));
    L4 = Link('d',0.21,'a',0,'alpha',pi/2,'offset',0,'qlim',deg2rad([-160 160]));
    L5 = Link('d',0.08535,'a',0,'alpha',-pi/2,'offset',0,'qlim',deg2rad([-120 120]));
    L6 = Link('d',0.07,'a',0,'alpha',0,'offset',0,'qlim',deg2rad([-360 360]));
    self.model = SerialLink([L1 L2 L3 L4 L5 L6],'name','VP6242');
end
%% PlotAndColourRobot
% Given a robot index, add the glyphs (vertices and faces) and
% colour them in if data is available 
function PlotAndColourRobot(self)%robot,workspace)
    for linkIndex = 1:self.model.n
        if self.useGripper && linkIndex == self.model.n
            [ faceData, vertexData, plyData{linkIndex+1} ] = plyread(['VP6242Link',num2str(linkIndex-1),'Gripper.ply'],'tri'); %#ok<AGROW>
        else
            display(['VP6242Link',num2str(linkIndex-1),'.PLY'])
            [ faceData, vertexData, plyData{linkIndex+1} ] = plyread(['VP6242Link',num2str(linkIndex-1),'.PLY'],'tri'); %#ok<AGROW>
        end
        self.model.faces{linkIndex+1} = faceData;
        self.model.points{linkIndex+1} = vertexData;
    end

    % Display robot
    self.model.plot3d(zeros(1,self.model.n),'noarrow','workspace',self.workspace);
    if isempty(findobj(get(gca,'Children'),'Type','Light'))
        camlight
    end  
    self.model.delay = 0;

    % Try to correctly colour the arm (if colours are in ply file data)
    for linkIndex = 0:self.model.n
        handles = findobj('Tag', self.model.name);
        h = get(handles,'UserData');
        try 
            h.link(linkIndex+1).Children.FaceVertexCData = [plyData{linkIndex+1}.vertex.red ...
                                                          , plyData{linkIndex+1}.vertex.green ...
                                                          , plyData{linkIndex+1}.vertex.blue]/255;
            h.link(linkIndex+1).Children.FaceColor = 'interp';
        catch ME_1
            disp(ME_1);
            continue;
        end
    end
end    
    end
end