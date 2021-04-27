clc;clear all;close all;

robot = VP6242(false);
q = zeros(1,6);
robot.model.base = transl(0,0,0);
robot.model.animate(q);

%%
workspace = [-5 5 -5 5 -1 5]; 
surf([-5,-5;5,5],[-5,5;-5,5],[0.01,0.01;0.01,0.01],'CData',imread('concrete.jpg'),'FaceColor','texturemap');
hold on;