%lab : Arif, Md Arifuzzaman
%2D kinematics lab

clear
clc
close all


P = imread('calibrationframe.jpeg');
image(P)
calp = ginput(4)

cal1=[0,0];cal2=[30,0];cal3=[30,30];cal4=[0,30];


%calp(1,1)*c1+calp(1,2)*c2+c3+0*c4+0*c6-calp(1,1)*cal1(1,1)*c7-calp(1,2)*cal1(1,1)*c8 = cal1(1,1)

A1 = [calp(1,1) calp(1,2) 1 0 0 0 0 0];
A2 = [0 0 0 calp(1,1) calp(1,2) 1 0 0];
A3 = [calp(2,1) calp(2,2) 1 0 0 0 -calp(2,1)*30 -calp(2,2)*30];
A4 = [0 0 0 calp(2,1) calp(2,2) 1 0 0];
A5 = [calp(3,1) calp(3,2) 1 0 0 0 -calp(3,1)*30 -calp(3,2)*30];
A6 = [0 0 0 calp(3,1) calp(3,2) 1 -calp(3,1)*30 -calp(3,2)*30];
A7 = [calp(4,1) calp(4,2) 1 0 0 0 0 0];
A8 = [0 0 0 calp(4,1) calp(4,2) 1 -calp(4,1)*30 -calp(4,2)*30];

A = [A1;A2;A3;A4;A5;A6;A7;A8];
B = [0;0;30;0;30;30;0;30];

C = linsolve(A,B);


    

