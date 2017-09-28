clear
clc
close all

load('CalCoeff.mat');
c1=C(1,1);c2=C(2,1);c3=C(3,1);c4=C(4,1);c5=C(5,1);c6=C(6,1);c7=C(7,1);c8=C(8,1);


for i=138:192
fstr = int2str(i);
p = imread(['P',fstr,'.jpeg']);
image(p)
uv = ginput(3);

for j=1:3
x(j)=(c1*uv(j,1)+c2*uv(j,2)+c3)/(1+c7*uv(j,1)+c8*uv(j,2));
y(j)=(c4*uv(j,1)+c5*uv(j,2)+c6)/(1+c7*uv(j,1)+c8*uv(j,2));
end

pank(i,:)=[x(1) y(1)];
pkne(i,:)=[x(2) y(2)];
phip(i,:)=[x(3) y(3)];
end

% pank= pank(138:end,:);
% pkne= pkne(138:end,:);
% phip= phip(138:end,:);

pank(1:137,:)=[];
pkne(1:137,:)=[];
phip(1:137,:)=[];
