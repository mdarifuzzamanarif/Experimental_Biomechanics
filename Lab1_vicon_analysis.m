%lab 1: 2D Kinematics
%Arif, Md Arifuzzaman


%Vicon data anlysis
clear
clc
close all

load('vicondata.mat');

pkneV = [vicondata(:,1) vicondata(:,2)];
pankV = [vicondata(:,3) vicondata(:,4)];
phipV = [vicondata(:,5) vicondata(:,6)];

% figure(5)
% plot(pankV(:,1))


for i = 1:length(pankV)
   
    TV = pkneV(i,:) - phipV(i,:);
    LV = pankV(i,:) - pkneV(i,:);
    
    uTV=TV/norm(TV);
    uLV=LV/norm(LV);
    
    mag_uTV = norm(uTV);
    mag_uLV = norm(uLV);
    
    thetaV(i) = acosd(dot(uTV,uLV)/(mag_uTV*mag_uLV));
end
Vtheta = thetaV';

MinVTheta = min(Vtheta);
MaxVTheta = max(Vtheta);
RangeofMotionV = MaxVTheta - MinVTheta;

% plot(pankV(:,1),pankV(:,2))

plot(Vtheta)
axis equal
