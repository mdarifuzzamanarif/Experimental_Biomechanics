%lab1:2D Kinematics
%Arif, Md. Arifuzzaman
%part 3
clear
clc
clear all

%load data
load('vicondata.mat');
load('pank2.mat');load('pkne2.mat');load('phip2.mat');

%Vicon data anlysis
pkne_v = [vicondata(:,1) vicondata(:,2)];
pank_v = [vicondata(:,3) vicondata(:,4)];
phip_v = [vicondata(:,5) vicondata(:,6)];

for i = 1:length(pank_v)
   
    t_v = pkne_v(i,:) - phip_v(i,:);
    l_v = pank_v(i,:) - pkne_v(i,:);
    
    ut_v=t_v/norm(t_v);
    ul_v=l_v/norm(l_v);
    
    mag_utv = norm(ut_v);
    mag_ulv = norm(ul_v);
    
    theta_v(i) = acosd(dot(ut_v,ul_v)/(mag_utv*mag_ulv));
end
theta_v = theta_v';
thetaVmin = min(theta_v);
thetaVmax = max(theta_v);
ROM_v = (thetaVmax -thetaVmin);

% figure(1)
% plot(theta_v)
% axis equal
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for p=1:55
pank_tx(p,1) = (pank(p,1)) + (mean(vicondata(:,3))- mean(pank(:,1)));
pank_ty(p,1) = (pank(p,2)) + (mean(vicondata(:,4))- mean(pank(:,2)));

pkne_tx(p,1) = (pkne(p,1)) + (mean(vicondata(:,1))- mean(pkne(:,1)));
pkne_ty(p,1) = (pkne(p,2)) + (mean(vicondata(:,2))- mean(pkne(:,2)));

phip_tx(p,1) = (phip(p,1)) + (mean(vicondata(:,5))- mean(phip(:,1)));
phip_ty(p,1) = (phip(p,2)) + (mean(vicondata(:,6))- mean(phip(:,2)));

end

pank_m = [pank_tx pank_ty];
pkne_m = [pkne_tx pkne_ty];
phip_m = [phip_tx phip_ty];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%manual data analysis
for i = 1:length(pank)
   
%     t_m = -(pkne_m(i,:) - phip_m(i,:));
%     l_m = pank_m(i,:) - pkne_m(i,:);
    
    t_m = -(pkne(i,:) - phip(i,:));
    l_m = pank(i,:) - pkne(i,:);
    
    ut_m=t_m/norm(t_m);
    ul_m=l_m/norm(l_m);
    
    mag_utm=norm(ut_m);
    mag_ulm= norm(ul_m);
    
    theta_m(i)=acosd(dot(ut_m,ul_m)/(mag_utm*mag_ulm));
end

theta_m = theta_m';

thetaMmin = min(theta_m);
thetaMmax = max(theta_m);
ROM_m = (thetaMmax -thetaMmin);

fs1= 50;
t1 = (1/fs1)*(1:length(pank_v));

fs2= 25;
t2 = (1/fs2)*(1:length(pank));

figure(2)
plot(t2,theta_m,'g')
hold on
plot(t1,theta_v,'b')
xlabel('time')
ylabel('flexion angle')
hold off
legend('Manual Data','Vicon Data')
% xlim([0 5])
% ylim([100 200])

figure(3)
subplot(3,2,1)

plot(t1,pank_v(:,1),t2,pank(:,1))
xlabel('Time (s)')
ylabel('Ankle position- horizontal(mm)')
legend('Vicon Data','Manual Data')

subplot(3,2,2)
plot(t1,pank_v(:,2),t2,pank(:,2))
xlabel('Time (s)')
ylabel('Ankle position- vertical (mm)')
legend('Vicon Data','Manual Data')

subplot(3,2,3)
plot(t1,pkne_v(:,1),t2,pkne(:,1))
xlabel('Time (s)')
ylabel('Knee position -horizontal(mm)')
legend('Vicon Data','Manual Data')

subplot(3,2,4)
plot(t1,pkne_v(:,2),t2,pkne(:,2))
xlabel('Time (s)')
ylabel('Knee position -vertical (mm)')
legend('Vicon Data','Manual Data')

subplot(3,2,5)
plot(t1,phip_v(:,1),t2,phip(:,1))
xlabel('Time (s)')
ylabel('Hip position-horizontal(mm)')
legend('Vicon Data','Manual Data')

subplot(3,2,6)
plot(t1,phip_v(:,2),t2,phip(:,2))
xlabel('Time (s)')
ylabel('Hip position - vertical (mm)')
legend('Vicon Data','Manual Data')

%RMS error calculation
%This has been done by using a m-file rms.m from mathworks, 
%m-file rmse.m didn't work


% RMSE_ankx = rmse(pank_v(:,1),pank(:,1));
% RMSE_anky= rmse(pank_v(:,2),pank(:,2));
% RMSE_knex = rmse(pkne_v(:,1),pkne(:,1));
% RMSE_kney = rmse(pkne_v(:,2),pkne(:,2));
% RMSE_hipx = rmse(phip_v(:,1),phip(:,1));
% RMSE_hipy = rmse(phip_v(:,2),phip(:,2));
% RMS_theta = rmse(theta_v,theta_m);

ankx_rmse = rms(pank_v(:,1)-rms(pank_m(:,1)))
anky_rmse = rms(pank_v(:,2)-rms(pank_m(:,2)))
knex_rmse = rms(pkne_v(:,1)-rms(pkne_m(:,1)))
kney_rmse = rms(pkne_v(:,2)-rms(pkne_m(:,2)))
hipx_rmse = rms(phip_v(:,1)-rms(phip_m(:,1)))
hipy_rmse = rms(phip_v(:,2)-rms(phip_m(:,2)))
theta_rmse = rms(theta_v(:,1)-rms(theta_m(:,1)))

figure(5)
subplot(3,2,1)

plot(t1,pank_v(:,1),t2,pank_m(:,1))
xlabel('Time (s)')
ylabel('Ankle position- horizontal(mm)')
legend('Vicon Data','Manual Data')

subplot(3,2,2)
plot(t1,pank_v(:,2),t2,pank_m(:,2))
xlabel('Time (s)')
ylabel('Ankle position- vertical (mm)')
legend('Vicon Data','Manual Data')

subplot(3,2,3)
plot(t1,pkne_v(:,1),t2,pkne_m(:,1))
xlabel('Time (s)')
ylabel('Knee position -horizontal(mm)')
legend('Vicon Data','Manual Data')

subplot(3,2,4)
plot(t1,pkne_v(:,2),t2,pkne_m(:,2))
xlabel('Time (s)')
ylabel('Knee position -vertical (mm)')
legend('Vicon Data','Manual Data')

subplot(3,2,5)
plot(t1,phip_v(:,1),t2,phip_m(:,1))
xlabel('Time (s)')
ylabel('Hip position-horizontal(mm)')
legend('Vicon Data','Manual Data')

subplot(3,2,6)
plot(t1,phip_v(:,2),t2,phip_m(:,2))
xlabel('Time (s)')
ylabel('Hip position - vertical (mm)')
legend('Vicon Data','Manual Data')

