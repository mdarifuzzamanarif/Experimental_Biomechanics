%Arif, Md Arifuzzaman, UMKC ID:16209626
%lab 2
%3D Kinermatics Lab
clear
clc
close all
%part 1: calibration trial
load ('data.mat'); %

% %filtering
% calibration_data=filter_data(calibration_data, 10, 100, 1:36);
% gait_data=filter_data(gait_data, 10, 100, 1:18);

caldata = mean(calibration_data); %
%step 1: joint centres and centre of mass
% Ankle Joint Center
pa_m = caldata(:,31:33); %Medial Ankle
pa_l = caldata(:,34:36); %Lateral Ankle
pa = pa_m + (pa_l-pa_m)/2; 

% Knee Joint Center
pk_m = caldata(:,25:27); %Medial Knee
pk_l = caldata(:,28:30); %Lateral Knee
pk = pk_m + (pk_l-pk_m)/2;

% Hip Joint Center
ph_LA = caldata(:,1:3); %Left ASIS
ph_RA = caldata(:,4:6); %Right ASIS
ph_amid = ph_LA + (ph_RA-ph_LA)/2;
dasis =  pdist2(ph_LA,ph_RA);
phip = [ph_amid(1) + 0.36.*dasis ph_amid(2) - 0.19.*dasis ph_amid(3) - 0.30.*dasis];

%Center of Mass
cm_leg = pa + 0.567*(pk-pa); %leg
cm_thi = pk + 0.567+(phip-pk); %thigh

%Step 2:
itemp_l = pk_l-pk; %for both leg and thigh

%leg
ktemp_l = pk-pa; 
jtemp_l = cross(ktemp_l,itemp_l);
itemp2_l = cross(jtemp_l,ktemp_l);

ilu = itemp2_l./norm(itemp2_l);
jlu = jtemp_l./norm(jtemp_l);
klu = ktemp_l./norm(ktemp_l);

%TA2G_L = [ilu' jlu' klu'];
TG2A_L = [ilu;jlu;klu];

%thigh
ktemp_t = phip-pk;
jtemp_t = cross(ktemp_t,itemp_l);
itemp2_t = cross(jtemp_t,ktemp_t);

itu = itemp2_t./norm(itemp2_t);
jtu = jtemp_t./norm(jtemp_t);
ktu = ktemp_t./norm(ktemp_t);

%TA2G_T = [itu' jtu' ktu'];
TG2A_T = [itu;jtu;ktu];

% Anatomical-Marker Transformation Matrices
L1 = caldata(:,16:18);
L2 = caldata(:,19:21);
L3 = caldata(:,22:24);

T1 = caldata(:,7:9);
T2 = caldata(:,10:12);
T3 = caldata(:,13:15);

pA_ank = (TG2A_L*[pa(1)-cm_leg(1);pa(2)-cm_leg(2);pa(3)-cm_leg(3)])'; 
pA_kne = (TG2A_L*[pk(1)-cm_leg(1);pk(2)-cm_leg(2);pk(3)-cm_leg(3)])';

PA_L1 = (TG2A_L*[L1(1)-cm_leg(1);L1(2)-cm_leg(2);L1(3)-cm_leg(3)])';
PA_L2 = (TG2A_L*[L2(1)-cm_leg(1);L2(2)-cm_leg(2);L2(3)-cm_leg(3)])';
PA_L3 = (TG2A_L*[L3(1)-cm_leg(1);L3(2)-cm_leg(2);L3(3)-cm_leg(3)])';

pA_hip = (TG2A_T*[phip(1)-cm_thi(1);phip(2)-cm_thi(2);phip(3)-cm_thi(3)])';

PA_T1 = (TG2A_T*[T1(1)-cm_thi(1);T1(2)-cm_thi(2);T1(3)-cm_thi(3)])';
PA_T2 = (TG2A_T*[T2(1)-cm_thi(1);T2(2)-cm_thi(2);T2(3)-cm_thi(3)])';
PA_T3 = (TG2A_T*[T3(1)-cm_thi(1);T3(2)-cm_thi(2);T3(3)-cm_thi(3)])';

%leg
jtempL = PA_L2 - PA_L1;
atempL = PA_L3 - PA_L1;
itempL = cross(atempL,jtempL);
ktempL = cross(itempL,jtempL);

iL = itempL./norm(itempL);
jL = jtempL./norm(jtempL);
kL = ktempL./norm(ktempL);

TM2A_L = [iL;jL;kL];
TA2M_L = TM2A_L';

%thigh
jtempT = PA_T2 - PA_T1;
atempT = PA_T3 - PA_T1;
itempT = cross(atempT,jtempT);
ktempT = cross(itempT,jtempT);

iT = itempT./norm(itempT);
jT = jtempT./norm(jtempT);
kT = ktempT./norm(ktempT);

TM2A_T = [iT;jT;kT];
TA2M_T = TM2A_T';

%PART 2: Gait Trial
%step 1: Global_Marker Transformation Matrices

thighM1 = gait_data(:,1:3);
thighM2 = gait_data(:,4:6);
thighM3 = gait_data(:,7:9);

legM1 = gait_data(:,10:12);
legM2 = gait_data(:,13:15);
legM3 = gait_data(:,16:18);

for n = 1:length(legM1)
    
   jtemp = legM2(n,:) - legM1(n,:);
   atemp = legM3(n,:)-legM1(n,:);
   itemp = cross(atemp,jtemp);
   ktemp = cross(itemp,jtemp);

   i = itemp./norm(itemp);
   j = jtemp./norm(jtemp);
   k = ktemp./norm(ktemp);
   
   TG2M = [i;j;k];
   TG2A = TM2A_L*TG2M;
   TA2G = TG2A';
   idenM = [1 0 0;0 1 0;0 0 1];
   uL = TA2G*idenM;
   
   jtemp = thighM2(n,:) - thighM1(n,:);
   atemp = thighM3(n,:)-thighM1(n,:);
   itemp = cross(atemp,jtemp);
   ktemp = cross(itemp,jtemp);

   i = itemp./norm(itemp);
   j = jtemp./norm(jtemp);
   k = ktemp./norm(ktemp);
   
   TG2M = [i;j;k];
   TG2A = TM2A_T*TG2M;
   TA2G = TG2A';
   uT = TA2G*idenM;
   
   TL2T = uT*uL;
   
   beta(n) = asin(TL2T(3,1));
   gama(n) = acos(TL2T(1,1)/cos(beta(n)));
   alpha(n) = asin(-TL2T(3,2)/cos(beta(n)));
   
end

fs = 100;
dt = 1/fs;
time = dt:dt:length(alpha)*dt;

figure(1)
subplot(3,1,1)
plot(time,alpha,'k')
title('Knee Flexion-Extension Angle')
xlabel('time(sec)')
ylabel('F-E Angle(rad)')

subplot(3,1,2)
plot(time,beta,'k')
title('Knee Abduction-Adduction Angle')
xlabel('time(sec)')
ylabel('Ab-Ad Angle(rad)')

subplot(3,1,3)
plot(time,gama,'k')
title('Knee Internal-External rotation Angle')
xlabel('time(sec)')
ylabel('I-E rotation Angle(rad)')

% percentage gait cycle
heel_contact = min(alpha); %minimum flexion angle occurs at heel strike
%minimum flexion occurs at 102th data point ~ -0.7198
alpha_g = alpha(1:102);beta_g = beta(1:102);gama_g = gama(1:102);

%based on the collected data, let's assume that 1st point @alpha_g occurs at approximately 20% of the gait cyle
per_g = [20:80/102:100]; per_g(1)=[];


%Plot: rotations around 3 axes
figure(2)
subplot(3,1,1)
plot(per_g,alpha_g,'k')
title('% Gait Cycle: Flexion-Extension')
xlabel('% gait cycle')
ylabel('F-E Angle(rad)')

subplot(3,1,2)
plot(per_g,beta_g,'k')
title('% Gait Cycle: Abduction-Adduction')
xlabel('% gait cycle')
ylabel('Ab-Ad Angle(rad)')

subplot(3,1,3)
plot(per_g,gama_g,'k')
title('% Gait Cycle: Internal-External Rotation')
xlabel('% gait cycle')
ylabel('I-E rotation Angle(rad)')

%angular velocity
[vel_alp,acc_alp] = central_diff(alpha_g,fs);
[vel_bet,acc_bet] = central_diff(beta_g,fs);
[vel_gam,acc_gam] = central_diff(gama_g,fs);

%plots :Angular velocity
per_v = per_g(1:101);
figure(3)
subplot(3,1,1)
plot(per_v,vel_alp,'k')
title('Angular Velocity: Flexion-Extension')
xlabel('% gait cycle')
ylabel('F-E angular velo(rad/s)')

subplot(3,1,2)
plot(per_v,vel_bet,'k')
title('Angular Velocity: Abduction-Adduction')
xlabel('% gait cycle')
ylabel('Ab-Ad angular velo(rad/s)')

subplot(3,1,3)
plot(per_v,vel_gam,'k')
title('Angular Velocity: Internal-External rotation')
xlabel('% gait cycle')
ylabel('I-E angular velo(rad/s)')

%results
%range of motion
FE_range = max(alpha_g)-min(alpha_g);
AbAd_range = max(beta_g)-min(beta_g);
IE_range = max(gama_g)-min(gama_g);

% maximal velocities@Percentage of gait cycle
%this gait cycle is assumed to start at 20% of the gait cycle, based on the data collected
FE_maxV = 20+ 100*(find(vel_alp == max(vel_alp))/length(alpha_g));
AbAd_maxV = 20+ 100*(find(vel_bet == max(vel_bet))/length(alpha_g));
IE_maxV = 20 + 100*(find(vel_gam == max(vel_gam))/length(alpha_g));

