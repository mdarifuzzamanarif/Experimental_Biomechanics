%lab 5: Inverse Dynamics Lab
%Arif, Md. Arifuzzaman, UMKC ID: 16209626

clear
clc
close all
%%
load('T1_data.mat');
fs=100;

fx = fdata(:,1);fy = fdata(:,2);fz = fdata(:,3); %GRF
IFx = IF(1);IFy = IF(2);IFz = IF(3); %Moment of Inertia:Foot
ILx = IL(1);ILy = IL(2);ILz = IL(3); %Moment of Inertia:Leg
%%

for i=1:length(fx)
%step 1
%Foot     
fx=-fx;fy=-fy;fz=-fz; F_foot_distG = [fx(i) fy(i) fz(i)];
F_foot_px = mfoot*afoot(i,1) - fx(i);
F_foot_py = mfoot*afoot(i,2) - fy(i);
F_foot_pz = mfoot*afoot(i,3) - fz(i) + (mfoot*g);

F_foot_proxG = [F_foot_px F_foot_py F_foot_pz];

%Leg
F_leg_distG = -F_foot_proxG;
F_leg_px = mleg*aleg(i,1) - F_leg_distG(1);
F_leg_py = mleg*aleg(i,2) - F_leg_distG(2);
F_leg_pz = mleg*aleg(i,3) - F_leg_distG(3) + (mleg*g);

F_leg_proxG = [F_leg_px F_leg_py F_leg_pz];

%step2
%foot
jtempF = F2G(i,:) - F1G(i,:);
atempF = F3G(i,:) - F1G(i,:);
itempF = cross(atempF,jtempF);
ktempF = cross(itempF,jtempF);

iuF = itempF./norm(itempF);
juF = jtempF./norm(jtempF);
kuF = ktempF./norm(ktempF);

T_M2G_F = [iuF' juF' kuF'];
T_G2M_F = T_M2G_F';
T_G2A_F = TF_MtoA*T_G2M_F;
T_A2G_F = T_G2A_F';

%leg
jtempL = L2G(i,:) - L1G(i,:);
atempL = L3G(i,:) - L1G(i,:);
itempL = cross(atempL,jtempL);
ktempL = cross(itempL,jtempL);

iuL = itempL./norm(itempL);
juL = jtempL./norm(jtempL);
kuL = ktempL./norm(ktempL);

T_M2G_L = [iuL' juL' kuL'];
T_G2M_L = T_M2G_L';
T_G2A_L = TL_MtoA*T_G2M_L;
T_A2G_L = T_G2A_L';

%step 3
%foot
F_foot_distA = T_G2A_F*[F_foot_distG]';
F_foot_proxA = T_G2A_F*[F_foot_proxG]';
%leg
F_leg_distA = T_G2A_L*[F_leg_distG]';
F_leg_proxA = T_G2A_L*[F_leg_proxG]';

%step 4
Mx_foot_proxA = IFx*alpha_FOOT(i,1)+(IFz-IFy)*omega_FOOT(i,3)*omega_FOOT(i,2)- F_foot_distA(3)*lfoot_dist(i)- F_foot_proxA(3)*lfoot_prox(i);
My_foot_proxA = IFy*alpha_FOOT(i,2)+ (IFx-IFz)*omega_FOOT(i,1)*omega_FOOT(i,3);
Mz_foot_proxA = IFz*alpha_FOOT(i,3)+(IFy-IFx)*omega_FOOT(i,2)*omega_FOOT(i,1)- F_foot_distA(1)*lfoot_dist(i)- F_foot_proxA(1)*lfoot_prox(i);

%converting into leg anatomical coordinates
M_foot_proxG = T_A2G_F*[Mx_foot_proxA; My_foot_proxA; Mz_foot_proxA];
M_leg_distG = -M_foot_proxG;
M_leg_distA = T_G2A_L*M_leg_distG;

% Reaction moments at leg
Mx_leg_proxA = ILx*alpha_LEG(i,1)+(ILz - ILy)*omega_LEG(i,3)*omega_LEG(i,2)- F_leg_distA(2)*lleg_dist(i)- F_leg_proxA(2)*lleg_prox(i)- M_leg_distA(1);
My_leg_proxA = ILy*alpha_LEG(i,2)+(ILx - ILz)*omega_LEG(i,1)*omega_LEG(i,3)- F_leg_distA(1)*lleg_dist(i)- F_leg_proxA(1)*lleg_prox(i)- M_leg_distA(2);
Mz_leg_proxA = ILz*alpha_LEG(i,3)+(ILy - ILx)*omega_LEG(i,2)*omega_LEG(i,1)- M_leg_distA(3);
               
M_leg_proxG = T_A2G_L*[Mx_leg_proxA;My_leg_proxA;Mz_leg_proxA];   
end

%step 5: power
for i = 1:length(omega_LEG)
P_ank(i,:) = M_foot_proxG.*(omega_LEG(i,:)-omega_FOOT(i,:))';
P_kne(i,:) = M_leg_proxG.*(omega_THIGH(i,:)-omega_LEG(i,:))';
end
%step 6: Energy
% Ankle - X
igen_ankx = find(P_ank(:,1)>0);
iabs_ankx = find(P_ank(:,1)<0);

Pgen_ankx = zeros(size(P_ank(:,1)));
Pgen_ankx(igen_ankx,:) = P_ank(igen_ankx,1);

Pabs_ankx = zeros(size(P_ank(:,1)));
Pabs_ankx(iabs_ankx,:) = P_ank(iabs_ankx,1);

Wgen_ankx = (1/fs)*trapz(Pgen_ankx);
Wabs_ankx = (1/fs)*trapz(Pabs_ankx);

% Ankle - Y
igen_anky = find(P_ank(:,2)>0);
iabs_anky = find(P_ank(:,2)<0);

Pgen_anky = zeros(size(P_ank(:,1)));
Pgen_anky(igen_anky,:) = P_ank(igen_anky,1);

Pabs_anky = zeros(size(P_ank(:,1)));
Pabs_anky(iabs_anky,:) = P_ank(iabs_anky,1);

Wgen_anky = (1/fs)*trapz(Pgen_anky);
Wabs_anky = (1/fs)*trapz(Pabs_anky);

% Ankle - Z
igen_ankz = find(P_ank(:,3)>0);
iabs_ankz = find(P_ank(:,3)<0);

Pgen_ankz = zeros(size(P_ank(:,1)));
Pgen_ankz(igen_ankz,:) = P_ank(igen_ankz,1);

Pabs_ankz = zeros(size(P_ank(:,1)));
Pabs_ankz(iabs_ankz,:) = P_ank(iabs_ankz,1);

Wgen_ankz = (1/fs)*trapz(Pgen_ankz);
Wabs_ankz = (1/fs)*trapz(Pabs_ankz);

%Knee - X
igen_knex = find(P_kne(:,1)>0);
iabs_knex = find(P_kne(:,1)<0);

Pgen_knex = zeros(size(P_kne(:,1)));
Pgen_knex(igen_knex,:) = P_kne(igen_knex,1);

Pabs_knex = zeros(size(P_kne(:,1)));
Pabs_knex(iabs_knex,:) = P_kne(iabs_knex,1);

Wgen_knex = (1/fs)*trapz(Pgen_knex);
Wabs_knex = (1/fs)*trapz(Pabs_knex);

%Knee - Y
igen_kney = find(P_kne(:,2)>0);
iabs_kney = find(P_kne(:,2)<0);

Pgen_kney = zeros(size(P_kne(:,1)));
Pgen_kney(igen_kney,:) = P_kne(igen_kney,1);

Pabs_kney = zeros(size(P_kne(:,1)));
Pabs_kney(iabs_kney,:) = P_kne(iabs_kney,1);

Wgen_kney = (1/fs)*trapz(Pgen_kney);
Wabs_kney = (1/fs)*trapz(Pabs_kney);

%Knee - Z
igen_knez = find(P_kne(:,3)>0);
iabs_knez = find(P_kne(:,3)<0);

Pgen_knez = zeros(size(P_kne(:,1)));
Pgen_knez(igen_knez,:) = P_kne(igen_knez,1);

Pabs_knez = zeros(size(P_kne(:,1)));
Pabs_knez(iabs_knez,:) = P_kne(iabs_knez,1);

Wgen_knez = (1/fs)*trapz(Pgen_knez);
Wabs_knez = (1/fs)*trapz(Pabs_knez);

%%
%plotting
t = (1/fs)*(1:1:length(P_ank));

subplot(2,1,1)
plot(t,P_kne(:,1),t,P_kne(:,2),t,P_kne(:,3))
xlabel('Time (s)','FontWeight','bold')
ylabel('Power - Knee Joint (W)','FontWeight','bold')
legend('ML(x)','AP(y)','Vertical(z)')
title('Dual Limb Squat','FontWeight','bold')
%axis([min(t) max(t) -4E4 4E4])

subplot(2,1,2)
plot(t,P_ank(:,1),t,P_ank(:,2),t,P_ank(:,3))
xlabel('Time (s)','FontWeight','bold')
ylabel('Power - Ankle Joint (W)','FontWeight','bold')
legend('ML(x)','AP(y)','Vertical(z)')
%axis([min(t) max(t) -1.5E4 1.5E4])