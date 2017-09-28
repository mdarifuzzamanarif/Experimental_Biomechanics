%lab 3: Kinetics:Force Plate Lab
%Arif, Md. Arifuzzaman
clear
clc
close all
fs=1000;
%stance: trial 1,2,3
%step 1: Determine foot positions
load('trial1.mat');
load footpositions.mat

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %left leg
% pankx_L=p_left_ank_med(:,1) + (p_left_ank_lat(:,1)-p_left_ank_med(:,1))./2;
% panky_L=p_left_ank_med(:,2) + (p_left_ank_lat(:,2)-p_left_ank_med(:,2))./2;
% pankz_L=p_left_ank_med(:,3) + (p_left_ank_lat(:,3)-p_left_ank_med(:,3))./2;
% 
% 
% ptoex_L=p_left_toe_med(:,1) + (p_left_toe_lat(:,1)-p_left_toe_med(:,1))./2;
% ptoey_L=p_left_toe_med(:,2) + (p_left_toe_lat(:,2)-p_left_toe_med(:,2))./2;
% ptoez_L=p_left_toe_med(:,3) + (p_left_toe_lat(:,3)-p_left_toe_med(:,3))./2;
% 
% cmx_L=pankx_L + (ptoex_L-pankx_L)./2;
% cmy_L=panky_L + (ptoey_L-panky_L)./2;
% cmz_L=pankz_L + (ptoez_L-pankz_L)./2;
% %right leg
% pankx_R=p_right_ank_med(:,1) + (p_right_ank_lat(:,1)-p_right_ank_med(:,1))./2;
% panky_R=p_right_ank_med(:,2) + (p_right_ank_lat(:,2)-p_right_ank_med(:,2))./2;
% pankz_R=p_right_ank_med(:,3) + (p_right_ank_lat(:,3)-p_right_ank_med(:,3))./2;
% 
% ptoex_R=p_right_toe_med(:,1) + (p_right_toe_lat(:,1)-p_right_toe_med(:,1))./2;
% ptoey_R=p_right_toe_med(:,2) + (p_right_toe_lat(:,2)-p_right_toe_med(:,2))./2;
% ptoez_R=p_right_toe_med(:,3) + (p_right_toe_lat(:,3)-p_right_toe_med(:,3))./2;
% 
% cmx_R=pankx_R + (ptoex_R-pankx_R)./2;
% cmy_R=panky_R + (ptoey_R-panky_R)./2;
% cmz_R=pankz_R + (ptoez_R-pankz_R)./2;
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %1. itemp
% itempx_L = pankx_L - p_left_ank_med(:,1); %for left foot,medial malleolus
% itempy_L = panky_L - p_left_ank_med(:,2);
% itempz_L = pankz_L - p_left_ank_med(:,3);
% 
% itempx_R = pankx_R - p_right_ank_lat(:,1); %for right foot,lateral malleolus
% itempy_R = panky_R - p_right_ank_lat(:,2);
% itempz_R = pankz_R - p_right_ank_lat(:,3);
% %2. jtemp
% jtempx_L = pankx_L - ptoex_L;
% jtempy_L = panky_L - ptoey_L;
% jtempz_L = pankz_L - ptoez_L;
% 
% jtempx_R = pankx_R - ptoex_R;
% jtempy_R = panky_R - ptoey_R;
% jtempz_R = pankz_R - ptoez_R;
% %3. temporary ktemp
% itemp_L = [itempx_L itempy_L itempz_L];
% jtemp_L = [jtempx_L jtempy_L jtempz_L];
% 
% itemp_R = [itempx_R itempy_R itempz_R];
% jtemp_R = [jtempx_R jtempy_R jtempz_R];
% 
% ktemp_L = cross(itemp_L,jtemp_L);
% ktemp_R = cross(itemp_R,jtemp_R);
% 
% %4. itemp2
% itemp2_L = cross(jtemp_L,ktemp_L);
% itemp2_R = cross(jtemp_R,ktemp_R);
% 
% %5. unit vectors
% itemp2U_L = itemp2_L./norm(itemp2_L);
% itemp2U_R = itemp2_R./norm(itemp2_R);
% 
% jtempU_L = jtemp_L./norm(jtemp_L);
% jtempU_R = jtemp_R./norm(jtemp_R);
% 
% ktempU_L = ktemp_L./norm(ktemp_L);
% ktempU_R = ktemp_R./norm(ktemp_R);
% 
% for a=1:3000
% %6. Transformation matrix
% TP2G_L = [(itemp2U_L(a,:))',(jtempU_L(a,:))',(ktempU_L(a,:))'];
% TP2G_R = [(itemp2U_R(a,:))',(jtempU_R(a,:))',(ktempU_R(a,:))'];
% 
% cm_L = [cmx_L(a,1) cmy_L(a,1) cmz_L(a,1)];
% 
% for b=1:980
%    temporary = cm_L' + (TP2G_L*pleft(b,:)');
%    PfG_R{a,1}(b,:) = temporary';
% % PfG_R{b,a} = cm_L' + (TP2G_L*pleft(b,:)');
% 
% % PfG_Ry = cmy_L + (TP2G_L.*pleft(a,2));
% % PfG_Rz = cmz_L + (TP2G_L.*pleft(a,3));
% 
% end
% end
for i = 1:length(p_left_ank_lat)
   %joint center: right
    pank_r = p_right_ank_med(i,:) + (p_right_ank_lat(i,:)-p_right_ank_med(i,:))/2; 
    ptoe_r = p_right_toe_med(i,:) + (p_right_toe_lat(i,:)-p_right_toe_med(i,:))/2;
    %center of mass: right
    cm_r = pank_r + (ptoe_r-pank_r)/2; 
    %joint center: left
    pank_l = p_left_ank_lat(i,:) + (p_left_ank_med(i,:)-p_left_ank_lat(i,:))/2; 
    ptoe_l = p_left_toe_lat(i,:) + (p_left_toe_med(i,:)-p_left_toe_lat(i,:))/2;
    %center of mass: left
    cm_l = pank_l + (ptoe_l-pank_l)/2; 
    %right foot anatomical coordinates
    itempR = p_right_ank_lat(i,:) - pank_r;
    jtempR = ptoe_r - pank_r;
    ktempR = cross(itempR,jtempR);
    itemp2R = cross(jtempR,ktempR);
    %unit vectors: right
    iR = itemp2R./norm(itemp2R); 
    jR = jtempR./norm(jtempR);
    kR = ktempR./norm(ktempR);
    %transformation mat: right
    TR = [iR;jR;kR]'; 
    %left foot anatomical coordinates
    itempL = -p_left_ank_lat(i,:) + pank_l; 
    jtempL = ptoe_l - pank_l;
    ktempL = cross(itempL,jtempL);
    itemp2L = cross(jtempL,ktempL);
    %unit vectors: left
    iL = itemp2L./norm(itemp2L); 
    jL = jtempL./norm(jtempL);
    kL = ktempL./norm(ktempL);
    %transformation mat: left
    TL = [iL;jL;kL]'; 
    %%%%%%%%%%%%%%%%%
    for j = 1:length(pleft)
        pfgR = cm_r'+ TR*[pright(j,:)'];
        pfgR1(i,j,:) = pfgR';
        pfgL = cm_l'+ TL*[pleft(j,:)'];
        pfgL1(i,j,:) = pfgL';
    end
end
i1 = find(pleft(:,1));
i2 = find(pleft(:,2));
i3 = find(pright(:,1));
i4 = find(pright(:,2));

pf1 = pfgR1(i1(1),:,:);
pf2 = pfgL1(i2(1),:,:);
pf3 = pfgR1(i3(3),:,:);
pf4 = pfgR1(i4(3),:,:);
%%step 2 :center of pressure
az = 3.8;
copxL = (M2(:,2)-(F2(:,1)*az))./-F2(:,3);   %left leg: FP 2: F2, M2
copyL = (M2(:,1) + (F2(:,2)*az))./F2(:,3); 
copxR = (M3(:,2)-(F3(:,1)*az))./-F3(:,3);   %right leg: FP 3: F3, M3 
copyR = (M3(:,1) + (F3(:,2)*az))./F3(:,3); 
% Offsets
fp_2x = -27; fp_2y = 895; fp_3x = 488; fp_3y = 895; %force plate 1: Left leg, force plate 2:right leg
copxL = copxL +fp_2x; copyL = copyL +fp_2y; copxR = copxR +fp_3x; copyR = copyR +fp_3y;
%COP
copL = [copxL copyL]; %force plate 2
copR = [copxR copyR]; %force plate 3
%net COP
cop1(1:3000,1:2)=0;
cop4(1:3000,1:2)=0;
cpnet = copnet4(cop1,copL,copR,cop4,F1(:,3),F2(:,3),F3(:,3),F4(:,3));

%step 3:center of mass
F = [(F2(:,1)+F3(:,1)),(F2(:,2)+F3(:,2))]; %medial lateral and anterior posterior respectively
W = abs(mean(F2(:,3) + F3(:,3)));
[CPout,CMout,RMSE] = glp(F,cpnet,W,fs);
RMSE = mean(RMSE);

%step 4: COP and COM movement quantification
%pathlengths
for i = 1:length(cpnet)-1
    %left foot
    P1L = copL(i+1,:);
    P2L = copL(i,:);
    dL = sqrt((P2L(1)-P1L(1))^2 + ((P2L(2)-P1L(2))^2));
    PLcopL(i) = dL;
    %right foot
    P1R = copR(i+1,:);
    P2R = copR(i,:);
    dR = sqrt((P2R(1)-P1R(1))^2 + ((P2R(2)-P1R(2))^2));
    PLcopR(i) = dR;
end
PLcopL = sum(PLcopL);
PLcopR = sum(PLcopR);
for j = 1:length(CMout)-1
    P1 = CMout(j+1,:)  ;
    P2 = CMout(j,:)    ; 
    dM = sqrt((P2(1)-P1(1))^2 + ((P2(2)-P1(2))^2));
    PLcom(j) = dM;
end
PLcom = sum(PLcom);
% Foot length
fl_l = max(pleft(:,2))-min(pleft(:,2));
fl_r = max(pright(:,2))-min(pright(:,2));
fl = (fl_l+fl_r)/2;

pl_COPL = 100*(PLcopL/fl_l);
pl_COPR = 100*(PLcopR/fl_r);
Pl_COM = 100*(PLcom/fl);

%Range of Motion:(A-P)~y
ROMcopL = 100*(max(copL(:,2))-min(copL(:,2)))/fl_l;
ROMcopR = 100*(max(copR(:,2))-min(copR(:,2)))/fl_r ;
ROM_COM =  100*(max(CMout(:,2))-min(CMout(:,2)))/fl;
% Distance from Base of Support
bos_l = max(pleft(:,2)) - max(copL(:,2)); bos_l = 100*(bos_l/fl_l);
bos_r = max(pright(:,2)) - max(copR(:,2)); bos_r = 100*(bos_r/fl_r);
bos_COM = ((max(pleft(:,2))+ max(pright(:,2)))/2) - (max(CMout(:,2))+mean(cpnet(:,2)));
bos_COM = 100*(bos_COM/fl);


%plotting
t = (1/1000)*(1:1:length(CPout));
figure(1)
%plot(cpnet(:,1),cpnet(:,2),CMout(:,1)+mean(cpnet(:,1)),CMout(:,2)+mean(cpnet(:,2)), copL(:,1),copL(:,2),copR(:,1),copR(:,2)    ,pleft(:,1),pleft(:,2),'k',pright(:,1),pright(:,2),'k')
plot(cpnet(:,1),cpnet(:,2),CMout(:,1)+mean(cpnet(:,1)),CMout(:,2)+mean(cpnet(:,2)), copL(:,1),copL(:,2),copR(:,1),copR(:,2)    ,pf1(:,:,1),pf1(:,:,2),pf2(:,:,1),pf2(:,:,2),pf3(:,:,1),pf3(:,:,2))
legend('net COP','COM','L COP','R COP')
xlabel('M-L Positions(mm)')
ylabel('A-P Positions (mm)')
%figure(2)
%plot(pleft(:,1),pleft(:,2),'k',pright(:,1),pright(:,2),'k')

