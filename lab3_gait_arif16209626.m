%lab 3: Kinetics:Force Plate Lab
clear
clc
close all
fs=1000;

%step 1: Determine foot positions
load('trial1.mat');
load footpositions.mat

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%step 1
%left leg
pankx_L=p_left_ank_med(:,1) + (p_left_ank_lat(:,1)-p_left_ank_med(:,1))./2;
panky_L=p_left_ank_med(:,2) + (p_left_ank_lat(:,2)-p_left_ank_med(:,2))./2;
pankz_L=p_left_ank_med(:,3) + (p_left_ank_lat(:,3)-p_left_ank_med(:,3))./2;


ptoex_L=p_left_toe_med(:,1) + (p_left_toe_lat(:,1)-p_left_toe_med(:,1))./2;
ptoey_L=p_left_toe_med(:,2) + (p_left_toe_lat(:,2)-p_left_toe_med(:,2))./2;
ptoez_L=p_left_toe_med(:,3) + (p_left_toe_lat(:,3)-p_left_toe_med(:,3))./2;

cmx_L=pankx_L + (ptoex_L-pankx_L)./2;
cmy_L=panky_L + (ptoey_L-panky_L)./2;
cmz_L=pankz_L + (ptoez_L-pankz_L)./2;
%right leg
pankx_R=p_right_ank_med(:,1) + (p_right_ank_lat(:,1)-p_right_ank_med(:,1))./2;
panky_R=p_right_ank_med(:,2) + (p_right_ank_lat(:,2)-p_right_ank_med(:,2))./2;
pankz_R=p_right_ank_med(:,3) + (p_right_ank_lat(:,3)-p_right_ank_med(:,3))./2;

ptoex_R=p_right_toe_med(:,1) + (p_right_toe_lat(:,1)-p_right_toe_med(:,1))./2;
ptoey_R=p_right_toe_med(:,2) + (p_right_toe_lat(:,2)-p_right_toe_med(:,2))./2;
ptoez_R=p_right_toe_med(:,3) + (p_right_toe_lat(:,3)-p_right_toe_med(:,3))./2;

cmx_R=pankx_R + (ptoex_R-pankx_R)./2;
cmy_R=panky_R + (ptoey_R-panky_R)./2;
cmz_R=pankz_R + (ptoez_R-pankz_R)./2;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%1. itemp
itempx_L = pankx_L - p_left_ank_med(:,1); %for left foot,medial malleolus
itempy_L = panky_L - p_left_ank_med(:,2);
itempz_L = pankz_L - p_left_ank_med(:,3);

itempx_R = pankx_R - p_right_ank_lat(:,1); %for right foot,lateral malleolus
itempy_R = panky_R - p_right_ank_lat(:,2);
itempz_R = pankz_R - p_right_ank_lat(:,3);
%2. jtemp
jtempx_L = pankx_L - ptoex_L;
jtempy_L = panky_L - ptoey_L;
jtempz_L = pankz_L - ptoez_L;

jtempx_R = pankx_R - ptoex_R;
jtempy_R = panky_R - ptoey_R;
jtempz_R = pankz_R - ptoez_R;
%3. temporary ktemp
itemp_L = [itempx_L itempy_L itempz_L];
jtemp_L = [jtempx_L jtempy_L jtempz_L];

itemp_R = [itempx_R itempy_R itempz_R];
jtemp_R = [jtempx_R jtempy_R jtempz_R];

ktemp_L = cross(itemp_L,jtemp_L);
ktemp_R = cross(itemp_R,jtemp_R);

%4. itemp2
itemp2_L = cross(jtemp_L,ktemp_L);
itemp2_R = cross(jtemp_R,ktemp_R);

%5. unit vectors
itemp2U_L = itemp2_L./norm(itemp2_L);
itemp2U_R = itemp2_R./norm(itemp2_R);

jtempU_L = jtemp_L./norm(jtemp_L);
jtempU_R = jtemp_R./norm(jtemp_R);

ktempU_L = ktemp_L./norm(ktemp_L);
ktempU_R = ktemp_R./norm(ktemp_R);

for a=1:3000
%6. Transformation matrix
TP2G_L = [(itemp2U_L(a,:))',(jtempU_L(a,:))',(ktempU_L(a,:))'];
TP2G_R = [(itemp2U_R(a,:))',(jtempU_R(a,:))',(ktempU_R(a,:))'];

cm_L = [cmx_L(a,1) cmy_L(a,1) cmz_L(a,1)];

for b=1:980
   temporary = cm_L' + (TP2G_L*pleft(b,:)');
   PfG_R{a,1}(b,:) = temporary';
% PfG_R{b,a} = cm_L' + (TP2G_L*pleft(b,:)');

% PfG_Ry = cmy_L + (TP2G_L.*pleft(a,2));
% PfG_Rz = cmz_L + (TP2G_L.*pleft(a,3));

end
end

%%step 2 : calculate center of pressure time story
%center of pressure
az = 3.8;
copxL = (M2(:,2)-(F2(:,1)*az))./-F2(:,3);   %left leg: FP 2: F2, M2
copyL = (M2(:,1) + (F2(:,2)*az))./F2(:,3); 
copxR = (M3(:,2)-(F3(:,1)*az))./-F3(:,3);   %right leg: FP 3: F3, M3 
copyR = (M3(:,1) + (F3(:,2)*az))./F3(:,3); 

% Offsets
FP2X = -27; FP2Y = 895; FP3X = 488; FP3Y = 895; %force plate 1: Left leg, force plate 2:right leg

copxL = copxL +FP2X; copyL = copyL +FP2Y; copxR = copxR +FP3X; copyR = copyR +FP3Y;

%filtering
copxL = filter_data(copxL,30,fs,[1]); 
copyL = filter_data(copyL,30,fs,[1]); 
copxR = filter_data(copxR,30,fs,[1]); 
copyR = filter_data(copyR,30,fs,[1]); 

copL = [copxL copyL]; %force plate 2
copR = [copxR copyR]; %force plate 3

cop1(1:3000,1:2)=0;
cop4(1:3000,1:2)=0;
cpnet = copnet4(cop1,copL,copR,cop4,F1(:,3),F2(:,3),F3(:,3),F4(:,3));

%step 3: calculate center of mass time history

F = [(F2(:,1)+F3(:,1)),(F2(:,2)+F3(:,2))]; %medial lateral and anterior posterior respectively
W = abs(mean(F2(:,3) + F3(:,3)));

[CPout,CMout,RMSE] = glp(F,cpnet,W,fs);
RMSE = mean(RMSE);

%step 4: COP and COM movement quantification
%pathlengths

for i = 1:length(cpnet)-1
    
    P1L = copL(i+1,:);
    P2L = copL(i,:);
    dL = sqrt((P2L(1)-P1L(1))^2 + ((P2L(2)-P1L(2))^2));
    PLcopL(i) = dL;
    
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
footlength_left = max(pleft(:,2))-min(pleft(:,2));
footlength_right = max(pright(:,2))-min(pright(:,2));
footlength = (footlength_left+footlength_right)/2;

PathLengthCOP_Left = 100*(PLcopL/footlength_left);
PathLengthCOP_Right = 100*(PLcopR/footlength_right);

PathLength_COM = 100*(PLcom/footlength);

%Range of Motion: anterior-posterior

LCOP_ROM = 100*(max(copL(:,2))-min(copL(:,2)))/footlength_left;
RCOP_ROM = 100*(max(copR(:,2))-min(copR(:,2)))/footlength_right ;
COM_ROM =  100*(max(CMout(:,2))-min(CMout(:,2)))/footlength;

% Distance from Base of Support

BOSLeft = max(pleft(:,2)) - max(copL(:,2)); BOSLeft = 100*BOSLeft/footlength_left;
BOSRight = max(pright(:,2)) - max(copR(:,2)); BOSRight = 100*BOSRight/footlength_right;
BOSCOM = ((max(pleft(:,2))+ max(pright(:,2)))/2) - (max(CMout(:,2))+mean(cpnet(:,2)));
BOSCOM = 100*BOSCOM/footlength;

%plotting
t = (1/1000)*(1:1:length(CPout));
figure(1)
plot(cpnet(:,1),cpnet(:,2),CMout(:,1)+mean(cpnet(:,1)),CMout(:,2)+mean(cpnet(:,2)), copL(:,1),copL(:,2),copR(:,1),copR(:,2))
legend('Net COP','COM','Left COP','Right COP')
xlabel('Medial Lateral Positions (mm)')
ylabel('Anterior Posterior Positions (mm)')


%figure(2)
%plot(pleft(:,1),pleft(:,2),'k',pright(:,1),pright(:,2),'k')