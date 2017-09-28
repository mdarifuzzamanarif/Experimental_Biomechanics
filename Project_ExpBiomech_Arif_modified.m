%project:Exp Biomechanics
%single leg drop jump

clear
clc
close all

fs=1000; %sampling frequency
dt = 1/fs; %time step

%load data
for n=1:43
filename=(['trial',num2str(n),'.csv']);   %filename to be defined according to the file index
temp{n}=csvread(filename,5,2);
end

%data orientation
%           6704                          5670
%fx fy fz mx my mz cx cy cz   fx fy fz mx my mz cx cy cz 
% 1  2  3  4  5  6  7  8  9   10 11 12 13 14 15 16 17 18

%           6703                          6700
%fx fy fz mx my mz cx cy cz   fx fy fz mx my mz cx cy cz 
%19 20 21 22 23 24 25 26 27   28 29 30 31 32 33 34 35 36

for m=1:43
data=temp{m};
% t = [1/fs : 1/fs : length(data)/fs];
cop1 = data(:,7:8);
cop2 = data(:,16:17);
avg_cop1(m,:)=mean(abs(cop1));
avg_cop2(m,:)=mean(abs(cop2));
%rms_cop1(m,:)= rms(cop1(:,1));
RMS(m,:)= rms(cop1);
max_cop1(m,:)=max(cop1);
min_cop1(m,:)=min(cop1);
var_max(m,:)=(max_cop1(m,:)-avg_cop1(m,:));
var_min(m,:)=(avg_cop1(m,:)-min_cop1(m,:));

%[vel_cop,acc_cop] = derivative(cop1,1,dt);
% for p=1:2
% [velocity,acceleration]=central_diff(cop1(:,p),fs);
% velocity=velocity';
% vel_cop(:,p)=velocity;
% end
[velo_x,acc_x]=central_diff(cop1(:,1),fs);
[velo_y,acc_y]=central_diff(cop1(:,2),fs);
velo=[velo_x' velo_y'];
v{m}=velo;
dat=v{m};
vel(m,:)=mean(abs(dat(:,:)));

%COP and COM movement quantification
%pathlengths

for i = 1:length(data)-1
    
    P1L = cop1(i+1,:);
    P2L = cop1(i,:);
    dL = sqrt((P2L(1)-P1L(1))^2 + ((P2L(2)-P1L(2))^2));
    PLcopL(i) = dL;
    %PLcopL(m) = sum(PLcopL);   
end
PLcopL=PLcopL';
PLcopL_sum = sum(PLcopL);
t=length(data)*dt;
COP_velocity(m) = PLcopL_sum/t;

% %power spectrum
[freq,Pow] = power_spectrum(cop1(:,1),fs);
freq=freq';
%MF = medianfreq(fx,Px);


end

%MF = medianfreq(freq,Pow);

%PLcopL_sum = (PLcopL_sum)';
COP_velocity = COP_velocity';
for a=2:41
    for b=1:2
    avgt_cop1(a,b)= (avg_cop1(a,b)+avg_cop1(a+1,b)+avg_cop1(a+2,b))/3;
    avgt_vel1(a,b)= (vel(a,b)+vel(a+1,b)+vel(a+2,b))/3;
    avgt_RMS1(a,b)= (RMS(a,b)+RMS(a+1,b)+RMS(a+2,b))/3;
    
    var_max_avg(a,b)= (var_max(a,b)+var_max(a+1,b)+var_max(a+2,b))/3;
    var_min_avg(a,b)= (var_min(a,b)+var_min(a+1,b)+var_min(a+2,b))/3;
    
    end
    COP_velocity_avg(a,1) = (COP_velocity(a,1)+COP_velocity(a+1,1)+COP_velocity(a+2,1))/3;
    %avgt_cop1(a,1)= (avg_cop1(a,1)+avg_cop1(a+1,1)+avg_cop1(a+2,1))/3;
%     a=a+3;
%     if a==43
%         break
%     else continue
%     end
%     end
end

avg_cop_trials=[avgt_cop1(2,:);avgt_cop1(5,:);avgt_cop1(8,:);avgt_cop1(11,:);avgt_cop1(14,:);avgt_cop1(17,:);avgt_cop1(20,:);avgt_cop1(23,:);avgt_cop1(26,:);avgt_cop1(29,:)];
distance=[15;20;25;30;35;40;45;50;54;58];
figure(1)
subplot(2,1,1)
%plot(distance(2:10,1),avg_cop_trials(2:10,1),'*')
x=distance(2:10,1);
y=avg_cop_trials(2:10,1);
% p = polyfit(x,y,4);
% x2 = 15:0.5:60;
% y2 = polyval(p,x2);
% plot(x,y,'o',x2,y2)
plot(x,y)
title('Average COP (M-L) Vs Distance');
xlabel('Distance(mm)');
ylabel('Avg. COP movement(M-L)(mm)');

subplot(2,1,2)
%plot(distance(2:10,1),avg_cop_trials(2:10,2),'*')
x=distance(2:10,1);
y=avg_cop_trials(2:10,2);
% p = polyfit(x,y,4);
% x2 = 15:0.5:60;
% y2 = polyval(p,x2);
% plot(x,y,'o',x2,y2)
plot(x,y)
title('Average COP (A-P) Vs Distance');
xlabel('Distance(mm)');
ylabel('Avg. COP movement(A-P)(mm)');

avg_vel_trials=[avgt_vel1(2,:);avgt_vel1(5,:);avgt_vel1(8,:);avgt_vel1(11,:);avgt_vel1(14,:);avgt_vel1(17,:);avgt_vel1(20,:);avgt_vel1(23,:);avgt_vel1(26,:);avgt_vel1(29,:)];
figure(2)
subplot(2,1,1)
%plot(distance(2:10,1),avg_vel_trials(2:10,1),'x')
x=distance(2:10,1);
y=avg_vel_trials(2:10,1);
% p = polyfit(x,y,4);
% x2 = 15:0.5:60;
% y2 = polyval(p,x2);
% plot(x,y,'o',x2,y2)
plot(x,y)
title('Average Velocity (M-L) Vs Distance');
xlabel('Distance(mm)');
ylabel('Avg. COP movement(M-L)(mm/s)');

subplot(2,1,2)
%plot(distance(2:10,1),avg_vel_trials(2:10,2),'x')
x=distance(2:10,1);
y=avg_vel_trials(2:10,2);
% p = polyfit(x,y,4);
% x2 = 15:0.5:60;
% y2 = polyval(p,x2);
% plot(x,y,'o',x2,y2)
plot(x,y)
title('Average Velocity (A-P) Vs Distance');
xlabel('Distance(mm/s)');
ylabel('Avg. COP movement(A-P)(mm/s)');


avg_RMS_trials=[avgt_RMS1(2,:);avgt_RMS1(5,:);avgt_RMS1(8,:);avgt_RMS1(11,:);avgt_RMS1(14,:);avgt_RMS1(17,:);avgt_RMS1(20,:);avgt_RMS1(23,:);avgt_RMS1(26,:);avgt_RMS1(29,:)];
figure(3)
subplot(2,1,1)
%plot(distance(2:10,1),avg_RMS_trials(2:10,1),'x')
x=distance(2:10,1);
y=avg_RMS_trials(2:10,1);
% p = polyfit(x,y,4);
% x2 = 15:0.5:60;
% y2 = polyval(p,x2);
% plot(x,y,'o',x2,y2)
plot(x,y)
title('Average RMS Amplitude COP (M-L) Vs Distance');
xlabel('Distance(mm)');
ylabel('Avg. COP movement(M-L)(mm)');

subplot(2,1,2)
%plot(distance(2:10,1),avg_RMS_trials(2:10,2),'x')
x=distance(2:10,1);
y=avg_RMS_trials(2:10,2);
% p = polyfit(x,y,4);
% x2 = 15:0.5:60;
% y2 = polyval(p,x2);
% plot(x,y,'o',x2,y2)
plot(x,y)
title('Average RMS Amplitude COP (A-P) Vs Distance');
xlabel('Distance(mm)');
ylabel('Avg. COP movement(A-P)(mm)');

COP_velocity_trials=[COP_velocity_avg(2,:);COP_velocity_avg(5,:);COP_velocity_avg(8,:);COP_velocity_avg(11,:);COP_velocity_avg(14,:);COP_velocity_avg(17,:);COP_velocity_avg(20,:);COP_velocity_avg(23,:);COP_velocity_avg(26,:);COP_velocity_avg(29,:)];
figure(4)
% plot(distance(:,1),COP_velocity_trials(:,1),'o')
% hold on
% coef_fit=polyfit(distance(:,1),COP_velocity_trials(:,1),4);
% y_fit = polyval(coef_fit,[0 60]);
% plot([0 60],y_fit,'r')
x=distance(:,1);
y=COP_velocity_trials(:,1);
% p = polyfit(x,y,4);
% x2 = 15:0.5:60;
% y2 = polyval(p,x2);
% plot(x,y,'o',x2,y2)
plot(x,y)
title('Average COP Velocity Vs Distance');
xlabel('Distance(mm)');
ylabel('Avg. COP movement(mm/s)');

%plotting: Eyes Open Vs Eyes Closed

avg_cop_EO=[avgt_cop1(2,:);avgt_cop1(5,:);avgt_cop1(8,:);avgt_cop1(11,:)];
avg_cop_EC=[avgt_cop1(32,:);avgt_cop1(35,:);avgt_cop1(38,:);avgt_cop1(41,:)];
distanceE=[15;20;25;30];
figure(5)
subplot(2,1,1)
%plot(distanceE,avg_cop_EO(:,1),'*',distanceE,avg_cop_EC(:,1),'o')
x=distanceE(:,1);
y=avg_cop_EO(:,1);
z=avg_cop_EC(:,1);
% p1 = polyfit(x,y,4);
% p2 = polyfit(x,z,4);
% x2 = 15:0.5:30;
% y2 = polyval(p1,x2);
% z2 = polyval(p2,x2);
% plot(x,y,'*',x2,y2,'o')
% hold on
%plot(x,z,'o',x2,z2,'r')
plot(x,y,x,z)
title('EO Vs. EC: Average COP(M-L)');
xlabel('Distance(mm)');
ylabel('Avg. COP movement(M-L)(mm)');
legend('EO','EC');

subplot(2,1,2)
%plot(distanceE,avg_cop_EO(:,2),'*',distanceE,avg_cop_EC(:,2),'o')
x=distanceE(:,1);
y=avg_cop_EO(:,2);
z=avg_cop_EC(:,2);
% p1 = polyfit(x,y,4);
% p2 = polyfit(x,z,4);
% x2 = 15:0.5:30;
% y2 = polyval(p1,x2);
% z2 = polyval(p2,x2);
% plot(x,y,'*',x2,y2,'o')
% hold on
% plot(x,z,'o',x2,z2,'r')
plot(x,y,x,z)

title('EO Vs. EC: Average COP (A-P)');
xlabel('Distance(mm)');
ylabel('Avg. COP movement(A-P)(mm)');
legend('EO','EC');


figure(6)
avg_vel_EO=[avgt_vel1(2,:);avgt_vel1(5,:);avgt_vel1(8,:);avgt_vel1(11,:)];
avg_vel_EC=[avgt_vel1(32,:);avgt_vel1(35,:);avgt_vel1(38,:);avgt_vel1(41,:)];
subplot(2,1,1)
%plot(distanceE(:,1),avg_vel_EO(:,1),'*',distanceE(:,1),avg_vel_EC(:,1),'o')
x=distanceE(:,1);
y=avg_vel_EO(:,1);
z=avg_vel_EC(:,1);
% p1 = polyfit(x,y,4);
% p2 = polyfit(x,z,4);
% x2 = 15:0.5:30;
% y2 = polyval(p1,x2);
% z2 = polyval(p2,x2);
% plot(x,y,'*',x2,y2,'o')
% hold on
% plot(x,z,'o',x2,z2,'r')
plot(x,y,x,z)

title('EO Vs. EC: Average COP Velocity (M-L)');
xlabel('Distance(mm)');
ylabel('Avg. Velocity(M-L)(mm/s)');
legend('EO','EC');


subplot(2,1,2)
%plot(distanceE(:,1),avg_vel_EO(:,2),'*',distanceE(:,1),avg_vel_EC(:,2),'o')
x=distanceE(:,1);
y=avg_vel_EO(:,2);
z=avg_vel_EC(:,2);
% p1 = polyfit(x,y,4);
% p2 = polyfit(x,z,4);
% x2 = 15:0.5:30;
% y2 = polyval(p1,x2);
% z2 = polyval(p2,x2);
% plot(x,y,'*',x2,y2,'o')
% hold on
% plot(x,z,'o',x2,z2,'r')
plot(x,y,x,z)

title('EO Vs. EC: Average COP Velocity (A-P)');
xlabel('Distance(mm)');
ylabel('Avg. Velocity(A-P)(mm/s)');
legend('EO','EC');

COP_velocity_EO=[COP_velocity_avg(2,:);COP_velocity_avg(5,:);COP_velocity_avg(8,:);COP_velocity_avg(11,:)];
COP_velocity_EC=[COP_velocity_avg(32,:);COP_velocity_avg(35,:);COP_velocity_avg(38,:);COP_velocity_avg(41,:)];
figure(7)
% plot(distance(:,1),COP_velocity_trials(:,1),'o')
% hold on
% coef_fit=polyfit(distance(:,1),COP_velocity_trials(:,1),4);
% y_fit = polyval(coef_fit,[0 60]);
% plot([0 60],y_fit,'r')
x=distanceE(:,1);
y=COP_velocity_EO(:,1);
z=COP_velocity_EC(:,1);
% p1 = polyfit(x,y,4);
% p2 = polyfit(x,z,4);
% x2 = 15:0.5:30;
% y2 = polyval(p1,x2);
% z2 = polyval(p2,x2);
% plot(x,y,'*',x2,y2,'o')
% hold on
% plot(x,z,'o',x2,z2,'r')
plot(x,y,x,z)

title('EO Vs. EC: Average COP Velocity');
xlabel('Distance(mm)');
ylabel('Avg. Velocity(mm/s)');
legend('EO','EC');
