%Frequency Analysis :Biceps Fatigue
%Arif, Md. Arifuzzaman
clear
clc
close all
fs = 1000;

% load data
name_list = {'trial1.txt' 'trial2.txt' 'trial3.txt' 'trial4.txt' 'trial5.txt' 'trial6.txt'};

for i = 1:length(name_list)
    placehold=importdata(name_list{i});
    %placehold_show = placehold;
    temp{i}=placehold.data;
end

data = temp{6}; %Biceps Fatigue


for i =1:2000:118000
    
emg_bi=data(i:i+2000,2);

% emg_bi=data(:,2);
%power spectrum analysis using power_spectrum.m
[f,p] = power_spectrum(emg_bi,fs);
p(1) = [];
f(1) = [];
Area = trapz(p)/2;
cum_Area = cumtrapz(p);
cond_mf = find(cum_Area>Area);
freq_ref = cond_mf(1);
m_f = f(freq_ref)
MedFreq(i,1)= m_f;
end
MedFreq(find(MedFreq==0)) = []; %only median frequencies as as sliding window of 2000 frames

%plot
figure(1)
plot(f,p) 
hold on
line ([m_f,m_f],ylim,'Color','r');
title('Biceps fatigue')
xlabel('Frequency(Hz)')
ylabel('Power(W)')
legend('Power Spectrum','Median Power Frequancy')
hold off

figure(2) %Median power frequnecy 
time_ref = [1:2000:118000];
plot(time_ref,MedFreq,'*')
Title('Median Frequency after every 2 seconds as muscle fatigues')
ylabel('Med Freq')
xlabel('Time Reference')

% figure(3)
% time_ref=time_ref';
% x=time_ref(:,1);
% y=MedFreq(:,1);
% p = polyfit(x,y,4);
% x2 = 0:0.20:12;
% y2 = polyval(p,x2);
% plot(x,y,'o',x2,y2)



