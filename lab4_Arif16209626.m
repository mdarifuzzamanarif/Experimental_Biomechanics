%EMG : Lab 4
%Arif,Md Arifuzzaman, UMKC ID: 16209626
clear
clc
close all

fs=1000;
dt=1/fs;
fc=10;
%load raw data

% % name_list = {'trial1.txt''trial2.txt''trial3.txt''trial4.txt''trial5.txt''trial6.txt'};
% name_list = {'trial1.xlsx','trial2.xlsx'};
% for i = 1:length(name_list)
%     placehold=xlsread(name_list{i});
%     %A{i}=placehold;
% end

name_list = {'trial1.txt' 'trial2.txt' 'trial3.txt' 'trial4.txt' 'trial5.txt' 'trial6.txt'};

for i = 1:length(name_list)
    placehold=importdata(name_list{i});
    %placehold_show = placehold;
    temp{i}=placehold.data;
end
%column 1:32, trial 1:5
%t Fx Fy Fz Mx My Mz Fx Fy Fz Mx My Mz Fx Fy Fz Mx My Mz Fx Fy Fz Mx My Mz
%1 2  3  4  5  6  7  8  9  10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25
%EMG1 EMG2 EMG3 EMG4 AcX AcY AcZ
% 26   27   28   29   30  31  32

data_raw = temp{2};
data_filt = filter_data(data_raw,fc,fs,[2:32]);
acc = data_filt(:,30:32);
t = (1/fs)*(1:1:length(acc));

emg_raw = data_raw(:,26:29);
%plot raw EMG data
figure(1)
subplot(2,2,1); plot(data_raw(:,1),emg_raw(:,1));title ('Trial 2: Left TA');
subplot(2,2,3); plot(data_raw(:,1),emg_raw(:,3));title ('Trial 2: Right TA');

subplot(2,2,2); plot(data_raw(:,1),emg_raw(:,2));title ('Trial 2: Left Soleus');
subplot(2,2,4); plot(data_raw(:,1),emg_raw(:,4));title ('Trial 2: Right Soleus');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%[row,column]=size(emg1_raw)
for j=1:4
    emg_dtrend(:,j)=emg_raw(:,j)-(mean(emg_raw(:,j))); %DC Offset
end

emg_rect=abs(emg_dtrend(:,1:4)); %rectification
emg_env=filter_data(emg_rect,fc,fs,[1:4]); %filtered: EMG linear envelope

%plot filtered(envelope)EMG data
figure(2)
subplot(2,2,1); plot(data_raw(:,1),emg_env(:,1));title ('Trial 2: Left TA');
subplot(2,2,3); plot(data_raw(:,1),emg_env(:,3));title ('Trial 2: Right TA');

subplot(2,2,2); plot(data_raw(:,1),emg_env(:,2));title ('Trial 2: Left Soleus');
subplot(2,2,4); plot(data_raw(:,1),emg_env(:,4));title ('Trial 2: Right Soleus');

%accelerometer data onset
ionset_temp = find(abs(acc)>10*1000);
ionset = ionset_temp(1);
tonset = t(ionset);

%plot accelerometer data
figure(3)
subplot(3,1,1); plot(t,acc(:,1));title ('Trial 2: AccX');line([tonset tonset],ylim,'Color','r') 
subplot(3,1,2); plot(t,acc(:,2));title ('Trial 2: AccY');line([tonset tonset],ylim,'Color','r') 
subplot(3,1,3); plot(t,acc(:,3));title ('Trial 2: AccZ');line([tonset tonset],ylim,'Color','r') 

%EMG threshold
threshold=mean(emg_env)+ 3*(std(emg_env)); %threshold: DiFabio Method

%EMG data onset
%EMG1
for k = 1:4
ionset_EMG1 = find(emg_env(:,k)>threshold(k));
ionset_EMG1 = ionset_EMG1(1);
tonset_EMG1(k) = t(ionset_EMG1);
end

%plot filtered(envelope)EMG data
figure(4)
subplot(2,2,1); plot(t,emg_env(:,1));title ('Trial 2: Left TA');line([tonset_EMG1(1) tonset_EMG1(1)],ylim,'Color','r') 
subplot(2,2,3); plot(t,emg_env(:,3));title ('Trial 2: Right TA');line([tonset_EMG1(3) tonset_EMG1(3)],ylim,'Color','r') 

subplot(2,2,2); plot(t,emg_env(:,2));title ('Trial 2: Left Soleus');line([tonset_EMG1(2) tonset_EMG1(2)],ylim,'Color','r') 
subplot(2,2,4); plot(t,emg_env(:,4));title ('Trial 2: Right Soleus');line([tonset_EMG1(4) tonset_EMG1(4)],ylim,'Color','r') 

for m = 1:4
Perturb_EMG_RTA(m) = tonset_EMG1(m) - tonset;
RMS_emg(m) = rms(emg_env(m));

%APA
emg_i=emg_env(:,m);
APA(m) = rms(emg_i(ionset-300:ionset));
PPR(m) = rms(emg_i(ionset+50:ionset+150));

end





