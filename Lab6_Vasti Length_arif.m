%lab 6 : Muscle Modelling
%Arif, Md Arifuzzaman; ID: 16209626

clear
clc
close all

v_length = csvread('vasti_length.csv',1);
moment_arm = 3e-2;
tv = v_length(:,1);
L_p = v_length(:,2);

%muscle.m
%program to simulate muscle contractions using 3-component Hill model


%muscle properties for human Vastus Lateralis

Lslack = 0.223; % slack length of SEE
Umax = 0.04;    % strain in SEE is 4% of Fmax
Lceopt = 0.093; % optimal length of CE
width = 0.63*Lceopt; % maximum length change of CE
Fmax = 7400;    %maximal isometric force
a = 0.25*Fmax;  % force-velocity parameter a
b = 0.25*10*Lceopt; % force-velocity parameter b (Nigg & Herzog, p. 174-175)

%set initial condition for state variable Lce and initialize ODE solver
Lce = 0.087; %this makes SEE just slack at t=0
t = tv(1); tend = tv(end); % start time and duration of experiment

%h = 0.001; % integration step size in seconds
%h=0.01;
h = (tv(end)-tv(1))/length(tv);
i = 1; % step counter
data = zeros(tend/h,2); %space to store time and force results 

%start simulation
while (t < tv(end))
    disp('hit')
    %prescribed ramp shortening profile
    
%     if (t<=1) Lm = 0.31; end %initial muscle+tendon length
%     if (t>1 & t<2) Lm = 0.31 - 0.04*(t-1); end % ramp shortening at 4 cm/s
    
    %calculate force in SEE from current SEE length
    Lm =L_p(i);
    Lsee = Lm - Lce; % length of SEE is total length minus CE length
    if (Lsee < Lslack)
        F = 0;  % SEE is slack, so no force
    else
        F = Fmax*((Lsee-Lslack)/(Umax*Lslack))^2; % SEE has quadratic force-length relationship
    end
    
    % calculate isometric force at this Lce from force-length relationship of CE
    Fo = max([0 Fmax*(1-((Lce-Lceopt)/width)^2)])
    % claculate Lcedot from Hill's equation assuming that muscle is
    % maximally activated
    
%     if (F > Fo)
%         disp('Error: program cannot do ecentric contractions'); return; end
    Lcedot = -b*(Fo-F)/(F+a); %note: Vce is negative for shortening!
    
    if Lcedot>0
        F = Fmax*(1.2*(Lsee-Lslack)/(Umax*Lslack))
    end
    
    % do one forward Euler integration step to calculate new Lc and store
    % data
    
    Lce = Lce + h*Lcedot;
    %t = t+h; 
    t=tv(i);
    i=i+1;
    data(i,:) = [t F];
    
end
data(1,:)=[];
figure(1)
plot (data(:,1),data(:,2));
xlabel('Time [s]');ylabel('Force [N]');
Title('Force predicted from muscle model')

figure(2)
plot (v_length(:,1),v_length(:,2));
xlabel('Time [s]');ylabel('Vasti length [N]');
Title('Vasti length change')
        
       
%knee moment from muscle model
kneeM_muscle = data(:,2).*moment_arm;
figure(3)
plot (v_length(:,1),kneeM_muscle(:,1));
xlabel('Time [s]');ylabel('Knee moment [N]');
Title('Knee moment predicted from muscle model')

load('kneeM_ID.mat');
fs=100;
t = (1/fs)*(1:1:length(M_leg_proxG));
figure(4)
plot (t,M_leg_proxG(:,1),t,M_leg_proxG(:,2),t,M_leg_proxG(:,3));
xlabel('Time [s]');ylabel('Knee moment [N]');
Title('Knee moment calculated from Inverse dynamics')
Legend('A-P','M-L','vertical');
   


