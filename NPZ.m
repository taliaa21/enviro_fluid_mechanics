%NPZ Model from Franks et al. 1986 
%Parameters from http://www.soest.hawaii.edu/oceanography/courses/OCN681/hw3.html#parameters
%Summer 2019
%Lily Engel

%% Simplified NPZ - no light dependence
clear all
close all

dt=1e-3; %[days]
tstart=0; %[days]
tend=50; %[days]
min_N=0.000001; %[ugN/L]
min_P=0.000001; %[ugN/L]
min_Z=0.000001; %[ugN/L]

N0=1.6; %[ugN L-1]
P0=0.3; %[ugN L-1] 
Z0=0.1; %[ugN L-1]
Vm=2.0; %[days-1] Max P growth rate
Ks=1.0; %[ugN L-1] Half-saturation constant
m=0.1; %[days-1] Mortality rate phytoplankton
gamma=0.3; %Messy eating 
Rm=1.5; %[days-1] Max Z grazing
ivlev=1.0; %[ugN L-1] Ivlev constant
g=0.2; %[days-1] Z mortality rate
grazing=1; %Type of grazing to use. grazing=1 means ivlev, 0 means MP

t=tstart:dt:tend;
Nt=length(t);
N=zeros(Nt,1);
P=zeros(Nt,1);
Z=zeros(Nt,1);

N(1) = N0;
P(1) = P0;
Z(1) = Z0;

alpha=1.0;

for n=2:Nt
    if grazing == 0  %MP condition
        alpha=ivlev*P(n-1);
    end
    
    P(n) = P(n-1)+dt*((Vm*N(n-1)*P(n-1))/(Ks+N(n-1))-m*P(n-1)-alpha*Rm*(1-exp(-ivlev*P(n-1)))*Z(n-1));
    Z(n) = Z(n-1)+dt*((1-gamma)*alpha*Rm*(1-exp(-ivlev*P(n-1)))*Z(n-1)-g*Z(n-1));
    N(n) = N(n-1)+dt*((-Vm*N(n-1)*P(n-1))/(Ks+N(n-1))+m*P(n-1)+g*Z(n-1)+gamma*alpha*Rm*(1-exp(-ivlev*P(n-1)))*Z(n-1));
    
    %Boundary conditions
    P(n) = max(P(n),min_P);
    Z(n) = max(Z(n),min_Z);
    N(n) = max(N(n),min_N);
end

first = t(1);
last = t(length(t)); %If want to plot whole time, do first:last

tplotstart = find(t==first); %First hour want to plot
tplotend = find(t==last); %Last hour want to plot 21.67 

figure(1)
plot(t(tplotstart:tplotend),N(tplotstart:tplotend),t(tplotstart:tplotend),P(tplotstart:tplotend),t(tplotstart:tplotend),Z(tplotstart:tplotend))
title('Basic NPZ')
ylabel('Value (ugN/L)')
xlabel('Time (days)')
legend('N','P','Z')

%% NPZ Model from Franks et al. 1986 merged with Roy 2016
%Parameters from http://www.soest.hawaii.edu/oceanography/courses/OCN681/hw3.html#parameters
%Summer 2019
%Lily Engel

clear all
% close all

dt=1e-3; %[days]
tstart=0; %[days]
tend=200; %[days]
min_N=0.000001; %[ugN/L]
min_P=0.000001; %[ugN/L]
min_Z=0.000001; %[ugN/L]

N0=1.6; %[ugN L-1]
P0=0.3; %[ugN L-1] 
Z0=0.1; %[ugN L-1]
Vm=2.0; %[days-1] Max P growth rate
Ks=1.0; %[ugN L-1] Half-saturation constant
m=0.1; %[days-1] Mortality rate phytoplankton
gamma=0.3; %Messy eating 
Rm=1.5; %[days-1] Max Z grazing
ivlev=1.0; %[ugN L-1] Ivlev constant
g=0.2; %[days-1] Z mortality rate
grazing=1; %Type of grazing to use. grazing=1 means ivlev, 0 means MP

t=tstart:dt:tend;
Nt=length(t);
N=zeros(Nt,1);
P=zeros(Nt,1);
Z=zeros(Nt,1);

N(1) = N0;
P(1) = P0;
Z(1) = Z0;

%Add irradiance term from Roy
V = 1.6; %Average body volume phytoplankton
Ec = 0.2; %Extinction coefficient water
d = 4; %Depth of water
Iopt = 680; %Optimum surface solar radiation for photosynthesis
Ir = 800; %Surface solar irradiance.
Pmax = 2.62; %Max growth rate phytoplankton
Temp = 28.1; %[degC] Water temperature.
alpha1 = 6.05; %Coefficient
gamma1 = 0.09; %Coefficient

Ext = 0.12*V^-0.33; %Self shading of phytoplankton due to body size
I = Ir.*exp((-Ec-Ext)*d); %Irradiance effect
Leff = I./Iopt.*exp(1-I./Iopt); %Light effect on photosynthesis
delta = Pmax*exp(gamma1.*Temp).*Leff./alpha1; %Temperature dependent phytoplankton growth rate

alpha=1.0;
% delta = 1; %Comment out when doing irradiance

for n=2:Nt
    if grazing == 0  %MP condition
        alpha=ivlev*P(n-1);
    end
    
    P(n) = P(n-1)+dt*(Vm*delta*N(n-1)/(Ks+N(n-1))*P(n-1)-m*P(n-1)-alpha*Rm*(1-exp(-ivlev*P(n-1)))*Z(n-1));
    Z(n) = Z(n-1)+dt*((1-gamma)*alpha*Rm*(1-exp(-ivlev*P(n-1)))*Z(n-1)-g*Z(n-1));
    N(n) = N(n-1)+dt*(-Vm*delta*N(n-1)/(Ks+N(n-1))*P(n-1)+m*P(n-1)+g*Z(n-1)+gamma*alpha*Rm*(1-exp(-ivlev*P(n-1)))*Z(n-1));
    
    %Boundary conditions
    P(n) = max(P(n),min_P);
    Z(n) = max(Z(n),min_Z);
    N(n) = max(N(n),min_N);
end  

figure(2)
plot(t,N,t,P,t,Z)
title('Basic NPZ with Irradiance')
ylabel('Value')
xlabel('Time')
legend('N','P','Z')

