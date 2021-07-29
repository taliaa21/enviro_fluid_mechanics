% ***************************************************************************
%  Initialize all profiles and closure parameters 
%      Call once before time loop
%      Sets all forcing: pressure gradients, stresses, etc.
%      Should be used to adjust initial temperature/salinity profiles
%      Velocity initialized to zero
%      Turbulence quantities initialized to "SMALL"; Lengthscale parabolic
% ***************************************************************************

% ***************************************************************************
%  Physical parameters 
% ***************************************************************************
z0=0.01; %bottom roughness [m]
zb=10*z0; %bottom height
g=9.81; %m^2/s - gravity
C_D = 0.0025; %friction coefficient
SMALL=1e-6;
kappa=0.4; %von Karman constant
nu=1e-6; %m^2/s kinematic viscosity

%Flags for certain parameters/plots
steady = 0; %Set to 1 for steady (T_Px = 0)
unstratified = 1; %Set to 1 to have unstratified (delC = 0)
irradiance = 1; %Set to 1 to perform irradiance 
pass_scalar = 0; %Set to 1 for passive scalar case (alpha=0)
do_coeff = 1; %Set to 1 to have freshwater creep
print = 0; %Set to 1 to print parameters used and Si
plot1 = 0; %Set to 1 to plot Salinity
plot2 = 0; %Set to 1 to plot Velocity
plot3 = 0; %Set to 1 to plot Turbulent Kinetic Energy
plot4 = 0; %Set to 1 to plot deltaC
plot5 = 0; %Set to 1 to plot Nutrients
plot6 = 0; %Set to 1 to plot Phytoplankton
plot7 = 0; %Set to 1 to plot Zooplankton
plot8 = 0; %Set to 1 to plot Progression each depth
plot9 = 0; %Set to 1 to plot Diffusivity
plot10 = 0; %Set to 1 to plot Density
plot11 = 0; %Set to 1 to plot NPZ surfplots
plot12 = 0; %Set to 1 to plot Progression each time
plot13 = 1; %Set to 1 to plot NPZ subplots

% ***************************************************************************
%  Forcing parameters - Need to modify to allow for time variable Px.
% ***************************************************************************
%Loop to find min Px0 before Si is 1
% while Si < 1
%     Px0 = Px0/1.00001;
%     Si = g*abs(alpha)*Cx*H^2/(C_D*(1/(rho0*2*pi)*Px0*1000*3600*T_Px)^2);
% end

Px0 = 0.00015; %[kPa*m^2/1000kg=m/s^2] Magnitude on pressure gradient forcing. 1/rho0 already included.  0.00015 0.001 0.00010799

if steady == 1
    T_Px = 0; %Set to 0 for steady
else
    T_Px = 12.4; %Period (hours) on pressure gradient forcing. 
end

Cx = 0.0005; %[psu/m] Along channel salinity gradient 0.0005 0.0001 
Nutx = 0; %Set these to 0 for now
Phyx = 0;
Zx = 0;

%Buoyancy & Swimming
wb = 0; %[m/s] Phytoplankton buoyancy
ws = 0; %[m/s] Zooplankton swimming

%From NPZ Model
min_C=0.000001; %[psu] 
max_C=35; %[psu]
min_Nut=0.000001; %[ugN/L]
min_P=0.000001; %[ugN/L]
min_Z=0.000001; %[ugN/L]

N0=1.6 * 100; %[ugN L-1]
P0=0.3 * 100; %[ugN L-1] 
Z0=0.1 * 100; %[ugN L-1]
Vm=2.0/86400; %[s-1] Max P growth rate 
Ks=1.0; %[ugN L-1] Half-saturation constant
mp=0.1/86400; %[s-1] Mortality rate phytoplankton
gamma=0.3; %Messy eating 
Rm=1.5/86400; %[s-1] Max Z grazing
ivlev=1.0; %[ugN L-1] Ivlev constant
gz=0.2/86400; %[s-1] Z mortality rate
grazing=1; %Type of grazing to use. grazing=1 means ivlev, 0 means MP (Mayzaud–Poulet)

if irradiance == 1 %Irradiance term-set to ones(1,N) when don't want I./Iopt.*exp(1-I./Iopt)
    %Irradiance - parameters from Roy
    Vol = 1.6; %Average body volume phytoplankton
    Ec = 0.2; %Extinction coefficient water
    Ext = 0.12*Vol^-0.33; %Self shading of phytoplankton due to body size
    Iopt = 680; %Optimum surface solar radiation for photosynthesis
    Ir = 800; %Surface solar irradiance.

    %Roy irradiance effect
    I = Ir.*exp((-Ec-Ext)*abs(z)); 

    %Mark irradiance effect
    % L = 2; %Lengthscale FIXME
    %I = Ir.*exp(-abs(z)/L); 
    
    f_I = I./Iopt.*exp(1-I./Iopt); 
else
    f_I = ones(1,N);
end

% Whale parameters
eat = 1; %Flag for eating
diam = 7.5; %Diameter of whale [m]
eat_ts_1 = M/2; %Top eating time step bound [s]
eat_ts_2 = M/2 + 60; %Bottom eating time step bound [s] = top + 1 hour
eat_z_1 = -10; %Top eating depth bound [m]
eat_z_2 = eat_z_1 - diam; %Bottom eating depth bound [m]

poo = 1; %Flag for pooping
n_deposit = 4.03; % Nitrogen addition [ug/L per s]
Nut_add = n_deposit * dt; %Nutrient addition [ug/L]
poo_ts = M/2; %Time step for pooping (instantaneous)
poo_z_1 = -10; %Top pooping depth bound [m]
poo_z_2 = poo_z_1 - diam; %Bottom pooping depth bound [m]

% ***************************************************************************
%  Turbulence closure parameters 
% ***************************************************************************
A1=0.92;
A2=0.74;
B1=16.6;
B2=10.1;
C1=0.08;
E1=1.8;
E2= 1.33;
E3=0.25;
Sq=0.2;

% ********************************************************************
%  Setup initial conditions for scalar and density
% ********************************************************************

if unstratified == 1
    delC=0; %set to zero for Unstratified Case FIXME-thermocline?
    delN=0; %set to zero for Unstratified Case
    delP=0; %set to zero for Unstratified Case
    delZ=0; %set to zero for Unstratified Case
else
    delC=1; %change in salinity at initial themocline [psu]; 
    delN=1; %change in nutrients at initial themocline [ugN/L]; 
    delP=1; %change in phytoplankton at initial themocline [ugN/L]; 
    delZ=1; %change in zooplankton at initial themocline [ugN/L]; 
end

zdelC = -5; %position of initial thermocline
dzdelC = -2; %width of initial thermocline
zdelN = -5; %position of initial thermocline
dzdelN = -2; %width of initial thermocline
zdelP = -5; %position of initial thermocline
dzdelP = -2; %width of initial thermocline
zdelZ = -5; %position of initial thermocline
dzdelZ = -2; %width of initial thermocline

if pass_scalar == 1 %Set to 1 for passive scalar case (alpha=0)
    alpha = 0;
else
    alpha = -7.5*10^-4; %thermal expansivity or coefficient haline contraction [psu-1]
end

if do_coeff == 1 %Set to 1 to have freshwater creep
    coeff = 0.3; %ratio of constant freshwater creep
else
    coeff = 0; 
end
rho0=1000; %kg/m^3 - water density
T = 15; %psu scalar salinity. 15 default case

%Salinity set up
for i=1:N
    C(i)=T; %Initially uniform salinity profile
    if z(i)<=zdelC-0.5*dzdelC
        C(i)=T;
    elseif z(i)>=zdelC+0.5*dzdelC
        C(i)=T+delC;
    else
        C(i)=T+delC*(z(i)-zdelC+0.5*dzdelC)/dzdelC;
    end
   rho(i)=rho0*(1-alpha*(C(i)-T)); % Single scalar, linear equation of state
end

%Brunt-Vaisala frequency from density profile
N_BV(1)=sqrt((-g/rho0)*(rho(2)-rho(1))/(dz));
for i =2:N-1
   N_BV(i)=sqrt((-g/rho0)*(rho(i+1)-rho(i-1))/(2*dz));
end
N_BV(N)=sqrt((-g/rho0)*(rho(N)-rho(N-1))/(dz));

%Nutrient set up
for i=1:N
    Nut(i)=N0; %Initially uniform salinity profile
    if z(i)<=zdelN-0.5*dzdelN
        Nut(i)=N0;
    elseif z(i)>=zdelN+0.5*dzdelN
        Nut(i)=N0+delN;
    else
        Nut(i)=N0+delN*(z(i)-zdelN+0.5*dzdelN)/dzdelN;
    end
end

%Phytoplankton set up
for i=1:N
    P(i)=P0; %Initially uniform salinity profile
    if z(i)<=zdelP-0.5*dzdelP
        P(i)=P0;
    elseif z(i)>=zdelP+0.5*dzdelP
        P(i)=P0+delP;
    else
        P(i)=P0+delP*(z(i)-zdelP+0.5*dzdelP)/dzdelP;
    end
end

%Zooplankton set up
for i=1:N
    Z(i)=Z0; %Initially uniform salinity profile
    if z(i)<=zdelZ-0.5*dzdelZ
        Z(i)=Z0;
    elseif z(i)>=zdelZ+0.5*dzdelZ
        Z(i)=Z0+delZ;
    else
        Z(i)=Z0+delZ*(z(i)-zdelZ+0.5*dzdelZ)/dzdelZ;
    end
end

% *******************************************************************
%  Set initial conditions for u, q^2, and other turbulence quantities
% *******************************************************************

t(1)=0;
for i=1:N
   U(i) = 0.0;
   V(i) = 0.0;
   Q2(i)=SMALL;  %"seed" the turbulent field with small values, then let it evolve
   Q2L(i)=SMALL;
   Q(i)=sqrt(Q2(i));
   L(i)=-kappa*H*(z(i)/H)*(1-(z(i)/H)); %Q2L(n,1)/Q2(n,1) = 1 at initialization;
   Gh= -(N_BV(i)*L(i)/(Q(i)+SMALL))^2; 
   Gh=min(Gh,0.0233);
   Gh=max(Gh,-0.28);
   num=B1^(-1/3)-A1*A2*Gh*((B2-3*A2)*(1-6*A1/B1)-3*C1*(B2+6*A1));
   dem=(1-3*A2*Gh*(B2+6*A1))*(1-9*A1*A2*Gh);
   Sm(i)=num/dem;
   Sh(i)=A2*(1-6*A1/B1)/(1-3*A2*Gh*(B2+6*A1)); 
   Kq(i)=Sq*Q(i)*L(i)+nu;  %turbulent diffusivity for Q2
   nu_t(i)=Sm(i)*Q(i)*L(i)+nu; %turbulent viscosity
   Kz(i)=Sh(i)*Q(i)*L(i)+ nu; %turbulent scalar diffusivity
end

% *******************************************************************
%  Pre-define Tridiagonal Arrays - Just in case
% *******************************************************************
aU = zeros(1,N);
bU = zeros(1,N);
cU = zeros(1,N);
dU = zeros(1,N);
aC = zeros(1,N);
bC = zeros(1,N);
cC = zeros(1,N);
dC = zeros(1,N);
aNut = zeros(1,N);
bNut = zeros(1,N);
cNut = zeros(1,N);
dNut = zeros(1,N);
aP = zeros(1,N);
bP = zeros(1,N);
cP = zeros(1,N);
dP = zeros(1,N);
aZ = zeros(1,N);
bZ = zeros(1,N);
cZ = zeros(1,N);
dZ = zeros(1,N);
aQ2 = zeros(1,N);
bQ2 = zeros(1,N);
cQ2 = zeros(1,N);
dQ2 = zeros(1,N);
aQ2L = zeros(1,N);
bQ2L = zeros(1,N);
cQ2L = zeros(1,N);
dQ2L = zeros(1,N);

% *******************************************************************
%  Save initial conditions as first columns in saved matrix
% *******************************************************************
Um(:,1) = U;
Cm(:,1) = C;
Nutm(:,1) = Nut;
Pm(:,1) = P;
Zm(:,1) = Z;
Q2m(:,1) = Q2;
Q2Lm(:,1) = Q2L;
rhom(:,1) = rho;
Lm(:,1) = L;
nu_tm(:,1) = nu_t;
Kzm(:,1) = Kz;
Kqm(:,1) = Kq;
N_BVm(:,1) = N_BV;
   