%********************************************************************
% Master code for 200B water column (vertical) cases
%    Original code from Lisa Lucas, modified by Tina Chow
%       Spring 2018: Mark Stacey and Michaella Chung
%********************************************************************

clear all;
close all;

%********************************************************************
%Define model set up - grid and timestep
%********************************************************************
N=80;%number of grid points
H=50; %depth (meters)
dz=H/N; %grid spacing - may need to adjust to reduce oscillations
dt=60; %(seconds) size of time step 
M=3600; %number of time steps 3600, 1440, 2400, 31205 (*2)
beta=dt/dz^2;
for i=1:N % Initialize grid   
   z(i)=-H+dz*(i-1/2); %bottom at z=-H, free surface at 0
end
isave=10; %increments for saving profiles. set to 1 to save all; 10 saves every 10th, etc. 
savecount=0;

% *******************************************************************
%  Call code to set up initial conditions and model parameters
% *******************************************************************
wc_setup

% *******************************************************************
%  Start of time loop
% *******************************************************************
Cplot = [];
times = [];

for m=2:M
   t(m)=dt*(m-1); %define time
   wc_advance %uses BGO/Mellor-Yamada 2-equation closure

% *******************************************************************
%  Saving Profiles - isave defines decimation
% ******************************************************************

   if mod(m,isave) == 1
      savecount = savecount+1;
      Um(:,savecount) = U;
      Cm(:,savecount) = C;
      Nutm(:,savecount) = Nut;
      Pm(:,savecount) = P;
      Zm(:,savecount) = Z;
      Q2m(:,savecount) = Q2;
      Q2Lm(:,savecount) = Q2L;
      rhom(:,savecount) = rho;
      Lm(:,savecount) = L;
      nu_tm(:,savecount) = nu_t;
      Kzm(:,savecount) = Kz;
      Kqm(:,savecount) = Kq;
      N_BVm(:,savecount) = N_BV;
      Cplot(savecount) = Cm(N,savecount)-Cm(1,savecount);
      times(savecount) = t(m);
   end
end

% *******************************************************************
%  End of time loop
% *******************************************************************

% *******************************************************************
%  PLOTTING FOLLOWS HERE
%   columns of variable matrices (Um, Cm, etc) vs. z array
% *******************************************************************

%Salinity Plots
%H=20, T = 15, Px0 = 0.00015, delC=0, add salinity baroclinic to pressure,
%alpha = -0.00075, T_Px = 12.4 

if print == 1 %Print all the following 
    Si = g*abs(alpha)*Cx*H^2/(C_D*(1/(2*pi)*Px0*3600*T_Px)^2); %Simpson Number new way - preferred
    Si_opt2 = g*abs(alpha)*Cx*H^2/(C_D*max(max(abs(Um)))^2); %Simpson number Monismith et. al (PS 4 w/C_D - includes baroclinic term already in Um)

    Un = 2*pi*H/(T_Px*3600*sqrt(C_D)*(1/(2*pi)*Px0*3600*T_Px)); %Unsteadiness Number
    Ro = ws/(kappa*sqrt(C_D)*(1/(2*pi)*Px0*3600*T_Px)); %Rouse Number

    fprintf(['Alpha = ', num2str(alpha),'\n'])
    fprintf(['\nT_Px = ', num2str(T_Px),'\n'])
    fprintf(['\nPx = ', num2str(Px0),'\n'])
    fprintf(['\nCx = ', num2str(Cx),'\n'])
    fprintf(['\nSi = ', num2str(Si),'\n'])
    %fprintf(['\nSi_opt2 = ', num2str(Si_opt2),'\n'])
    fprintf(['\nUn = ', num2str(Un),'\n'])
    %fprintf(['\nRo = ', num2str(Ro),'\n'])
end

first = times(1)/3600;
last = times(length(times))/3600; %If want to plot whole time, do first:last

tstart = find(times/3600==first); %First hour want to plot
tend = find(times/3600==last); %Last hour want to plot

if plot1 == 1
    figure(1)
    contourf(times(tstart:tend)/3600,z,Cm(:,tstart:tend),'Linestyle','none');
    ylabel('Depth (m)')
    xlabel('Time (hrs)')
    title(['Salinity (psu) Cx = ',num2str(Cx),', coeff= ',num2str(coeff)])
    colorbar
end

if plot2 == 1
    figure(2)
    contourf(times(tstart:tend)/3600,z,Um(:,tstart:tend),'Linestyle','none');
    ylabel('Depth (m)')
    xlabel('Time (hrs)')
    title('Velocity (m/s)')
    colorbar
end

if plot3 == 1
    figure(3)
    contourf(times(tstart:tend)/3600,z,Q2m(:,tstart:tend));
    ylabel('Depth (m)')
    xlabel('Time (hrs)')
    title('Turbulent Kinetic Energy (m^2/s^2)')
    colorbar
end

if plot4 == 1
    figure(4)  %plot delta C at those two spots
    plot(times(tstart:tend)/3600,Cplot(:,tstart:tend),'o-');
    ylabel('Salinity (psu)')
    xlabel('Time (hrs)')
    title('Scalar Concentration Difference')
end

if plot5 == 1
    figure(5)
    contourf(times(tstart:tend)/3600,z,Nutm(:,tstart:tend));
    ylabel('Depth (m)')
    xlabel('Time (hrs)')
    title(['Nutrients (ugN/L), Nutx = ',num2str(Nutx)])
    colorbar
end

if plot6 == 1 
    figure(6)
    contourf(times(tstart:tend)/3600,z,Pm(:,tstart:tend));
    ylabel('Depth (m)')
    xlabel('Time (hrs)')
    title(['Phytoplankton (ugN/L), Phyx = ',num2str(Phyx),', wb= ',num2str(wb)])
    colorbar
end

if plot7 == 1
    figure(7)
    contourf(times(tstart:tend)/3600,z,Zm(:,tstart:tend));
    ylabel('Depth (m)')
    xlabel('Time (hrs)')
    title(['Zooplankton (ugN/L), Zx = ',num2str(Zx),', ws= ',num2str(ws)])
    colorbar
end

if plot8 == 1
  %  for i=N:-3:1 %See the progression each depth
        i = N;  %comment out if doing full loop
        figure(8)
        plot(times(tstart:tend)/3600,Nutm(i,tstart:tend),times(tstart:tend)/3600,Pm(i,tstart:tend),times(tstart:tend)/3600,Zm(i,tstart:tend));
        legend('N','P','Z')
        title(['NPZ at z = ',num2str(z(i)),' m'])
        xlabel('Time (hrs)')
        ylabel('Concentration (ugN/L)')
 %   end
end

if plot9 == 1
    figure(9)
    contourf(times(tstart:tend)/3600,z,Kzm(:,tstart:tend));
    ylabel('Depth (m)')
    xlabel('Time (hrs)')
    title('Diffusivity (m^2/s)')
    colorbar
end

if plot10 == 1
    figure(10)
    contourf(times(tstart:tend)/3600,z,rhom(:,tstart:tend));
    ylabel('Depth (m)')
    xlabel('Time (hrs)')
    title('Density (kg/m^3)')
    colorbar
end

if plot11 == 1
    figure(11)
    surf(times/3600,z,Nutm)
    hold on
    surf(times/3600,z,Pm)
    surf(times/3600,z,Zm)
    shading interp
    legend('N','P','Z')
    hold off
    colorbar
    xlabel('Time (hrs)')
    ylabel('Depth (m)')
    zlabel('Concentration (ugN/L)')
    title('NPZ')
end

if plot12 == 1
   % for i=tstart:5:tend %See the progression each time
        i = tend;  %comment out if doing full loop
        figure(12)
        plot(Nutm(:,i),z,Pm(:,i),z,Zm(:,i),z);
        legend('N','P','Z')
        title(['NPZ at t = ',num2str(times(i)/3600),' hrs'])
        xlabel('Concentration (ugN/L)')
        ylabel('Depth (m)')
  %  end
end

if plot13 == 1
    figure(13)
    subplot(3,1,1)
    contourf(times(tstart:tend)/3600,z,Nutm(:,tstart:tend));
    ylabel('Depth (m)')
    xlabel('Time (hrs)')
    title(['Nutrients (ugN/L), Nutx = ',num2str(Nutx)])
    colorbar
    caxis([130 200])
    subplot(3,1,2)
    contourf(times(tstart:tend)/3600,z,Pm(:,tstart:tend));
    ylabel('Depth (m)')
    xlabel('Time (hrs)')
    title(['Phytoplankton (ugN/L), Phyx = ',num2str(Phyx),', wb= ',num2str(wb)])
    colorbar
    caxis([0 55])
    subplot(3,1,3)
    contourf(times(tstart:tend)/3600,z,Zm(:,tstart:tend));
    ylabel('Depth (m)')
    xlabel('Time (hrs)')
    title(['Zooplankton (ugN/L), Zx = ',num2str(Zx),', ws= ',num2str(ws)])
    colorbar
    caxis([0 40])
end

% *******************************************************************
%  End of Main Program
% *******************************************************************