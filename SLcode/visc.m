% Solving for viscous and elastic part of SL code - calculates and returns
% curly SL
% USE 1:time is last time step to present for delI_inc, i.e delI_inc(1) is
% last time step
function savevisc = visc(N,maxdeg,time,timestep,delI_inc_lm,delS_lm_save,spoles,h_amp,k_amp)
addpath SLFunctions
addpath ReadSaveLN
%load LN_l90_VM2
%load LN_l120_ump5_lm5
%load LN_l120_ump5lm5_s
% load 120p55
% h_amp = h;
% k_amp =k; 

long_0 = linspace(0,360,2*N+1);
[x,w] = GaussQuad(N);
lat = acos(x)*180/pi - 90;
colat = lat + 90;
long = long_0(1:end-1);

[longitude, latitude] = meshgrid(long, lat);
LON = longitude;
LAT = latitude;
rho_ice = 920;
rho_water = 1000;

%save values for each time step
savevisc = zeros(1,length(delI_inc_lm{1}));

spoles_1 = spoles;
spoles_1(spoles==0) = 1;% use spoles with 1's for 0's for division to avoid nans


%initialize time steps
ti = timestep;

%begin time loop
for t = 1:ti-1
    dL_lm = rho_ice*delI_inc_lm{t}+ rho_water*delS_lm_save{t};%for now JUST ICE
    tj = time(ti); %current time step
    tn = time(t); % %time steps from -inf to prev time step
    %viscous component
    exp_part = exp(-spoles*(abs(tj - tn)));
    hovers = (k_amp-h_amp)./spoles_1;
    sumoverk = sum(hovers.*(1-exp_part),2);
    %add l = 0 part
    sumoverk_l0 = [0, sumoverk'];
    
    %multiply by load
    dLRvisc = zeros(size(dL_lm));
    for l = 0:maxdeg
        pos_alm = l*(l+1)/2;%find positions in vector of sph coordinates
        dLRvisc(pos_alm+1:pos_alm+l+1) = dL_lm(pos_alm+1:pos_alm+l+1)*sumoverk_l0(l+1);
    end
    % sum over time steps
    savevisc = savevisc + dLRvisc;
end
end



