%% Sea Level elastic shorelines
%Solve the sea level equation for elastic effects and include shoreline
%migration



addpath SLFunctions
%maximum degree to which spherical transformations are carried out
maxdeg = 64;



rho_ice = 916;
rho_water = 1000;
g = 9.81;

%%--------
% input
%---------


%make mesh grid for lat long

N = 512;
long_0 = linspace(0,360,2*N+1);
[x,w] = GaussQuad(N);
lat = acos(x)*180/pi - 90;
colat = lat + 90;
long = long_0(1:end-1);

[longitude, latitude] = meshgrid(long, lat);
LON = longitude;
LAT = latitude;

%%--------------
%ICE
%%------------------

%load ice data - gives height of ice over coordinates where there is ice
load ice5g_data;

%ice_end takes the matrix at the last time point (21 ya LGM)
ice_end = squeeze(ice5g(end,:,:));

% ice beg takes matrix at the first time point (now)
ice_beg = squeeze(ice5g(1,:,:));

%change in ice over time period
del_I_LGM = ice_end - ice_beg;

[ICELON, ICELAT] = meshgrid(ice_lon, ice_lat);

%interpolate onto meshgrid of lat long that we made
del_I = interp2(ICELON, ICELAT, del_I_LGM, LON, LAT);

%clean up nans
del_I(isnan(del_I))= 0;



%load topography present
load Topo_20;

%move half of matrix to last half
Z = [Z(:, length(Z)/2+1: length(Z)) Z(:, 1:length(Z)/2)];


%topo lat is divided by number of items in first row
topolat = linspace(-90,90, length(Z(:,1)));
topolong = linspace(0,360, length(Z(1,:)));


%latitute is column, longitude is row, want to put these on our mesh grid
%we made

%Specify T_0 and C_0
topo_0 = interp2(topolong, topolat, Z, LON, LAT);
%make ocean function
c_0 = zeros(size(topo_0));
c_0(topo_0<0) = 1;

%% ------
% Love numbers
%% ------

load LN_l70_ump5_lm5
T_lm = get_tlm(maxdeg);
h_lm = love_lm(h_el, maxdeg);
k_lm = love_lm(k_el, maxdeg);
E_lm = 1 + k_lm - h_lm;

%%------
% start solving sea level equation
%%--------



%setup topo/ocean function 
c_j = c_0;
topo_j = topo_0; %+iceincrement+sedincrement; %if need to correct for sed or ice


%expand ice function into spherical coordinates
delI_lm = spa2sph(del_I, maxdeg, long, colat);



% max number of iterations to be performed
imax = 10;
%initialize conv
conv = 'not converged yet';

%initialize epsilon criteria
epsilon = 0.01;

for i = 1:imax
    
    %expand topo into spherical
    t_j_lm = spa2sph(topo_j, maxdeg, long, colat);
    %expand ocean into spherical
    c_j_lm = spa2sph(c_j, maxdeg, long, colat);
    %calculate topography correction - 1 where ocean is gained
    topo_corr_j = topo_0.*(c_j-c_0);
    TO_lm = spa2sph(topo_corr_j,maxdeg, long,colat);
    
    
    if i ==1
        %set initial guess for delS_lm
        delS_lm = (c_j_lm./c_j_lm(1)).*((-rho_ice/rho_water).*delI_lm(1)+TO_lm(1))-TO_lm;
    end
    
    %solve curly SL
    delL_lm = rho_water*delS_lm + rho_ice*delI_lm;
    delCurlySL_lm = T_lm.*E_lm.*delL_lm; %PUT IN LOVE NUMBERS
    delCurlySL = sph2spa_tp1(delCurlySL_lm, maxdeg, long, colat);
    
    disp('after SL')
    %term for curly SL times ocean fucntion
    delCurlySL_timesOc = delCurlySL.*c_j;
    delCurlySL_timesOc_lm = spa2sph(delCurlySL_timesOc, maxdeg, long, colat);
   
    %solve delPhi/g
    delPhig = 1/c_j_lm(1)*((-rho_ice/rho_water)*delI_lm(1)...
        - delCurlySL_timesOc_lm(1)+TO_lm(1));
    
    %update delta SL
    delSL = delPhig + delCurlySL;
    
    
    %update tophography and ocean function
    topo_j = topo_0 - delSL;% is it plus or minus delSL???
    c_j = zeros(size(topo_j));
    c_j(topo_j<0) =1;
    %or
    %          c_j = zeros(size(topo_j));
    %          [cr,cc]=find(topo_j < 0);
    %          c_j(cr,cc)=1;
    
    %calculate new delta S
    delS_new = delSL.*c_j- topo_0.*(c_j - c_0);
    delS_new_lm = spa2sph(delS_new, maxdeg, long, colat);
    
    disp('got delS!!')
    %check convergence
    chi = abs((sum(abs(delS_new_lm)) - sum(abs(delS_lm)))./sum(abs(delS_lm)));
    
    if chi < epsilon
        conv = 'converged';
        disp(['Converged with chi of ' num2str(chi) 'after ' num2str(i) ' iterations'])
        break
        
    else
        conv = 'not converged yet';
        disp(['Not converged yet! Chi is' num2str(chi) 'after' num2str(i) 'iterations'])
    end
    %update delS ocean height
    delS_lm = delS_new_lm;
    
    
end

%plot results

figure
pcolor(LON, LAT, delSL)
hold on

%plot topography contour
[xc, yc] = getcon(LON, LAT, topo_0, 0);

plot(xc,yc,'-k')

shading flat
axis image
colorbar

