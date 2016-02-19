%Solve the sea level equation for elastic and viscous effects and include shoreline
%migration

addpath SLFunctions
addpath ReadSaveLN
addpath Love_numbers
%ICE load
load ice5g_griddata;
%load ice5g_gl; %sam's ice grids and lat long and time
%load LN_l120_ump5_lm5
%load VM2_tp
load l71Cump3lm8
%load LN_l90_VM2 %this is different VM2!!
%only replace these love numbers with sam's format
%load VM2_sam
% load 120p55
% h_el = he;
% k_el = ke;
% h_amp = h;
% k_amp = k;
% spoles = s;


%maximum degree to which spherical transformations are carried out
maxdeg = 256;
rho_ice = 920;
rho_water = 1000;
g = 9.80665;
tic;
%%--------
% input
%---------


%make mesh grid for lat long
%
  N =256;
long_0 = linspace(0,360,2*N+1);
[x,w] = GaussQuad(N);
lat = acos(x)*180/pi - 90;
colat = lat + 90;
long = long_0(1:end-1);
% LON = longitude;
% LAT = latitude;
% 
% long = longitude(1:end-1);%use this for sam's lat long
% lat = latitude;
% colat = lat' + 90;
 [LON, LAT] = meshgrid(long, lat);



timefilename = 'zxx';
fid = fopen(timefilename);
time_increment = fscanf(fid,'%f');
fclose(fid);
timeintervals = length(time_increment);
time = flipud(time_increment);
%time = times; %use this for sam's time
%%--------------
%ICE FOR ELASTIC - ice5g from LGM to pres
%%------------------
%load topography present
%load Topo_20;
%load toposam
topofilename = 'topo_3.i.txt';
fid = fopen(topofilename);
%find time
initialtime = fscanf(fid, '%f', 1);

%make topography matrix
r = 256;
c = 512;
topo_i = zeros(r,c);
for n = 1:r
    %read in data by row
    topo_i(n,:) = transpose(fscanf(fid, '%f', c));
end

%close file
fclose(fid);

% %move half of matrix to last half
% Z = [Z(:, length(Z)/2+1: length(Z)) Z(:, 1:length(Z)/2)];
%
% %topo lat is divided by number of items in first row
% topolat = linspace(-90,90, length(Z(:,1)));
% topolong = linspace(0,360, length(Z(1,:)));

%latitute is column, longitude is row, want to put these on our mesh grid
%we made

%% Specify T_0 and C_0
%Topo_0 = interp2(topolong, topolat, Z, LON, LAT);
Topo_0 = topo_i; %topobr
%Topo_0 = topobr;
%make ocean function
C_0 = zeros(size(Topo_0));
C_0(Topo_0<0) = 1;
%% ------
% Love numbers
%% ------

%load LN_l70_ump5_lm5
T_lm = get_tlm(maxdeg);
h_lm = love_lm(h_el, maxdeg);
k_lm = love_lm(k_el, maxdeg);
E_lm = 1 + k_lm - h_lm;

%%------
%% Start solving sea level equation
%%--------

%% Compute ice load increment
%expand ice function into spherical coordinates
%INSERT LOADS
% %set ice load
iceload = ice5g_grid; %these are in opposite order 
%iceload = ice5g_gl;

% %making delta load increments
% %for ice increments
delI_inc_lm = cell(size(time));%spherical ice increments
delI_inc = cell(size(time)); %spatial grids of ice increments
I = cell(size(time));
delI_lm = cell(size(time)); %total ice in spherical - used in elastic

%% Initiate cells to store topo, c, delS, DELS, SL

topo_j_cell = cell(length(time),1);
c_j_cell = cell(length(time),1);
delS_lm_save =  cell(length(time),1);
DelS_lm_save =  cell(length(time)+1,1);

SLmat =  cell(length(time),1);
kmax = 3;
imax = 10;
%initialize epsilon criteria
epsilon = 1E-4;
for k = 1:kmax
    
    if k ==1 %Initialize first loop to be initial topography
        for ti = 1:length(time)
            topo_j_cell{ti} = Topo_0;
            c_j_cell{ti} = C_0;
        end
    end
    %% Create beta using ice grids and check for floating ice
    
    
    Istar = cell(size(time));
    beta =  cell(length(time),1);
    for t = 1:length(time)
        
        Istar{t} = squeeze(iceload(end-t+1,:,:)); %USE THIS WITH MY
        %ICE5G GRIDS
        %Istar{t} = squeeze(iceload(:,1:end-1,t)); %use these with SAMS ICE5G
        
%         %for no floating ice check comment this out and fill beta's with
%         %ones
%         groundice = squeeze(iceload(end-t+1,:,:)) - (rho_water/rho_ice)*abs(topo_j_cell{t});
%         %groundice = squeeze(iceload(:,1:end-1,t)) - (rho_water/rho_ice)*abs(topo_j_cell{t});
%         groundice(c_j_cell{t} ==0) = 0;%only look at marine ice - don't look at continents
%         Istar{t}(groundice<0) = 0; %floating ice!
%         beta{t} = zeros(size(groundice));
%         beta{t}(Istar{t}==0)= 1; %ones where no ice, 0 where there is grounded ice
%         
        beta{t} = ones(size(Istar{1})); %uncomment this to get rid of
        %floating ice check!! 
        beta{t}(Istar{t} > 0) = 0;
        
        I{t} = spa2sph(Istar{t},maxdeg,long,colat);
        delI_lm{t} = I{t};
        if t == 1
            delI_inc_lm{t} = zeros(size(delI_lm{1}));
        else
            delI_inc{t} = squeeze(Istar{t} - Istar{t-1});
            delI_inc_lm{t} = spa2sph(delI_inc{t},maxdeg,long,colat);
        end
    end
    
    viscous = zeros(length(time),length(delI_lm{1}));% save viscous component at each time
    DelS_lm_save{1} = zeros(size(delI_lm{1}));
    
    for j =1:length(time)
        
        topo_0 = topo_j_cell{1};
        topo_j = topo_j_cell{j};
        c_0 = c_j_cell{1};
        c_j = c_j_cell{j};
        c_jstar = c_j.*beta{j}; % c_j is c_jstar, includes grounded ice
        c_0star = c_0.*beta{1};
        %expand ocean into spherical
        c_jstar_lm = spa2sph(c_jstar, maxdeg, long, colat);
        
        %calculate topography correction, - 1 where ocean is gained
        TO_lm = spa2sph(topo_0.*(c_jstar-c_0star),maxdeg, long,colat);
        
        
        for i =1:imax
            %% Set initial guess for dS for i ==1
            if i ==1 && k ==1 %if first iteration and first topo iteration
                if j ==1 %if j = 1 make prev beta same as beta
                    prevbeta = beta{j};
                else
                    prevbeta = beta{j-1};
                end
                Clm = spa2sph(c_0.*prevbeta,maxdeg,long,colat);
                Tcurl_lm = spa2sph(topo_0.*c_0.*(beta{j} - prevbeta),maxdeg,long,colat);
                %set initial guess for delS_lm
                delS_lm = (Clm./Clm(1)).*((-rho_ice/rho_water).*delI_inc_lm{j}(1)+Tcurl_lm(1))-Tcurl_lm;
                
            else if i ==1 % if k>1 than use from previous iteration
                    delS_lm = delS_lm_save{j};% if k ! = 1, then use delS from previous topography round
                end
            end
            %% Define loads
            delL_lm = rho_ice*(delI_lm{j}-delI_lm{1}) + ... %Total elastic load
                rho_water*(DelS_lm_save{j}+delS_lm);% j for DelS is actually
            %j-1, for j==1 this is 0
            delS_lm_save{j} = delS_lm;
            if j ==1
                viscous(j,:) = zeros(1,length(delI_lm{1}));
            else
                viscous(j,:) = visc(N,maxdeg,time,j,delI_inc_lm,delS_lm_save,spoles,h_amp,k_amp);
            end
            %% Calculate curly SL
            delCurlySL = sph2spa_tp1(T_lm.*(squeeze(viscous(j,:))+ E_lm.*delL_lm),maxdeg,long,colat);
            
            %% term for curly SL times ocean function
            ROlm = spa2sph(delCurlySL.*c_jstar,maxdeg,long,colat);
            
            %% Compute delPhi/g
            delPhig = 1/c_jstar_lm(1)*((-rho_ice/rho_water)*(delI_lm{j}(1)-delI_lm{1}(1))...
                - ROlm(1)+TO_lm(1));
            
            delSL = delPhig + delCurlySL;
            %delSLCstar_lm = spa2sph((delPhig + delCurlySL).*c_jstar,maxdeg,long,colat);
            %DelS_lm_save{j+1} = delSLCstar_lm - TO_lm;
            DelS_lm_save{j+1} = spa2sph((delSL).*c_jstar - topo_0.*(c_jstar-c_0star),maxdeg,long,colat);
            %alternatively
            %DelS_lm_save{j+1} = ROlm + phiC - TO_lm;
            
            delS_lm_new = DelS_lm_save{j+1} - DelS_lm_save{j};
            SLmat{j} = delSL;
            
            %% Check convergence
            chi = abs((sum(abs(delS_lm_new)) - sum(abs(delS_lm)))./sum(abs(delS_lm)));
            if j ==1
                chi = 0;
            end
            if chi < epsilon
                conv = 'converged';
                %disp(['Converged with chi of ' num2str(chi) 'after ' num2str(i) ' iterations'])
                break
            else
                conv = 'not converged yet';
                %disp(['Not converged yet! Chi is' num2str(chi) 'after' num2str(i) 'iterations'])
            end
            %update delS ocean height
            delS_lm = delS_lm_new;
            delS_lm_save{j} = delS_lm;
            j
            
            
        end
    end
    
    chi_t = abs((sum(abs(spa2sph(Topo_0,maxdeg,long,colat))) - sum(abs(spa2sph(topo_j_cell{end-1},maxdeg,long,colat))))./sum(abs(spa2sph(topo_j_cell{end-1},maxdeg,long,colat))))
    %     if chi_t < epsilon
    %         conv = 'converged';
    %         disp(['Converged with chi_t of ' num2str(chi_t) 'after ' num2str(k) ' iterations'])
    %         %break
    %     else
    %         conv = 'not converged yet';
    %         disp(['Not converged yet! Chi_t is' num2str(chi_t) 'after' num2str(k) 'iterations'])
    %     end
    %% Shift Shorelines
    %update tophography and ocean function
    topo_j_cell{1} = Topo_0+SLmat{end};
    for n =1:length(time)
        topo_j_cell{n} = topo_j_cell{1} - SLmat{n};% is it plus or minus delSL???
        c_j_cell{n} = zeros(size(topo_j_cell{n}));
        c_j_cell{n}(topo_j_cell{n}<0) =1;
    end
    
    k
end
toc;