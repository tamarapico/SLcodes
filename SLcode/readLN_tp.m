% Read in love numbers from the code maxwell.m output from JXM
% J. Austermann 2013
% in file make sure to replace non-zero numbers in jerry's output txt file
% replace -306 with E-300

addpath Love_numbers
% import file
love = importdata('Love_numbers/prem.l71C.ump3.lm8',' ',1);
maxdegree = 256;

% initiate data vector
data_sc = zeros(2*length(love.data),1);

% put data in one colum vector
i = 1;
for j = 1:length(love.data)
    data_sc(i) =  love.data(j,1);
    data_sc(i+1) =  love.data(j,2);
    i = i+2;
end


% intiate variables

deg = zeros(maxdegree,1);
mode_found = zeros(maxdegree,1);
spoles1 = cell(maxdegree,1);
h_el = zeros(maxdegree,1);
l_el = zeros(maxdegree,1);
k_el = zeros(maxdegree,1);
h_fl = zeros(maxdegree,1);
l_fl = zeros(maxdegree,1);
k_fl = zeros(maxdegree,1);
h_amp1 = cell(maxdegree,1);
l_amp1 = cell(maxdegree,1);
k_amp1 = cell(maxdegree,1);

h_el_tide = zeros(maxdegree,1);
l_el_tide = zeros(maxdegree,1);
k_el_tide = zeros(maxdegree,1);
h_fl_tide = zeros(maxdegree,1);
l_fl_tide = zeros(maxdegree,1);
k_fl_tide = zeros(maxdegree,1);
h_amp_tide1 = cell(maxdegree,1);
l_amp_tide1 = cell(maxdegree,1);
k_amp_tide1 = cell(maxdegree,1);

data_sc(isnan(data_sc)) = [];
ind = 1;
for i = 1:maxdegree
    % read in degree and number of modes found
    deg(i) = data_sc(ind);
    ind = ind+1;
    mode_found(i) = data_sc(ind);
    ind = ind+1;
    if i == maxdegree
        degree_vec = data_sc(ind:end);
    else if i ==1
            degree_vec = data_sc(ind:ind +10+4*mode_found(i)-1);
        ind = ind+ 10+4*mode_found(i);
       else 
        degree_vec = data_sc(ind:ind +20+6*mode_found(i)+ mode_found(i)-1);
        ind = ind+ 20+6*mode_found(i)+ mode_found(i);
        end
    end
  

         j = mode_found(i);
        if i ==1
            max_mode = j;  
        end
        if j > max_mode
            max_mode = j;
        end
    %     ind = 2;
    
    % decay constants s
    vec_ind = 1;
    spoles1{i} = degree_vec(vec_ind:j);

    vec_ind = vec_ind+j;
    
    % elastic LN
    h_el(i) = degree_vec(vec_ind);
    l_el(i) = degree_vec(vec_ind+1);
    k_el(i) = degree_vec(vec_ind+4)/deg(i);
    
    % fluid LN
    h_fl(i) = degree_vec(vec_ind+5);
    l_fl(i) = degree_vec(vec_ind+6);
    k_fl(i) = degree_vec(vec_ind+9)/deg(i);
    vec_ind = vec_ind+9;
    
    % amplitudes of viscous LN
    h_amp1{i} = degree_vec(vec_ind+1:vec_ind+j);
    vec_ind = vec_ind+j;
    l_amp1{i} = degree_vec(vec_ind+1:vec_ind+j);
    vec_ind = vec_ind+j;
    k_amp1{i} = degree_vec(vec_ind+1:vec_ind+j)/deg(i);
    
    if i ==1
        % tidal LN only start at degree 2
    else
        
        vec_ind = vec_ind+j;
        % elastic tidal LN
        h_el_tide(i) = degree_vec(vec_ind+1);
        l_el_tide(i) = degree_vec(vec_ind+2);
        k_el_tide(i) = degree_vec(vec_ind+5)/deg(i);
        
        % fluid tidal LN
        h_fl_tide(i) = degree_vec(vec_ind+6);
        l_fl_tide(i) = degree_vec(vec_ind+7);
        k_fl_tide(i) = degree_vec(vec_ind+10)/deg(i);
        
        vec_ind = vec_ind+10;
        
        % amplitudes of viscous tidal LN
        h_amp_tide1{i} = degree_vec(vec_ind+1:vec_ind+j);
        vec_ind = vec_ind+j;
        l_amp_tide1{i} = degree_vec(vec_ind+1:vec_ind+j);
        vec_ind = vec_ind+j;
        k_amp_tide1{i} = degree_vec(vec_ind+1:vec_ind+j)/deg(i);
        
    end
end

spoles = zeros(maxdegree,max_mode);
h_amp = zeros(maxdegree,max_mode);
l_amp = zeros(maxdegree,max_mode);
k_amp = zeros(maxdegree,max_mode);
h_amp_tide = zeros(maxdegree,max_mode);
l_amp_tide = zeros(maxdegree,max_mode);
k_amp_tide = zeros(maxdegree,max_mode);

for n =1:maxdegree
    spoles(n,1:length(spoles1{n})) = spoles1{n};
    h_amp(n,1:length(spoles1{n})) = h_amp1{n};
    l_amp(n,1:length(spoles1{n})) = l_amp1{n};
    k_amp(n,1:length(spoles1{n})) =k_amp1{n};
    if n ==1 %tide ones are empty for degree ==1
    else
    h_amp_tide(n,1:length(spoles1{n})) = h_amp_tide1{n};
    l_amp_tide(n,1:length(spoles1{n})) = l_amp_tide1{n};
    k_amp_tide(n,1:length(spoles1{n})) = k_amp_tide1{n};
    end
end



save('prem.l71C.ump3.lm8.mat','spoles','mode_found','h_el','l_el','k_el','h_fl','l_fl','k_fl', ...
    'h_amp','l_amp','k_amp',...
    'h_el_tide','l_el_tide','k_el_tide', 'h_fl_tide','k_fl_tide','l_fl_tide', ...
    'h_amp_tide','k_amp_tide','l_amp_tide');


% fid = fopen('LoveNum.txt', 'w');
% fprintf(fid,'%f ', h_el, k_el, h_fl, k_fl);
% fclose(fid);
