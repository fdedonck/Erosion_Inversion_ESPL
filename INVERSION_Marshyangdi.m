clear all; close all;

% The model adapted to the Marsyandi Catchment w/ mineral concentration data
% Last adapted on 24/6/2020

% LINEAR INVERSE MODEL
% G is the tracer concentration matrix
% d the measured tracer distribution at the outlet of the catchment
% edot is the erosion rate (expected mean)
% in a forward case G*edot=d

% Here, the source areas are subcatchments of the Marsyandi catchment

%Input:
% 1) Source area map                        --->(source_areas)
%   (geological units or tributaries)
%   (geotiff file, projected coordinate system <-- specify this)
%   ('nanval' (e.g. zero) outside of source areas (nanval can be changed))

% 2) Sample data                            --->(SampleData, source, d)
%   (.csv file, first column: sample ID, other columns: different tracers,
%   every row corresponds to tracer information for a different source area)

% 3) Shapefile with sample locations        --->(source_samples)
%   (attribute table: ID corresponds to source area e.g. sample 3 is
%   representative for source area 3)
%   (contains n_datasamples samples for which the data will be inverted,
%   default = 1)

% 4) ID of sample w/ data to be inverted    --->(detrital_sample_tag)

% 5) Q_s                                    --->(Q_s)
%   (sediment discharge for catchment outlet, m^3/y)

%% Flags
flag_geol = 1;  %which type of source area is used? If it's geological data,
                %change this flag to 1, for tributary source areas, leave
                %the flag equal to 0

flag_w_sc = 0;  %with subcatchment samples == 1; note that subcatchment
                %samples can only be used with the dry season data (the
                %only additional subcatchment data have been sampled in
                %the dry season and there is no such dataset available for
                %the monsoon season

%% Data variable names: change to fit your data
% there are two possibilities of source areas:
% 1) geological source areas
% 2) tributary source areas
% when one of these options is chosen, variable input for the other option
% is not required. Be sure to adapt the flags to your choice.

% filename of geological source area map
tifname_geol =  'geology2_epsg32644.tif';

% filename of tributary source area map
tifname_trib = 'M_sc_epsg32644.tif';

% filename of subcatchments for additional data (leave empty if no
% subcatchment data is available)
% more downstream subcatchments should have a larger ID
tifname_subc = 'SubCatchmentsSamples.tif';

% NaN value for pixels outside of source areas
nanval = 0;

% projection system of .tif files
ProjectedCSTypeGeoKey   = 32644;

% filename source areas
%   this file is structured as follows:
%   each row represents a new source area
%   the first column contains the ID's of the different source areas, which
%   correspond to the ID's of the tif file of the source area map
%   columns 2 to end contain tracer concentration data
%   ! this file should also contain the detrital data with the
%   corresponding detrital sample ID in the first column (in the examples,
%   the detrital sample data corresponds to the last line of the .csv
%   files)
%       -> filename geological source areas
geol_data_name = 'Geol_Data_monsoon.csv';
%       -> filename tributary source areas
trib_data_name = 'SampleData_dryseason.csv';


%  additional subcatchment data (for geological source area setting)
subc_data_geol_name = 'SubCatchmentSamples.csv';

%  additional subcatchment data (for tributary source area setting)
subc_data_trib_name = 'SubCatchmentSamples_dry_catchm.csv';

% yearly sediment load
Q_s_yr = 16 * 10^6; %m^3/y

% scaling factor: which fraction of the entire catchment is covered by
% tributary catchments?
scaling_Qs = 0.7911;

% ID in .csv file (geol_data_name or trib_data_name) of detrital sample
%   for geological source areas
detrital_sample_ID_geol  = 9;
%   for tributary source areas
detrital_sample_ID_trib  = 11;

% shapefile name of sampling points for tributary source data; ! this 
% shapefile should also contain the location of the detrital sample whose
% data will be inverted
source_samples_points_name = 'samples_epsg32644.shp';

%% Inverse parameters
%prior erosion rate estimate
edot_prior_set  = 'mean';   % standard: equal to mean (derived from Qs)
                            % you can set this to any constant value that
                            % will be taken as prior erosion rate for the
                            % whole prior map. 
%model covariance
sigma_cov   = 3;        %standarddeviation model covariance (std:3)

%smoothing distance (expressed in pixelsize units)
L_set          = 5;     %critical distance covariance matrix (std: 5*dx)

% data uncertainty
d_sigma     = 0.1;

%% No need to change anything past this line

%% The inverse model
%1. Loading the catchment parameters
%1.1. Load source area map
if flag_geol == 1
    tifname =  tifname_geol;
else
    tifname = tifname_trib;
end
[source_areas, Rsource] = geotiffread(tifname);
info      = geotiffinfo(tifname);

if flag_w_sc == 1
    [subc, Rsubc] = geotiffread(tifname_subc);
end

source_areas      = double(source_areas);
mask = source_areas;
mask(mask~=nanval) = 1;
source_areas(source_areas==nanval) = NaN;


%1.2. Load data
%c1 : sample number, c2:cend : tracer data
if flag_geol == 1
    SampleData = table2array(readtable(geol_data_name));
else
    SampleData = table2array(readtable(trib_data_name));
end

SampleData(isnan(SampleData)) = 0;

if flag_geol == 1
    SubcData = table2array(readtable(subc_data_geol_name));
else
    SubcData = table2array(readtable(subc_data_trib_name));
end

if flag_geol == 1
    Q_s = Q_s_yr; %m^3/y
else
    Q_s = Q_s_yr * scaling_Qs; %m^3/y
end

%1.3. Link subcatchment index to sample index
n_datasamples       = 1;

if flag_geol == 1
    detrital_sample_ID  = detrital_sample_ID_geol;
else
    detrital_sample_ID  = detrital_sample_ID_trib;
end

if flag_geol == 1
    source_data = SampleData;
    ind = find(SampleData(:,1) == detrital_sample_ID);
    source_data(ind,:) = [];
    S = source_data(:,1)';
else
    source_samples  = shaperead(source_samples_points_name);
    n_samples       = length(source_samples) - n_datasamples;
    S               = zeros(n_samples,1);
    S_ind           = zeros(n_samples,1);
    
    j = 0;
    for i = 1 : n_samples + n_datasamples
        if source_samples(i).id ~= detrital_sample_ID
            j = j+1;
            S(j,1)      = source_samples(i).id;
            S_ind(j,1)  = find(SampleData(:,1) == S(i,1));
        end
    end
    
end
clear source samples

%1.4. Source data
if flag_geol == 1
    source = source_data(:,2:end);
else
    source = SampleData(S_ind,2:end);
end

n_tracers = size(source,2);

%1.5. Detrital data
% d is a vector containing tracer concentrations found at outlet
%(rows:concentration of tracer i)
detrital_sample_num = find(SampleData(:,1)==detrital_sample_ID);
d = SampleData(detrital_sample_num,2:end)';
d = d .* (Q_s * ...
    (1/(Rsource.CellExtentInWorldX*Rsource.CellExtentInWorldY))) * 1000;

if flag_w_sc == 1
    area_catchm = sum(sum(mask));
    area_pixel = Rsource.CellExtentInWorldX*Rsource.CellExtentInWorldY;
    subc_nums = unique(subc);
    subc_nums(subc_nums == 0) = [];
    n_subc = length(subc_nums);
    Q_s_0 = Q_s;
    for i = 1 : n_subc
        catchment_i = subc_nums(i);
        catchm{catchment_i} = [catchment_i];   %number in catchment map
        area_sc = sum(sum(mask(subc <= catchment_i)))*area_pixel;
        Q_s_i = Q_s_0 * (area_sc/area_catchm) ;
        d(end+1:end+n_tracers)  = SubcData(i,2:end)*Q_s_i*(1/area_pixel);
    end
end

% 2. Parameters of the model
%DIMENSIONS
[ny,nx] = size(source_areas); %nx= number of cells in x-direction (#columns)
nn = nx*ny;     %number of gridcells of map
dx = Rsource.CellExtentInWorldX;      %resolution in x direction (column)
dy = Rsource.CellExtentInWorldY;
Ly = dy*(ny-1);
Lx = dx*(nx-1);


edot_v  = zeros(nn,1);  %erosion model to be filled in (vectorial form)

index   = 0;            %needed to calculate grid spacing
mask_v  = reshape(mask,[nn,1]);
mask_size = sum(sum(mask));

%INVERSE PARAMETERS
%prior erosion rate estimate
if edot_prior_set == 'mean'
    edot_prior  = (Q_s / (mask_size*dx*dy) ) * 1000; %equal to mean erosion rate
else 
    edot_prior = edot_prior_set;
end

%smoothing distance
L          = L_set*dx;     %critical distance covariance matrix (std: 5*dx)


%GRID SPACING
index=0;
for i = 1:nx
    for j = 1:ny
        index    = index+1;
        x(index) = dx*(j-1);
        y(index) = dy*(i-1);
    end
end
x(mask_v==0) = [];
y(mask_v==0) = [];


%3. Create the tracer matrix
G_m = zeros(ny,nx,n_tracers);   %empty tracer matrix to be filled in (stacked maps)
G_v = zeros(n_tracers, mask_size);     % "" (stacked vectors)
for i=1:n_tracers           %tracer
    G_mi = source_areas;
    for j = 1:length(S)         %source area
        sample_idx = S(j);           %sample
        G_mi(G_mi==sample_idx)=source(j,i); %fill this source area with tracer(i) data
    end
    G_m(:,:,i) = G_mi;  %update tracer map matrix
    G_v_i  = reshape(G_mi,[nn,1]); %update vector tracer matrix
    G_v_i(mask_v == 0)  = []; %update vector tracer matrix
    G_v(i,:) = G_v_i;
end
G = G_v; %tracer matrix A is matrix with stacked vectors
clear A_m A_mi A_v A_vi

if flag_w_sc == 1
    G1 = G;
    subc_small = reshape(subc, [nx*ny,1]);
    subc_small(mask_v==0) = [];
    for i = 1:n_subc
        %select coordinates subcatchment
        catchment_i = subc_nums(i);
        catchm_i = catchm{catchment_i};
        [r,~] = size(G);
        for j = 1:max(size(catchm_i))
            coord = find(subc_small == catchm_i(j));
            %fill forward operator only for subcatchment (the rest remains 0)
            G(r+1:r+n_tracers,coord) = G1(:,coord);
        end
    end
end


% 4. Inverse model
%4.1. Prior erosion rate equal to mean erosion rate forward model
edot_prior=ones(nn,1)*edot_prior;
edot_prior(mask_v==0) = [];


%4.2. Covariance matrix calculation
%4.2.1. Model covariance
dist = zeros(mask_size,mask_size);
cov = zeros(mask_size,mask_size);
for i = 1:mask_size
    for j = i:mask_size
        dist(i,j)   = sqrt((x(i)-x(j))^2+(y(i)-y(j))^2);
        cov(i,j)    = sigma_cov^2*exp(-(dist(i,j)/L)^2);
        cov(j,i)    = cov(i,j); %model covariance
    end
end
disp('Model covariance computed')

%4.2.2. Data covariance
sigma_d = d.*d_sigma;
cov_d = diag(sigma_d.^2); %data covariance

%4.3. Solve the problem
%4.3.1. Inversion
H       = cov*G'*pinv((G*cov)*G'+cov_d);
edot_post = edot_prior+H*(d-G*edot_prior);


%4.3.2. Handle negative values
neg_post = abs( sum (edot_post(edot_post<0) ) ); %total negative erosion
total = Q_s * (1/(dx*dy)) * 1000;
perc = neg_post / total;

edot_post(edot_post<0) = 0;
edot_post(edot_post>0) = edot_post(edot_post>0) - ...
    ( perc * edot_post(edot_post>0) );

disp('edot post computed')

%4.4. Posterior uncertainties
C_M     = cov - H*G*cov;                %posterior covariance
R       = H*G;                          %resolution
res_dot = edot_post*0; res_dot(20) = 1; %show for which dot resolution was computed
spread  = edot_post*0;

for i = 1:length(edot_post)
    for j = 1:length(edot_post)
        spread(i) = spread(i)+(dist(i,j)*R(i,j).^2);
    end
end
disp('spread computed')

sigma_post   = sqrt(diag(C_M));         %posterior variance
sigma_prior  = sqrt(diag(cov));         %prior variance

mask_v(mask_v==0) = NaN;

sigma_post_v = mask_v; sigma_post_v(mask_v>0) = sigma_post;
sigma_xy     = reshape((sigma_post_v),[ny,nx]); %_xy denotes map form instead of vectorial form
sigma_prior_v = mask_v; sigma_prior_v(mask_v>0) = sigma_prior;
red_sigma_xy = reshape((sigma_post_v./sigma_prior_v),[ny,nx]);
clear sigma_post_v sigma_prior_v sigma_post sigma_prior

res20_v = mask_v; res20_v(mask_v>0) = R(20,:);
res_xy       = reshape(res20_v,[ny,nx]);
clear res20_v R

res_dot_v = mask_v; res_dot_v(mask_v>0) = res_dot;
res_dot      = reshape(res_dot_v,[ny,nx]);
clear res_dot_v

spread_v = mask_v; spread_v(mask_v>0) = spread;
spread_xy    = reshape(spread_v,[ny,nx]);
clear spread_v spread

edot_post_v = mask_v; edot_post_v(mask_v>0) = edot_post;
edot_esti = reshape(edot_post_v,[ny,nx]);
clear edot_post_v edot_post

%5. Plotting

figure()
if flag_geol == 1
    mapshow(edot_esti,Rsource,'DisplayType','surface');
else
    mapshow(zeros(size(source_areas)),Rsource,'CData',edot_esti,'DisplayType','surface');
    hold on; mapshow(source_samples)
end
hcb = colorbar; title('\bf Result $$\bf e_{inverse model}$$','Interpreter','latex')
colorTitleHandle = get(hcb,'Title');
titleString = 'Erosion rate (mm/y)';
set(colorTitleHandle ,'String',titleString);

figure()
mapshow(sigma_xy,Rsource,'DisplayType','surface');
colorbar; title('\bf Posterior variance','Interpreter','latex');

figure()
mapshow(red_sigma_xy,Rsource,'DisplayType','surface');
axis square
colorbar; title('\bf Reduced variance','Interpreter','latex');

figure()
mapshow(res_xy,Rsource,'DisplayType','surface');
axis square
colorbar; title('\bf Resolution','Interpreter','latex');

figure()
mapshow(res_dot,Rsource,'DisplayType','surface');
axis square
colorbar; title('\bf Resolution dot','Interpreter','latex');

figure()
mapshow(spread_xy,Rsource,'DisplayType','surface');
axis square
colorbar; title('\bf Spread','Interpreter','latex');

figure(); hist(edot_esti(:));
title('$$\bf e_{post}$$','Interpreter','latex');


%% Maps
%6. Exporting maps
key = struct( ...
    'GTModelTypeGeoKey',1, ...
    'GTRasterTypeGeoKey',1, ...
    'ProjectedCSTypeGeoKey',[ProjectedCSTypeGeoKey]);
geotiffwrite('e_dot_post_catchm_monsoon.tif',edot_esti,Rsource,'GeoKeyDirectoryTag',key);
geotiffwrite('red_sigma_catchm_monsoon.tif',red_sigma_xy,Rsource,'GeoKeyDirectoryTag',key);
geotiffwrite('spread_catchm_monsoon.tif',spread_xy,Rsource,'GeoKeyDirectoryTag',key);