clear all; close all;

% Script to run synthetic tests with synthetic catchment
% Last adapted on 26/6/2020 by F. De Doncker

% Variables that can be changed:
% 1) Forward parameters:
%   - flags
%   - settings G matrix (number of tracers, number of geological units)
%   - number of subcatchment samples for additional data
%   - number of analysed ages (needed to construct tracer concentrations)
% 2) Inverse parameters
%   - L: smoothing distance
%   - sigma_cov: model variance
%   - forw_err: data uncertainty


%% Parameters and dimensions
%DIMENSIONS
load('synthetic_catchment_setting.mat'); %contains setting of synthetic
                                         %catchment, output of FastScape
                                         %model, X and Y  dimensions of
                                         %the catchment
x  = linspace(botX,topX,xNodes);
y  = linspace(botY,topY,yNodes);
dx = topX/(xNodes-1);
dy = topY/(yNodes-1);
[xgrid,ygrid] = meshgrid(x,y);
[nx,ny] = size(xgrid);

for j=1:ny
    for i=1:nx
        ij=i+(j-1)*nx;
        xv(ij)=xgrid(i,j);
        yv(ij)=ygrid(i,j);
        h(i,j)=hv(ij);
        catchx(i,j)=catchment(ij);
    end
end
[r,c] = find(catchx==1);
bound_x = [min(c) max(c)];
bound_y = [min(r) max(r)];

% FORWARD PARAMETERS
% 1. Flags
flag_river_samples= 1;  %if flag == 1: river samples;
                        %if flag ~= 1: extract geology for upstream
                        %catchment and append to A
                        % --> change n_datapoints to change the number of
                        % added samples
flag_blocks       = 1; %if flag == 1: true er. pattern = block pattern; else: gaussian bump
flag_checker      = 0; %if flag == 1: true erosion pattern = checkerboard pattern
connectivity_experiment = 0; %if ==1: only pixels with more than 5 donors totally exported to outlet (others: 50%)

% 2. Parameters for if true erosion rate map == gaussian bump
sigma_bump  = 10000;
xc          = ((max(c)+min(c))/2)*dx;
yc          = ((max(r)-min(r))/4)*dy;
amplitude   = 3.e-3;

% 3. Settings G matrix
ntracer      = 15; %number of tracers
igeol_unit   = 10; %number of geological units

% 4. Number of subcatchment samples
n_subc_datapoints = 2; %number of datapoints (catchment sediment samples)

% 5. Parameters of tracer concentrations
nage         = 1e2; %assume 100 ages are measured per geological unit

% INVERSE PARAMETERS
sigma_cov   = 1e-2;
L           = 2*dx;
forw_err     = 0.1; %error on data

%% Create geological map and tracer concentration for the entire domain
% Initiate G matrix
G = zeros(ntracer,nx*ny);

%Create catchment mask
mask_catchm = find(catchment ~= 0);

% Randomly select centre points of geological polygons
rng('shuffle');
sel_shuffle = mask_catchm(1,randperm(length(mask_catchm),igeol_unit));
x_vo = xv(sel_shuffle);
y_vo = yv(sel_shuffle);
tracer_conc_geol=rand(ntracer,igeol_unit);


% Pick a random age distribution for each geological unit
rng(10)
for igeol=1:igeol_unit
    mean_age=20.+abs(randn)*100;
    sigma_age=50+abs(randn)*100;
    [concentration]=generate_logNdistribution(mean_age,sigma_age,ntracer+1,nage);
    tracer_conc_geol(:,igeol)=concentration';
end

% Initialise geological map
geol_map = zeros(1,nx*ny);

% Fill G matrix
for itracer=1:ntracer
    for ineig=1:nx*ny
        if catchment(ineig) > 0
            %find closest neighbour
            disti=sqrt((10*dx)^2 + (10*dy)^2);
            jcell=1;
            for jneig=1:igeol_unit
                dist=sqrt((x_vo(jneig)-xv(ineig))^2+(y_vo(jneig)-yv(ineig))^2);
                if dist < disti
                    disti=dist;
                    jcell=jneig;
                end
            end
            if G(itracer,ineig) == 0
                G(itracer,ineig)=tracer_conc_geol(itracer,jcell);
                geol_map(1,ineig) = jcell;
            end
        end
    end
end
G = G(1:ntracer,:);

% Create G matrix only for catchment pixels
G_xy_Catch=zeros(ntracer,length(mask_catchm));
for itracer=1:ntracer
    index=0;
    for ij=1:ny*nx
        if catchment(ij) > 0 
            index=index+1;
            G_xy_Catch(itracer,index)=G(itracer,ij);
        end
    end
end

%% Create erosion rate map
% Extract erosion vector and concentration matrix within the catchment
% 1. Initialization
index=0;
index_map_catchm = zeros(1,nx*ny);

% 2. Create erosion map with block pattern
if flag_blocks == 1
    for ij=1:ny*nx
        if catchment(ij) > 0 
            index=index+1;
            ijk2index(ij)=index;
            xe(index)=xv(ij);
            ye(index)=yv(ij);
            if     yv(ij) <  topY/4
                eros=1.e-3;
            elseif yv(ij) >= topY/4 && yv(ij) < (topY*3/4)
                eros=3.e-3;
            elseif yv(ij) >= (topY*3/4)
                eros=0;
            end
            erosion(index,1)=eros;
            index_map_catchm(ij) = index; %create map containing indices of catchment
        end
    end
    
% 3. OR: create erosion map with gaussian pattern
else
    for ij=1:ny*nx
        if catchment(ij) > 0 %&& yv(ij)< (topY*3/4)
            index=index+1;
            ijk2index(ij)=index;
            xe(index)=xv(ij);
            ye(index)=yv(ij);
            exponent = ((xe(index)-xc).^2 + (ye(index)-yc).^2)./(2*sigma_bump^2);
            erosion(index,1) = amplitude  * exp(-exponent);
            index_map_catchm(ij) = index; %create map containing indices of catchment
        end
    end
end

% 4. OR: create erosion map with checkerboard pattern
if flag_checker == 1
    b = 15;
    block1 = ones(b,b).*1e-3;
    block2 = ones(b,b).*3e-3;
    block_large = [block1 block2; block2 block1];
    erosion = repmat(block_large,ceil(nx/(2*b)), ceil(ny/(2*b)));
    erosion(nx+1:end,:) = [];
    erosion(:,ny+1:end) = [];
    erosion = reshape(erosion, [nx*ny,1]);
    erosion(index_map_catchm==0) = [];
end



%% Connectivity experiment
% Create mask to assess impact of river transport: eroded material only
% to outlet when number of water donating pixels > 5 
% => for data multiplication: d = A*(e*mask)

if connectivity_experiment == 1
    n_donors_mask = n_donors;
    n_donors_mask(n_donors_mask>1) = 1;
    n_donors_mask(n_donors_mask<=1) = 0.7;
else
    n_donors_mask = ones(size(n_donors));
end
n_donors_mask_catchm = n_donors_mask;
n_donors_mask_catchm(index_map_catchm<1) = [];

%% Compute data
% Create data with error eps
d       = G_xy_Catch * (erosion .* n_donors_mask_catchm') + eps;
sigma_d = zeros(size(d));
sigma_d(:) = d.*forw_err;

% Choose prior erosion rate equal to the mean rate
edot_prior=ones(size(erosion))*harmmean(erosion);

%% Add subcatchment data
% This part relies on output of FastScape stored in 
% 'synthetic_catchment_setting.mat'. 
% For more information on 'stack' and 'n_donors' and other variables, 
% see DOI: 10.1016/j.geomorph.2012.10.008 ("A very efficient O(n), implicit 
% and parallel method to solve the stream power equation governing fluvial 
% incision and landscape evolution")

if flag_river_samples == 1
    
    %1. Distribute random points
    catchm_ind  = find(index_map_catchm > 0);
    w_don       = catchm_ind(n_donors(catchm_ind)>10);
    rng('shuffle')
    selection_random = datasample(w_don,n_subc_datapoints); %coordinates LEM
    selection_stack = stack(selection_random);
    
    % 2. Create data with a 10% error
    % 2.1. Initialization
    d_i             = zeros(ntracer,1);
    [stack2,I]      = sort(stack); %I represents the stack order for every cell
    stack_i         = zeros(1,n_subc_datapoints);
    n_don_i         = zeros(1,n_subc_datapoints);
    G_catch_old     = G_xy_Catch;
    
    % 2.2. Find coordinates of upstream catchment of every randomly selected
    %      datapoint and search catchments and subcatchments
    for i = 1:n_subc_datapoints
        ind     = selection_random(i);  %index in LEM map of selected datapoint (coordinate)
        
        n_don_i(i)  = n_donors(ind);    %total number of donors of selected point
        stack_i(i)  = I(ind);           %stack number for datapoint i
        
        catch_i     = stack(stack_i(i):stack_i(i)+n_don_i(i));    %LEM map coordinates of upstream catchment
        catch_index = index_map_catchm(catch_i);        %catchment coordinates
        
        row_i = length(d)+1;
        row_end = row_i + ntracer - 1;
        row_array = row_i : row_end;
        
        d(row_array,:)   = G_catch_old(:,catch_index) * ...
            (erosion(catch_index).*n_donors_mask(catch_i)')...
            + eps;
        
        G_xy_i = zeros(ntracer,length(erosion));
        G_xy_i(:,catch_index) = G_catch_old(:,catch_index);
        G_xy_Catch(row_array,:) = G_xy_i(:,:);
    end
    
    %3. Add 10% error on d
    sigma_d     = zeros(size(d));
    sigma_d     = d.*forw_err;
end


%% Inversion assuming a spatial covariance
% Building the model and data covariance matrices
for i=1:length(erosion)
    for j=i:length(erosion)
        dist(i,j)=sqrt((xe(i)-xe(j))^2+(ye(i)-ye(j))^2);
        cov(i,j)=sigma_cov^2*exp(-(dist(i,j)/L));
        cov(j,i)=cov(i,j);
        dist(j,i) = dist(i,j);
    end
end
cov_d=diag(sigma_d.^2);

% Solve the inverse problem
H       = cov*G_xy_Catch'*pinv((G_xy_Catch*cov)*G_xy_Catch'+cov_d);
edot_post = edot_prior+H*(d-G_xy_Catch*edot_prior);

% Estimate the posterior covariance
cov_post=cov-H*G_xy_Catch*cov;

% Estimate resolution
R=H*G_xy_Catch;

% Compute spread
spread = erosion*0;
for i=1:length(erosion)
    for j=1:length(erosion)
        spread(i)=spread(i)+(dist(i,j)*R(i,j).^2);
    end
end


% Posterior and prior variances
sigma_post   = sqrt(diag(cov_post).');   %posterior variance
sigma_prior  = sqrt(diag(cov).');        %prior variance

%% Plotting
%1. Catchment with samples

figure()
plot_streams = n_donors;
plot_streams2 = reshape(plot_streams,[nx, ny]);
plot_streams2(catchx==0) = NaN;
surface(xgrid,ygrid,plot_streams2); hold on
if flag_river_samples == 1
    scatter3(xgrid((selection_random)),ygrid((selection_random)),ones(size(selection_random)),'r')
end
cm = brewermap([max(max((plot_streams2)))],'Blues');
for i = 1:2
    cm(i,:) = [];
end
for i = 3:20
    cm(i,:) = cm(300,:);
end
for i = 20:50
    cm(i,:) = cm(500,:);
end
for i = 50:100
    cm(i,:) = cm(700,:);
end
for i = 100: max(max((plot_streams2)))
    cm(i,:) = cm(end,:);
end
colormap(cm);
title('samples')
axis equal; axis tight


%2. Prior erosion rate estimation
plot_prior = index_map_catchm;
plot_prior(plot_prior>0) = erosion;
plot_prior2 = reshape(plot_prior,[nx,ny]);
plot_prior2(catchx==0) = NaN;

figure()
surface(xgrid,ygrid,plot_prior2);
title('prior erosion rate')
colormap(flipud(brewermap([],'YlGnBu')))
axis equal; axis tight
caxis([0 3e-3])
colorbar


%3. Posterior erosion rate estimation
plot_post = index_map_catchm;
plot_post(plot_post>0) = edot_post;
figure()
plot_post2 = reshape(plot_post,[nx,ny]);
plot_post2(catchx==0) = NaN;
surface(xgrid,ygrid,plot_post2)
colormap(flipud(brewermap([],'YlGnBu')))
title('posterior erosion rate')
axis equal; axis tight
caxis([0 3e-3])
colorbar

%4. Residuals
figure()
d_m = G_xy_Catch * edot_post;
errorbar(d_m,d,sigma_d,'o'); hold on
xlabel('data (m^3/yr)')
ylabel('model data (m^3/yr)')
axis square
set(gca,'FontSize',12)

plot_geol = reshape(geol_map,[nx,ny]);
plot_geol(catchx==0) = NaN;
figure()
surface(xgrid,ygrid,plot_geol);
cm = brewermap([],'Set3');
cm = cm(1:igeol_unit,:);
colormap(cm)
title('geological map')
axis equal; axis tight
colorbar
unique(plot_geol(~isnan(plot_geol)))

figure()
diff_p_f = plot_post-plot_prior;
diff_p_f(catchx==0) = NaN;
surface(xgrid,ygrid,reshape(diff_p_f,[nx,ny]))
colormap(flipud(brewermap([],'RdBu')))
title('difference posterior - forward')
axis equal; axis tight
colorbar
caxis([-0.003 0.003])
disp(['Sum difference =', num2str(sum(abs(diff_p_f(~isnan(diff_p_f))))) ' m/y']);

plot_covp = index_map_catchm;
plot_covp(plot_covp>0) = sigma_post./sigma_prior;
plot_covp = reshape(plot_covp,[nx,ny]);
plot_covp(catchx==0) = NaN;
disp(['Normalised variance = ' mean(mean( (diag(cov_post).')/(diag(full(cov)).') ))])

plot_R = index_map_catchm;
plot_R(plot_R>0) = R(500,:)';
plot_R = reshape(plot_R,[nx,ny]);
plot_R(catchx==0) = NaN;

plot_covp2 = index_map_catchm;
plot_covp2(plot_covp2>0) = (diag(cov_post).');
plot_covp2 = reshape(plot_covp2,[nx,ny]);
plot_covp2(catchx==0) = NaN;

figure()
surface(xgrid,ygrid,plot_covp);
title('Posterior variance (normalized)')
colormap(brewermap([],'Reds'))
axis equal; axis tight
caxis([0.5 1])
colorbar

figure()
surface(xgrid,ygrid,plot_R);
colormap(flipud(brewermap([],'RdBu')))
title('Resolution centre point')
axis equal; axis tight
caxis([-max(max(plot_R)) max(max(plot_R))])
colorbar
hold on
scatter3(xgrid((index_map_catchm==500)),ygrid((index_map_catchm==500)),1,'MarkerEdgeColor','k',...
    'MarkerFaceColor',[0 .6 .75])

figure()
surface(xgrid,ygrid,plot_covp2);
colormap(brewermap([],'Reds'))
title('Posterior variance')
axis equal; axis tight
colorbar

figure()
oops=(erosion-edot_post)*1e3;
rhist(oops(:),10)
xlabel('e_{true}-e_{post} (mm/yr)')
axis square
set(gca,'FontSize',12)

[r,c] = size(dist);
plot_dist = reshape(dist,[1,r*c]);
plot_resolution = reshape(abs(R),[1,r*c]);

plot_dist_c = dist(500,:);
plot_resolution_c = abs(R(500,:));

[sort_dist,I] = sort(plot_dist);
sort_resolution = plot_resolution(I);

unique_dist = unique(sort_dist);
mean_resolution = unique_dist*0;
median_resolution = mean_resolution;

for i = 1:length(unique_dist)
    index_sort_dist = find(sort_dist==unique_dist(i));
    values_sort_resolution = sort_resolution(index_sort_dist);
    mean_resolution(i) = mean(values_sort_resolution);
    median_resolution(i) = median(values_sort_resolution);
end

[rown,coln] = size(R);

j = ones(1,rown);
r = 1:rown;
r2 = r.^2;

n = j*R*j';

sum_x = r*R*j';
sum_y = j*R*r';

sum_x2 = r2*R*j';
sum_y2 = j*R*r2';

sum_xy = r*R*r';

sample_corr_coef = ( (n * sum_xy) - (sum_x * sum_y) ) / ...
    ( sqrt( (n * sum_x2) - (sum_x).^2 ) * sqrt( (n * sum_y2) - (sum_y)^2 ) );


disp(['sample correlation coefficient = ', num2str(sample_corr_coef)])

figure();
plot(plot_dist,plot_resolution,'o','MarkerSize',2,'MarkerEdgeColor',[0.7 0.7 0.7])
hold on
plot(plot_dist_c,plot_resolution_c,'o','MarkerSize',4,'MarkerEdgeColor',[189 21 68]/255)
hold on
plot(unique_dist,median_resolution,'Color', [80, 103, 138]/255,'LineWidth',1)
hold on
plot(unique_dist,mean_resolution,'k','LineWidth',1)
hold on; lala = gca; lala_YLim = lala.YLim;
rectangle('Position',[0 min(lala_YLim) L max(lala_YLim)+abs(min(lala_YLim))],'FaceColor', [23 150 189 255*0.2]/255 , 'EdgeColor', [23 150 189 255*0.7]/255, 'LineWidth',0.75)
legend(['Resolution; sample correlation coefficient =' num2str(sample_corr_coef) ],'Resolution of catchment centre point','Median resolution', 'Mean resolution')
xlabel('Distance'); ylabel('Resolution')

plot_spread = index_map_catchm;
plot_spread(plot_spread>0) = spread;
plot_spread = reshape(plot_spread,[nx,ny]);
plot_spread(catchx==0) = NaN;

figure()
surface(xgrid,ygrid,plot_spread);
colormap(brewermap([],'Reds'))
title('Spread')
axis equal; axis tight
%caxis([0 3e-3])
colorbar

plot_relief = h;
[cmap,~] = demcmap(h);
plot_relief(catchx==0) = NaN;
figure(); surface(xgrid,ygrid,plot_relief);
colormap(cmap);
axis equal; axis tight;
title('Relief of the synthetic catchment')


