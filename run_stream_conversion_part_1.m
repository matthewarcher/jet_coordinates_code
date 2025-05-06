% RUN STREAM CONVERSION PART (1)
% Find core of Western Boundary Current using maximum velocity values, interactively
% pick region of good data, find direction of core at chosen locations.

clear;clc

% Pathways
cd /home/nfs/z3519950/sci-maths-ocean/Data/RADAR/VECTORS
addpath /home/nfs/z3519950/hdrive/hf_radar/Meandering

% Data
load CoffsHarbour_2012_2016_Prepped
load core_data_4yr_2

%% Pick a Test Time

% ts = find(t == datenum(2013,02,1,3,0,0)) % Choose start time
% te = find(t == datenum(2013,03,5,3,0,0)) % Choose end time

%% Run Code

% Initialize variables   
core = nan(length(lat),length(t(ts:te)));
dir = nan(length(lat),length(t(ts:te)));
core_orig = nan(length(lat),length(t(ts:te)));
n = 0;

for k = ts:te%1:length(t(ts:te))
    
    n = n + 1;
    
    % Timer
    display([num2str(n) ' / ' num2str(length(t(ts:te)))])
    
    % Initiate Variables at Time Step
    up = u(:,:,k);
    vp = v(:,:,k);
    
    % Run Function
    [coret,dirt,core_origt] = stream_conversion_p1(LON,LAT,up,vp);
    
    core(:,n) = coret;
    core_orig(:,n) = core_origt;
    dir(:,n) = dirt;
        
end

%% Manually Check&Change

figure('position',[100 60 900 900])
n=0;
for k = 12054:length(t)%te %2993   %172 is where I started?!
    
    n = k;%n + 1;
    
    if sum(~isnan(core(:,n)))>0
    
    % Initiate Variables at Time Step
    up = u(:,:,k);
    vp = v(:,:,k);
    s = sqrt(up.^2+vp.^2);
    
    % Plot Figure to Check Data
    clf
    pcolor(LON,LAT,vp);shading flat
    colormap(jet);caxis([-1.5 0]);hold on
    clo = core(:,n);clo(isnan(clo)) = [];
    cla = find(~isnan(core(:,n)'));
    plot(lon(clo),lat(cla),'ko','markerfacecolor','m');hold on
    plot(core_orig(:,n),lat,'ks','markerfacecolor','g')
    ylim([-30.85 -29.8]);title(k)
    
    % Left click = accept; 
    % Right click = reject and re-do interactively
    
    [~,~,button] = ginput(1);
    
    if button == 3 % (i.e left click = reject and re-do)
        
            % Run Function
            [coret,dirt,core_origt] = istream_conversion_p1(LON,LAT,up,vp);
            
            core(:,n) = coret;
            core_orig(:,n) = core_origt;
            dir(:,n) = dirt;
            
            clf
            pcolor(LON,LAT,vp);shading flat
            colormap(jet);caxis([-1.5 0]);hold on
            clo = core(:,n);clo(isnan(clo)) = [];
            cla = find(~isnan(core(:,n)'));
            plot(lon(clo),lat(cla),'ko','markerfacecolor','m');hold on
            pause
            
    elseif button == 2 % reject off the bat because clearly bad data
            
            core(:,n) = nan(size(lat))';
            core_orig(:,n) = nan(size(lat))';
            dir(:,n) = nan(size(lat))';
    end
    
    end
    
end

%save core_data_4yr_3 core core_orig dir k

%% Plot Information to Check Quality of Algorithm

figure('position',[100 60 900 900])

for n = 1:size(core,2)

    clf
    plot(LON,LAT,'k.');hold on
    clo = core(:,n);clo(isnan(clo)) = [];
    cla = find(~isnan(core(:,n)'));
    plot(lon(clo),lat(cla),'ro','markerfacecolor','r');hold on
    plot(core_orig(:,n),lat,'ks','markerfacecolor','g')
    pause(.3)

end

% %% Plot dir and core isnans
% 
% figure;
% plot(t(ts:te),sum(~isnan(dir),1))
% hold on
% plot(t(ts:te),sum(~isnan(core),1))
% 
% %% 
% 
% clo = round(nanmedian(core,2));
% clo(isnan(clo)) = [];
% cla = find(~isnan(nanmedian(core,2)));
% figure
% plot(lon(clo)',lat(cla),'ro','markerfacecolor','r');hold on

