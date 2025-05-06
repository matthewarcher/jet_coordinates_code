% RUN STREAM CONVERSION PART (2)
% Get gridpoints along a perpendicular line from each identifed core point
% located in "stream_conversion_p1.m".  Get the u and v at each
% perpendicular gridpoint and get distance from core. 

clear;clc

% Pathways
cd /home/nfs/z3519950/sci-maths-ocean/Data/RADAR/VECTORS
addpath /home/nfs/z3519950/hdrive/hf_radar/Meandering

% Data
load CoffsHarbour_2012_2016_Prepped
load core_data_2012_2016; U=u; V=v;

% %% Test Subset
% nn = 1000;
% t = t(1:nn);
% U = u(:,:,1:nn);
% V = v(:,:,1:nn);
% core = core(:,1:nn);
% core_orig = core_orig(:,1:nn);
% dir = dir(:,1:nn);

%% Run Time Loop here

% Initialize Variables to Save
longest = max([length(lat) length(lon)]);
first = longest;
second = longest;
third = length(t); 
R = nan(first,second,third);US = nan(first,second,third);VS = nan(first,second,third);
perplon = nan(first,second,third);perplat = nan(first,second,third);

tic
for KK = 1:third
    
    display([num2str(KK) '/' num2str(third)])
    
    %  Get values for timestep
    c = core(:,KK); % Core location
    d = dir(:,KK); % Direction of core
    u = U(:,:,KK); % u velocity
    v = V(:,:,KK); % v velocity
    
    if (sum(~isnan(c)) > 0)
        
    % Run Function
    [r,us,vs,PLO,PLA] = stream_conversion_p2(LON,LAT,u,v,c,d);
    
    % Prepare Variables in Time
    ind = find(~isnan(c));
    R(1:length(r),ind,KK) = r; 
    VS(1:length(r),ind,KK) = vs; clear vs
    US(1:length(r),ind,KK) = us; clear us
    perplon(1:length(r),ind,KK) = PLO; clear PLO
    perplat(1:length(r),ind,KK) = PLA; clear PLA r
    
    end    
end
toc

VS(VS==0)=NaN;
US(US==0)=NaN;

save stream_data_EAC_2012_2016 R VS US perplon perplat LON LAT lon lat U V core dir core_orig

%% Test Plot
figure('position',[60 100 900 900])
along = 1:1:size(VS,2);
AS = repmat(along,[length(lat), 1]);

for n = 1:1:third
    if sum(sum(~isnan(VS(:,:,n)))) ~= 0 
    clf
    subplot(1,2,1);pcolor(LON,LAT,sqrt(V(:,:,n).^2));
    shading flat;axis square;colormap(jet);caxis([0 1.75])
    hold on
    subplot(1,2,2);pcolor(R(:,:,n),AS,sqrt(VS(:,:,n).^2));
    shading flat;axis square;colormap(jet);caxis([0 1.75])
    %xlim([-50 90]);ylim([1 80])
    pause
    end
end  
close

% Plot just gridpoints
figure('position',[60 100 900 900])
along = 1:1:size(VS,2);
AS = repmat(along,[length(lat), 1]);

for n = 1:third
    if sum(sum(~isnan(VS(:,:,n)))) ~= 0 
        plot(R(:,:,n),AS,'.');pause
    end
end
close
    
%% Plot in Increments
figure

for n = 1:third % Choose Time
    
    % Initialize Variables
    RR = nan(size(R,1),size(R,2));
    ASS = nan(size(R,1),size(R,2));
    VSS = nan(size(R,1),size(R,2));
    along = size(VS,2):-1:1;
    AS = repmat(along,[length(lon), 1]);
    
    % Plot Figure
    figure('position',[60 100 900 900])
    subplot(1,2,1);
    pcolor(LON,LAT,sqrt(V(:,:,n).^2+U(:,:,n).^2));
    shading interp;axis square;colormap(jet);caxis([0 1.75])

    hold on
    ccore = core(:,n);
    pplon = perplon(:,:,n);
    pplat = perplat(:,:,n);
    ind = find(~isnan(ccore)); % Index of core points

    for kk = ind(end):-1:ind(1)
        
        % Left Side
        subplot(1,2,1)
        hold on
        plot(lon(ccore(kk)),lat(kk),'ko');
        hold on
        plot(pplon(:,kk),pplat(:,kk),'k.')  
        hold on
        quiver(lon(ccore(kk)),lat(kk),...
            U(ccore(kk),kk,n),...
            V(ccore(kk),kk,n),0.0005,'k','linewidth',2)
        
        % Right Side
        subplot(1,2,2);hold on;ylim([1 length(along)])
        RR(:,kk) = R(:,kk,n);
        ASS(:,kk) = AS(:,kk);
        VSS(:,kk) = VS(:,kk,n);
        
        pcolor(RR,ASS,sqrt(VSS.^2));
        shading interp;axis square;colormap(jet);caxis([0 1.75])
        
        pause
    end
    close
end

% Plotting Vectors
k = find(t == datenum(2013,2,2,1,0,0))
figure;n=3;
quiver(LON(1:n:end,1:n:end),LAT(1:n:end,1:n:end),U(1:n:end,1:n:end,k),V(1:n:end,1:n:end,k),5)

figure;n=3;
quiver(R(1:n:end,1:n:end,k),AS(1:n:end,1:n:end),US(1:n:end,1:n:end,k),VS(1:n:end,1:n:end,k),5)

figure;n=3;
quiver(R(1:n:end,1:n:end,k),AS(1:n:end,1:n:end),-US(1:n:end,1:n:end,k),-VS(1:n:end,1:n:end,k),5)
