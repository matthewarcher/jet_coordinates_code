function [core,dir,core_lon] = stream_conversion_p1(LON,LAT,up,vp,datatype)

% AUTOMATED VERSION 
% UNSW Australia 2016

% INPUT
% LON/LAT = matrices of longitude/latitude
% U/V = 3D matrix of u and v components of velocity (3rd dim = time)

% OUTPUT
% core = for each timestep - index vector of smoothed core lon at each lat
% dir = direction of each identified core point (same size matrix as core)
% core_lon = raw identified core for archival/QC purposes

% FUNCTION DESCRIPTION
% Find core of Florida Current using maximum vecloity values, automatically
% pick region of good data, find direction of core at chosen locations.

% REQUIRES: "findcore.m"

% Written by M. Archer July 2015 ------------------------------------------

if exist('datatype','var')==0
    datatype = 1; % default to HF radar
end

% keyboard
% Get lon/lat vectors
lon = LON(:,1);
lat = LAT(1,:);

s = (up.^2+vp.^2).^0.5; % Get magnitude

% Identify the maximum
[~,core_lon,~]=findcore(lon,lat,up,vp,1,datatype); % CHECK THRESHOLDS WITHIN THIS SCRIPT
core_lon(core_lon == lon(1)) = NaN;

if  (sum(~isnan(core_lon)) > 5)  %(sum(abs(diff(core_lon))) ~= 0) &&
    
% Create datafile of interactively chosen core position indices
icore = core_lon;
icore(icore==lon(1))=NaN;

% Spline fit to interactively chosen core position
[score,~] = csaps(lat(~isnan(icore)),icore(~isnan(icore)),0.995,lat); % Choose roughness parameter
sCORE = nan(size(lat));
sCORE(~isnan(icore)) = score(~isnan(icore)); % Take only points where core identified

% Get a longitude index for each sCORE value
for i = 1:length(lat)
    [~,index(i)]=min(abs(lon-sCORE(i))); % I.e. get closest points
end

core=index;
core(core == 1) = NaN;

%%%%%%%%%%%%%% YOU HAVE NOW IDENTIFIED THE CORE REGION %%%%%%%%%%%%%%%%
%%%%%%%%%%%%%% Take Away: core %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Determine direction of downstream flow
d = nan(size(lat));
ind = find(~isnan(core));
for i = ind(1):ind(end)
    if ~isnan(core(i))
        [th,~] = cart2pol(up(:,i),vp(:,i));
        D = th(core(i)-1:core(i)+1); % Vector of directions
        dr = nanmedian(D); % median direction in radians
        d(i) = math2true(rad2deg(dr));
%     else
%         [th,~] = cart2pol(up(:,i),vp(:,i));
%         dr = th(core(i));
%         d(i) = math2true(rad2deg(dr));
    end
end

dir = d;

%%%%%%%%%%%%% YOU HAVE NOW IDENTIFIED THE CORE DIRECTION %%%%%%%%%%%%%%
%%%%%%%%%%%%% Take Away: dir %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

else 
   core = nan(size(lat));
   dir = nan(size(lat));
end

end
