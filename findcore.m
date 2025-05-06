function [core_idx,core_lon,mx] = findcore(lon,lat,up,vp,spar,datatype)

% <<< INPUT
% LON/LAT = longitude/latitude vectors
% u/v = u/v 2D matrices (at a single time frame)
% spar = smoothing parameter (# of gridpoints)
% datatype ---> if HF radar = 1, if ADCP = 2

% >>> OUTPUT
% core = vector (length of lat) with location of max value along longitude

% This script loops through all gridpoints in longitude to identify
% the FC core for each step in latitude.

% ---------------------------------- Written by Matt Archer; 22nd June 2015 
%

if exist('datatype','var')==0
    datatype = 1;
end
    
% Prepare Variables
sp = (up.^2 + vp.^2).^0.5;  % Calculate magnitude 
[LAT,LON] = meshgrid(lat,lon); % Create matrix of lat/lon
core_idx = nan(length(lat),1); % Initialize
mx = nan(length(lat),1); % Initialize
p1 = spar; p2 = spar-1; % Smoothing parameter

% For each latitude: Smooth Lateral Profile and identify core 
%!- why am I smoothing here and not before hand properly? -!
SH = nan(size(sp));
for j = 1:length(lat)
    
    s = sp(:,j);
    
    for i = p1:length(lon)-p2
        
        tS = squeeze(s(i-p2:i+p2));
        div = sum(~isnan(tS)); % To divide by
        
        if div > 0
            SH(i,j) = sum(tS)/div;
        end       

        if datatype == 1
        if vp(i,j) > -.5%-0.75
            SH(i,j) = NaN;
        end
        end
        
    end
    
    %SH(:,vp(:,j)>0) = NaN;

    [mx(j),core_idx(j)] = nanmax(SH(:,j));

end
core_lon = lon(core_idx);

end
