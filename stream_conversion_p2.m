function [R,Us,Vs,perplon,perplat] = stream_conversion_p2(LON,LAT,up,vp,core,dir)

% INPUT
% LON/LAT = 2D matrices of longitude/latitude
% U/V = 2D matrix of u and v components of velocity
% core = for each timestep - index vector of core lon at each lat
% dir = direction of each identified core point (same size matrix as core)

% OUTPUT
% R = distance from core for each gridpoint
% Us/Vs = matrices of u/v in stream coordinates

% FUNCTION DESCRIPTION
% Get gridpoints along a perpendicular line from each identifed core point
% located in "stream_conversion_p1.m".  Get the u and v at each
% perpendicular gridpoint and get distance from core. 

% REQUIRES: N/A

% Written by M. Archer July 2015 ------------------------------------------

% Get lon/lat vectors
lon = LON(:,1);
lat = LAT(1,:);

ind = find(~isnan(core)); % Find indices of core (i.e. remove NaN values)
nn = 0;

%keyboard

% % Check Data
% figure
% pcolor(LON,LAT,vp);shading flat
% hold on
% plot(lon(core(~isnan(core))),lat(~isnan(core)),'ko','markerfacecolor','k')
% 
% figure;n=3;
% quiver(LON(1:n:end,1:n:end),LAT(1:n:end,1:n:end),up(1:n:end,1:n:end),vp(1:n:end,1:n:end),3)
% hold on;plot(lon(core(~isnan(core))),lat(~isnan(core)),'ko','markerfacecolor','r')

if (~isempty(ind))
    
%% -------------------------------------------------------------> Stage 1
% For each core point identified, get the perpenicular line through it -
% (using dir of core variable 'dir').  Finish up with the perp lines' 2 end points.

    for ii = 1:length(ind)%ind(1):ind(end) % From first identified core point to last
        
        i = ind(ii);
               
        nn = nn + 1; % the core point # identified, counting from 1
        
        % Get core grid point at i - copy location to every grid point
        cy = lat(i).*ones(size(LON,1),size(LON,2)); % Value over entire matrix
        cx = lon(core(i)).*ones(size(LON,1),size(LON,2)); % Value over entire matrix
        
        % Get angle of each gridpoint in matrix to core(i) (figure;pcolor(LON,LAT,azy);shading interp)
        azy = azimuth(cy,cx,LAT,LON,'degrees'); % Calculate angle between every point and core point
       
        % Get desired perpendicular angle value 
        ep = (dir(i)-90).*ones(size(LON,1),size(LON,2)); % What is the desired angle? To EAST
        ep(ep>360) = ep(ep>360)-360; EP(nn) = ep(1,1); % for debugging purposes
        %
        wm = (dir(i)+90).*ones(size(LON,1),size(LON,2)); % What is the desired angle? To WEST
        wm(wm<0) = wm(wm<0)+360; WM(nn) = wm(1,1);
        
        % Get difference matrix - will be looking for zero (i.e. get same)
        east = abs(azy-ep);west = abs(azy-wm);
        EAST(:,:,nn) = east; WEST(:,:,nn) = west;
        
        % East gridpoint
        k=1;
        [m(k),in(k)] = min(east(end,:));k = 2; % i.e. right hand wall
        [m(k),in(k)] = min(east(:,end));k = 3; % i.e. bottom wall
        [m(k),in(k)] = min(east(:,1)); % i.e. top wall
        [~,index] = min(m);
        if index == 1
            egpi(:,nn) = [length(lon) in(index)]; % lon/lat indices
            egp(:,nn) = [lon(end) lat(egpi(2,nn))]; % lon/lat values
        elseif index == 2
            egpi(:,nn) = [in(index) length(lat)];
            egp(:,nn) = [lon(egpi(1,nn)) lat(end)];
        else
            egpi(:,nn) = [in(index) 1];
            egp(:,nn) = [lon(egpi(1,nn)) lat(1)];
        end
        clear k m in index
        
        % West gridpoint
        k=1;
        [m(k),in(k)] = min(west(1,:));k = 2; % i.e. left wall
        [m(k),in(k)] = min(west(:,1));k = 3; % i.e. top wall
        [m(k),in(k)] = min(west(:,end)); % i.e. bottom wall
        [~,index] = min(m);
        if index == 1
            wgpi(:,nn) = [1 in(index)]; % lon/lat indices
            wgp(:,nn) = [lon(1) lat(wgpi(2,nn))]; % lon/lat values
        elseif index == 2
            wgpi(:,nn) = [in(index) 1];
            wgp(:,nn) = [lon(wgpi(1,nn)) lat(1)];
        else
            wgpi(:,nn) = [in(index) length(lat)];
            wgp(:,nn) = [lon(wgpi(1,nn)) lat(end)];
        end
        clear k m in index

    end
    %clear k n m ind index indy i ii nn east west ep wm azy ang cy cx
    
%     % Now you have wgp and egp (2*length(core(ind))) - i.e. lon/lat for
%     each core point
%     figure;
%     plot(egp(1,:),egp(2,:),'ro');hold on
%     plot(wgp(1,:),wgp(2,:),'bo');hold on
%     plot(lon(core(~isnan(core))),lat(~isnan(core)),'ko','markerfacecolor','k');hold on
%     quiver(lon(core(~isnan(core))),lat(~isnan(core)),...
%         up(core(~isnan(core)),~isnan(core)),vp(core(~isnan(core)),~isnan(core)))
    
%% -------------------------------------------------------------> Stage 2 
% Least square fit the line from one edge point to the other

    n = 0;
    LSQlon = nan(length(lon),length(egp)); 
    LSQlat = nan(length(lat),length(egp)); 
    latloc = nan(length(lon),length(egp)); % Initialize
    lonloc = nan(length(lat),length(egp));
    mv = nan(length(egp),1);mav = nan(length(egp),1); % this is a northernmost/southernmost longitude index
    wi = nan(length(egp),1);ei = nan(length(egp),1);% this is a westernmost/easternmost longitude index
    
    for i = 1:length(egp) %egp = east gridpoint, identified far E edge of perpendicular line for each core point
        n = n+1; % counting from 1
        WI = wgpi(1,i); EI = egpi(1,i); wni = wgpi(2,i); esi = egpi(2,i); 
        lonnumber(n) = length(WI:EI); % how many longitude points between ends of line?
        latnumber(n) = length(wni:esi) + length(esi:wni); % how many latitude points between ends of line?      
        lo(:,n) = [wgp(1,i); egp(1,i)];
        la(:,n) = [wgp(2,i); egp(2,i)];
        
        if (latnumber(n) > lonnumber(n)) %then use lat to regress
            % Switch lo/la because other way around to if using unqiue lon
            bls = regress(lo(:,n),[ones(size(la(:,n))) la(:,n)]); % Regress line 
            [mv(n),~] = min([wni esi]);[mav(n),~] = max([wni esi]);
            lsq = bls(1)+bls(2)*lat(mv(n):mav(n)); % Least-Squares
            LSQlat(mv(n):mav(n),n) = lsq; % Save out of loop
            NN = 0;
            for ii = mv(n):mav(n)
                NN = NN + 1;
                lonnt = lsq(NN); %lat not true
                diflon = lon-lonnt*(ones(length(lon),1)); %diff between not true lsq and lat
                [~,indx] = min(abs(diflon)); % find closest lat to lsq not true lat
                lonloc(ii,n) = indx; % True Latitudes for the lsq line
            end
            
        else
            ei(n) = EI;wi(n) = WI;  
            bls = regress(la(:,n),[ones(size(lo(:,n))) lo(:,n)]); % Regress line
            lsq = bls(1)+bls(2)*lon(wi(n):ei(n)); % Least-Squares
            LSQlon(wi(n):ei(n),n) = lsq; % Save out of loop
           % Now get true latitude points
            %for ii = 1:length(lsq)
            NN = 0;
            for ii = wi(n):ei(n)
                NN = NN + 1;
                latnt = lsq(NN); %lat not true (i.e. not a real grid point)
                diflat = lat-latnt*(ones(1,length(lat))); %diff between not true lsq and lat
                [~,indx] = min(abs(diflat)); % find closest lat to lsq not true lat
                latloc(ii,n) = indx; % True Latitude indices for the lsq line
            end
            
        end
        
    end
    
%     figure;plot(lon,LSQlon,'b.');hold on;plot(LSQlat,lat,'r.') 
%     %latloc (or lonloc) is a proxy for number of gridpoints along perp line
%     figure;plot(lon(core(~isnan(core))),lat(~isnan(core)),'ro')
%     % - note the points arent EXACTLY aligned with core since not using latloc
%     hold on
%     plot(lon,LSQlon,'b.')
% %     hold on
% %     plot(lon,latloc)
%     
%% -------------------------------------------------------------> Stage 3
% Extract the actual lon/lat of every gridpoint in relation to its' core
% point

    % i.e. get lon and lat indices of core points removing any nans
    core_loni = core; core_loni(isnan(core_loni)) = [];
    core_lati = find(~isnan(core)); % i.e same as ind
    
    %Initialize
    perplon = nan(length(lat),length(core_loni)); % Value
    perplat = nan(length(lat),length(core_loni)); % Value
    perploni = nan(length(lat),length(core_loni));% Index
    perplati = nan(length(lat),length(core_loni));% Index
    
    for i = 1:length(core_loni) % from 1 to number of identified core points
        
        clon(i) = lon(core_loni(i)); % core longitude value
        clat(i) = lat(core_lati(i)); % core latitude value
        
        if sum(~isnan(latloc(:,i))) == 0 %i.e. if lonloc needs to be used
        perpvec = lonloc(:,i); 
        ind_p = find(~isnan(perpvec));
        % Longitude
        perplon(1:length(ind_p),i) = lon(perpvec(~isnan(perpvec)));
        perploni(1:length(ind_p),i) = perpvec(~isnan(perpvec));
        % Latitude
        perplat(1:length(ind_p),i) = lat(ind_p);
        perplati(1:length(ind_p),i) = ind_p;
        % Check the core lines with the perpvec LON only
        [~,minind] = min(abs(perplat(:,i)-(clat(i)*ones(length(perplat(:,i)),1))));
        perplon(minind,i) = clon(i);    
                
        else
        perpvec = latloc(:,i); 
        ind_p = find(~isnan(perpvec));
        % Longitude
        perplon(1:length(ind_p),i) = lon(ind_p); %should I switch to perplon(ind_p,i) = ?
        perploni(1:length(ind_p),i) = ind_p;
        % Latitude
        perplat(1:length(ind_p),i) = lat(perpvec(~isnan(perpvec)))';
        perplati(1:length(ind_p),i) = perpvec(~isnan(perpvec));
        % Check the core lines with the perpvec LAT only
        [~,minind] = min(abs(perplon(:,i)-(clon(i)*ones(length(perplon(:,i)),1))));
        perplat(minind,i) = clat(i);
        
        end
        
    end
    %clear i ind_p clon clat
    perploni(perploni==0)= NaN;
    perplati(perplati==0)= NaN;
    
%     figure;plot(lon(core(~isnan(core))),lat(~isnan(core)),'ro')
%     hold on
%     plot(perplon,perplat,'.')
    
%% -------------------------------------------------------------> "stage3.m"
    
    % Get u and v at perp gps - initialize
    U = nan(size(perploni,1),size(perploni,2));
    V = nan(size(perploni,1),size(perploni,2));
    Us = nan(size(perploni,1),size(perploni,2));
    Vs = nan(size(perploni,1),size(perploni,2));
       
    d = deg2rad(true2math(dir(~isnan(core)))); % convert direction back to math/radians
        
    for i = 1:size(perploni,1)
        for j = 1:size(perploni,2)
            
            if ~isnan(perploni(i,j))
                U(i,j) = up(perploni(i,j),perplati(i,j));
                V(i,j) = vp(perploni(i,j),perplati(i,j));
                Us(i,j) = -(U(i,j).*sin(d(j)) - V(i,j).*cos(d(j))); % convert to down/cross stream
                Vs(i,j) = -(U(i,j).*cos(d(j)) + V(i,j).*sin(d(j))); % convert to down/cross stream
            end
            
        end
    end
    
    % Get distance between perpgps and core gp
    E = nan(size(perploni,1),size(perploni,2));
    N = nan(size(perploni,1),size(perploni,2));
    R = nan(size(Us,1),size(Us,2));
    
    for i = 1:size(Us,1)
        for j = 1:size(Us,2)
            
            if ~isnan(Us(i,j))
             [E(i,j),N(i,j)] = lonlat2km(lon(core_loni(j)),lat(core_lati(j)),lon(perploni(i,j)),lat(perplati(i,j)));
            end
            
            r = sqrt(E(i,j).^2+N(i,j).^2); % r = distance from core
            r(E(i,j)<0) = -r;
            R(i,j) = r;
            
        end
    end

else
    first = length(lat);
    second = length(lat);
    R = nan(first,second);
    Us = nan(first,second);
    Vs = nan(first,second);
end

end    
