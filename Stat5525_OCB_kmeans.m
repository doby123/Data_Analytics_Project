% Stat 5525 Course Project: find OCB via machine learning.
clear; clc; close all;
path0 = [pwd '/data/'];
stations = ['THL'; 'SVS'; 'KUV'; 'UPN'; 'UMQ'; 'GDH'; 'ATU'; 'STF'; 'SKT'; 'GHB'; 'FHB'; 'NAQ']; % DTUS station names
latDTUm = [84.40, 82.68, 80.36, 78.57, 75.99, 74.82, 73.54, 72.14, 70.93, 69.49, 66.92, 65.23];
lonDTUm = [27.48, 31.23, 40.28, 38.71, 41.22, 38.15, 37.09, 39.96, 36.43, 37.12, 38.43, 42.61];
sday = datenum('20140101','yyyymmdd');
eday = datenum('20140101','yyyymmdd');
Bcomp = 'Bx';

tresl = 10; nrecs = 86400/tresl; 
nstat = length(stations);
nstep = 300/tresl; % 300s segements as one point to be clustered.
nsegs = nrecs/nstep; % # of segments with length of nstep
ndays = eday-sday+1;
nclts = nsegs*nstat; % number of folded data points for clustering
% B_set = zeros(nday*24*12,360);

rstat = cosd(latDTUm')*ones(1,nsegs);
tstat = ones(nstat,1)*(0:nsegs-1)*360/nsegs + lonDTUm'*ones(1,nsegs);
xstat = rstat.*cosd(tstat); ystat = rstat.*sind(tstat);

rstat0 = cosd(latDTUm')*ones(1,nrecs);
tstat0 = ones(nstat,1)*(0:nrecs-1)*360/nrecs + lonDTUm'*ones(1,nrecs);
xstat0 = rstat0.*cosd(tstat0); ystat0 = rstat0.*sind(tstat0);

% Construct data_set
for iday = sday:eday
    iday-sday+1
    ymd = datestr(iday,'yyyymmdd');
    Bx = zeros(nrecs,nstat)+NaN; By = Bx; Bz = Bx;
    for istat=1:nstat
        station = stations(istat,:); 
        pathtxt = [path0 station '_txt/' station '_2014/'];
        staymd = [lower(station), ymd];
        if (exist([pathtxt, staymd, 'Bx.txt'],'file'))
            Bx(:,istat) = load([pathtxt, staymd, 'Bx.txt']); 
            By(:,istat) = load([pathtxt, staymd, 'By.txt']);
            Bz(:,istat) = load([pathtxt, staymd, 'Bz.txt']);  
        end
    end
    Bt = sqrt(Bx.^2+By.^2+Bz.^2);
    if(strcmp(Bcomp,'Bx'))
        B_ana = Bx; % B_ana is the variable to be analyzed, can be Bt, Bx, By, or Bz etc.
    end
%     for itime = 1:24
%         set_index = (iday-sday)*288 + (itime-1)*12 +istat;
%         B_set(set_index,:) = Bt((itime-1)*360+1:itime*360);
%     end
    B_clt = B_ana-ones(nrecs,1)*mean(B_ana); % remove mean of each station
    B_set = reshape(B_clt(:),nstep,nclts)'; % nclts = nsegs*ndays*nstat

    %Normalization
    B_max = max(B_set(:));
    B_min = min(B_set(:));
    B_set = (B_set-B_min)./(B_max-B_min);
    %kmeans
    kmeans_index = kmeans(B_set,3);
    kmeans_plot = reshape(kmeans_index,nsegs,nstat)'; % (nclts/nstat) columns x nstat rows

    % plot cluster distribution.
    fig = figure('position',[50,50,800,380],'paperpositionmode','auto'); 
    subplot(1,2,1); hold on; box on; axis equal;
    pcolor(xstat0,ystat0,B_clt'); shading flat; colorbar;
    contour(xstat0,ystat0,B_clt',[-50,0,50],'k');
    title([ymd,' ',Bcomp]);

    subplot(1,2,2); hold on; box on; axis equal;
    pcolor(xstat,ystat,kmeans_plot); shading flat; colorbar;
    contour(xstat,ystat,kmeans_plot,[1.5,2.5],'k');
    title([ymd,' Cluster']);

    print(fig,['3_cluster_',Bcomp,'_Nstep',num2str(nstep,'%4.4d'),'_',ymd],'-dpng');
    close(fig);
end
%------------------ THE END --------------------