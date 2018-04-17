% Stat 5525 Course Project: find OCB via machine learning.
clear; clc; close all;
stations = ['THL'; 'SVS'; 'KUV'; 'UPN'; 'UMQ'; 'GDH'; 'ATU'; 'STF'; 'SKT'; 'GHB'; 'FHB'; 'NAQ']; % DTUS station names
comps = ['Bx'; 'By'; 'Bz'];

tresl = 10; nrecs = 86400/tresl; 
ttime = (1:nrecs)*tresl;
nstat = length(stations);

% Specify date range 'yyyymmdd'
sday = datenum('20140101','yyyymmdd');
eday = datenum('20140110','yyyymmdd');
day = eday-sday+1;
B_set = zeros(day*24*12,360);

% Construct data_set
for iday = sday:eday
    iday-sday+1
    ymd = datestr(iday,'yyyymmdd');
    Bx = zeros(nrecs,nstat)+NaN; By = Bx; Bz = Bx; B = Bx;
    for istat=1:nstat
        station = stations(istat,:); 
        path0 = [pwd '\data\'];
        pathtxt = [path0 station '_txt\' station '_2014\'];
        staymd = [lower(station), ymd];
        if (exist([pathtxt, staymd, 'Bx.txt'],'file'))
            Bx = load([pathtxt, staymd, 'Bx.txt']); 
            By = load([pathtxt, staymd, 'By.txt']);
            Bz = load([pathtxt, staymd, 'Bz.txt']);  
            B = sqrt(Bx.^2+By.^2+Bz.^2);
        else
            B = zeros(nrecs,nstat)+NaN;
        end
        for itime = 1:24
            set_index = (iday-sday)*288 + (itime-1)*12 +istat;
            B_set(set_index,:) = B((itime-1)*360+1:itime*360);
        end
    end
end

%Normalization
B_max = max(max(B_set));
B_min = min(min(B_set));
B_set = (B_set-B_min)/(B_max-B_min);

%kmeans
kmeans_index = kmeans(B_set,3);d
kmeans_index = reshape(kmeans_index,12,day*24);



