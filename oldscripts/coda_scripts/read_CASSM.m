clear; close all; clc;
addpath(genpath('/home/spri902/scripts/coda_scripts'));
addpath(genpath('/home/spri902/EGS_Collab/joint_inv/'));
addpath(genpath('/home/spri902/Collab_metadata'));

%% Figure out the geometry for the data Parker wants me to analyze.

src_geom = readtable('SigmaV_Channel_Monitoring_Rev_1.0_PetrovMod.xlsx','sheet','CASSM_Channels');
rec_geom = readtable('SigmaV_Channel_Monitoring_Rev_1.0_PetrovMod.xlsx','sheet','Geode');

% now grab the specific sources and receivers of interest from entire set 
src_geom = src_geom(1:6,:);
rec_geom = rec_geom(13:23,:);
% plots connecting all sources and receivers...later we will color by dv/v
figure('Color','w');
hold on;
for i_src = 1 : size(src_geom,1)
    for i_rec = 1 : size(rec_geom,1)
        plot3(...
            [src_geom.Easting_m(i_src), rec_geom.Easting_m_(i_rec)],...
            [src_geom.Northing_m(i_src), rec_geom.Northing_m_(i_rec)],...
            [src_geom.Depth_m(i_src), rec_geom.Depth_m_(i_rec)],...
            'k');
        
    end
    plot3(src_geom.Easting_m(i_src), src_geom.Northing_m(i_src), ...
        src_geom.Depth_m(i_src), 'rp', 'MarkerFaceColor', 'r',...
        'MarkerSize', 12);
    text(src_geom.Easting_m(i_src), src_geom.Northing_m(i_src)-1,...
        src_geom.Depth_m(i_src), sprintf('%d',i_src))
end

for i_rec = 1 : size(rec_geom,1) 
    plot3(rec_geom.Easting_m_(i_rec), rec_geom.Northing_m_(i_rec),...
        rec_geom.Depth_m_(i_rec),'bv','MarkerFaceColor','b',...
        'MarkerSize', 12);
    text(rec_geom.Easting_m_(i_rec)+0.03, rec_geom.Northing_m_(i_rec),...
        rec_geom.Depth_m_(i_rec), sprintf('%d',i_rec));
end
xlabel('Easting [m]'); ylabel('Northing [m]'); zlabel('Depth [m]');
grid on; axis('xy');


%% Read in CASSM Geometry

fid    = fopen('Geode_01.txt'); % create file id
header = textscan(fid, '%s %s %s %s %s %s %s', 1, 'delimiter', ' ', ...
    'MultipleDelimsAsOne', 1, 'CollectOutput', 1); %scan for header
header = [header{:}]; % pull out the header

% read in the text file that has the test bed sensor geometry
geodeMap = readtable('Geode_01.txt', 'Format', '%u %s %s %f %f %f %f');
fclose(fid);

%% read in CASSM source metadata

sheet        = 'CASSM_Channels';
cassmGeom    = readtable('SigmaV_Channel_Monitoring_Rev_1.0_PetrovMod.xlsx', ...
    'Sheet', sheet);
src='PSB';
PSB_src_dep   = cassmGeom.Depth_m(    strcmp(src, cassmGeom.location));
PSB_src_enc   = cassmGeom.CytekCh(    strcmp(src, cassmGeom.location)); % use cytek channel not encoder ch
PSB_src_num   = cassmGeom.SourceNum(  strcmp(src, cassmGeom.location));
PSB_src_east  = cassmGeom.Easting_m(  strcmp(src, cassmGeom.location));
PSB_src_north = cassmGeom.Northing_m( strcmp(src, cassmGeom.location));
PSB_src_elv   = cassmGeom.Elev_m(     strcmp(src, cassmGeom.location));


%% Search for file number to locate OT CASSM sources in the data folders

% Use the function findCollabSource here
files = dir(fullfile('~/Desktop/Projects/EGS_Collab/collab_seismic/Data/CASSM/2018*'));% grab all data
datFolders = {files.name}'; % grab file numbers
datFolders = datFolders(~(endsWith({files.name}','.txt')));
dataStartTimes = datetime(datFolders, 'InputFormat', 'yyyyMMddHHmmss');
srchDir = '~/Desktop/Projects/EGS_Collab/collab_seismic/Data/CASSM/';% yours will be different

% srcNumstr  = ['00';'01';'02';'03';'04';'05'];
% thresh = 0.5;
% cytkCh = decodeCASSMChannel(t,data(:,94),thresh,0);
% srcChLoc = cassmGeom(find(cytkCh == cassmGeom.CytekCh),:);

%% Read in CASSM data and select only data from columns that are OT sensors

% OT hydrophone channels (13:23) in original data are 1:11 now
% OT accelerometer channels (70:78) in original data 12:20 now

% each cell contains the OT receiver traces from a single source near the
% middle of the OB bore hole directly beneath OT

srcFamily = ['PSB1'; 'PSB2'; 'PSB3'; 'PSB4'; 'PSB5'];
%data = cell(length(srcFamily) + 1, length(datFolders)); for OB plus 1 is
%bc 1 source is bad
data = cell(length(srcFamily) , length(datFolders));
t = data;
hydata = data;
acdata = data;

for k = 1:length(datFolders) % loop through the data folders
    
    [fileStruct,shotblocks] = findCollabSource2(srchDir, datFolders{k}); % find files for specific sources
    combFileStruct{k} = fileStruct; % #ok<SAGROW>
    sfileStruct = fileStruct( ismember([fileStruct.srcName], srcFamily) );
    shotblocks = shotblocks( ismember([fileStruct.srcName], srcFamily) );
    %[~,idx] = unique([fileStruct.fileName]','stable');
    %fileStruct = fileStruct(idx);
    files2use = [sfileStruct.fileName];
    
    if length(unique([sfileStruct.fileName]', 'stable')) < length(files2use)
        b = cell2mat(files2use');
        [~, ~, ridx] = unique(b, 'rows');
        rptFiles{k} = files2use(accumarray(ridx, 1) == 2);
        [~, ~, notes{k}] = readSeg2([srchDir, datFolders{k}, '/', char(rptFiles{k})], 'noplot'); %#ok<SAGROW>
        realTime = cellstr(strtrim(notes{k}.fileStrings{2}(end-8:end)));
        
        for ii = 1 : length(files2use)
            for jj = 1 : length(shotblocks{ii})
                shotTimes{ii, jj}  = strsplit(shotblocks{1, ii}{jj});
                shotTimes{ii, jj}  = shotTimes{ii, jj}{2};
            end
        end
        
        [~, loc] = ismember(shotTimes, realTime);
        [r, c] = find(loc);
        rptInds = find(r == ridx);
        files2use(rptInds(rptInds ~= r)) = [];
        sfileStruct(rptInds(rptInds ~= r)) = [];
    end
    
    %             times2srch = {fileStruct(ismember([fileStruct.fileName],rptFiles{j,k})).shotTime};
    %             times2srch=[times2srch{1,1},times2srch{1,2}];
    %             [~,loc]=ismember(times2srch,realTime);
    for j = 1 : length(sfileStruct) % loop through the different source files
        
        [data{sfileStruct(j).srcNum,k} , t{sfileStruct(j).srcNum,k}] = ...
            readSeg2([srchDir,datFolders{k},'/',char(files2use{j})], 'noplot');
        fs = floor(1/t{sfileStruct(j).srcNum, k}(2)-t{sfileStruct(j).srcNum, k}(1));
        dt = 1 / fs;
        datTrim = 0.03 * fs;
        
        % pull out the data
        data{sfileStruct(j).srcNum, k}   = data{sfileStruct(j).srcNum,k}(1:datTrim,[13:23]); % remember to change this %%%%%%
        t{sfileStruct(j).srcNum, k}      = t{sfileStruct(j).srcNum,k}(1:datTrim);
        % detrend and reassign the data (e.g. from 13:23 to 1:11 and 70:78 to 12:20)
        %         hydata{sfileStruct(j).srcNum, k} = detrend(data{sfileStruct(j).srcNum,k}(1:datTrim,1:11));
        %acdata{sfileStruct(j).srcNum, k} = detrend(data{sfileStruct(j).srcNum,k}(1:datTrim,12:20));
        
        hydata{sfileStruct(j).srcNum, k} = data{sfileStruct(j).srcNum,k}(1:datTrim,1:11);
        
    end
end

 hydata=hydata(any(~cellfun('isempty',hydata),2),:);
 t=t(any(~cellfun('isempty',t),2),:);
% OT/OB wells
% hydata (OB sources in bottom well x date) , 53 date files if we want
% instead of 3 I have now.

% Grab the individual hydrophone from the 11 working hydrophones

% May 22 is the first day
% May 24 is when the fracture broke through from injector to producer
% Pressure decreases after the last file a few hours later.

% hydata(row,col): row = srcNum, col = day (from the data folders)

% so hydata{1,1} is 1440x11 (npts, nrec) for source 1
% There are 5 sources:
% 'OB1'
% 'OB2'
% 'OB3'
% 'OB5'
% 'OB6'

%% Let's process each source

% Allocate some matrices and build a taper
npts      = datTrim;
ndays     = length(datFolders);
trc_data  = zeros(npts, ndays);
tim_data  = trc_data;
taper_win = tukeywin(npts, 0.05); % cosine taper on both ends

f_nyq = 1 / 2 / dt;
Wn = [2000, 3000] ./ f_nyq; % [Hz] low and high frequencies
[bf, af] = butter(4, Wn);

[n_src,n_days] = size(hydata);
[n_pts,n_rec] = size(hydata{1,1});

% allocat the dv/v matrix
dv_v_mat = NaN(n_src, n_rec, n_days);

% the reference will be the first day
maxLag = 200; % max nuber of points to search forward and backward (can be npts, just takes longer and is unrealistic)
b = 4; % b-value to limit strain

t_start = 0.01; % [s] time to start the polynomial fitting
[~,n_start] = min(abs(tim_data(:,1) - t_start));

for i_src = 1 : n_src
    for i_rec = 1 : n_rec
        for i_day = 1 : n_days
            rec_data = hydata{i_src,i_day};
            if ~isempty(rec_data)
                trc_data(:,i_day) = detrend(rec_data(:,i_rec)) .* taper_win; % apply detrend & taper
                tim_data(:,i_day) = t{i_src,i_day}; % get time data
            end
        end
        if isempty(rec_data); continue; end
        
        trc_data = filtfilt(bf, af, trc_data);
        
        for i_day = 1 : n_days
            err = computeErrorFunction( trc_data(:,i_day), trc_data(:,1), npts, maxLag ); % compute error function over lages
            
            % accumuluate error in FORWARD direction
            direction = 1; % direction to accumulate errors (1=forward, -1=backward)
            dist  = accumulateErrorFunction( direction, err, npts, maxLag, b ); % forward accumulation to make distance function
            stbar = backtrackDistanceFunction( -1*direction, dist, err, -maxLag, b ); % find shifts in backward direction
            stbarTime = stbar .* dt; % convert from samples to time
            %             tvec2 = tim_data(:,test) + stbarTime'; % make the warped time axis for the trace you are testing
            % linear polynomial fitting
            p = polyfit(tim_data(n_start:end,1), stbarTime(n_start:end), 1);
            dv_v_mat(i_src, i_rec, i_day) = -p(1);
        end
        fprintf('Finished source-receiver pair: %d-%d\n', i_src, i_rec);
    end
end

dv_v_perc = dv_v_mat*100;

% plot the percent velocity change

c_lim = max([abs(min(dv_v_perc(:))), max(dv_v_perc(:))]); % make a symmetric colorbar using largest positive or negative value

color_intensity = linspace(-c_lim, c_lim, 11); % create the color range
colors = jet(numel(color_intensity)); % make RGB color matrix

% plots connecting all sources and receivers...later we will color by dv/v
for i_day = 1 : n_days
    figure('Color','w');
    hold on; title(sprintf('Day: %d', i_day))
    for i_src = 1 : size(src_geom,1)
        for i_rec = 1 : size(rec_geom,1) 
            
            % find the correct color index for plotting this line
            [~,c_idx] = min(abs(color_intensity - dv_v_perc(i_src, i_rec, i_day)));
            % skip plotting if this was a NaN dv/v value
            if isnan(dv_v_perc(i_src, i_rec, i_day)); continue; end
            plot3(...
                [src_geom.Easting_m(i_src), rec_geom.Easting_m_(i_rec)],...
                [src_geom.Northing_m(i_src), rec_geom.Northing_m_(i_rec)],...
                [src_geom.Depth_m(i_src), rec_geom.Depth_m_(i_rec)],...
                'color', colors(c_idx,:)); % set the correct color
            
        end
        plot3(src_geom.Easting_m(i_src), src_geom.Northing_m(i_src), ...
            src_geom.Depth_m(i_src), 'rp', 'MarkerFaceColor', 'r',...
            'MarkerSize', 12);
        text(src_geom.Easting_m(i_src), src_geom.Northing_m(i_src)-1,...
            src_geom.Depth_m(i_src), sprintf('src_%d',i_src));
    end
    
    for i_rec = 1 : size(rec_geom,1) 
        plot3(rec_geom.Easting_m_(i_rec), rec_geom.Northing_m_(i_rec),...
            rec_geom.Depth_m_(i_rec),'bv','MarkerFaceColor','b',...
            'MarkerSize', 12);
        text(rec_geom.Easting_m_(i_rec)+0.03, rec_geom.Northing_m_(i_rec),...
            rec_geom.Depth_m_(i_rec), sprintf('%d',i_rec));
    end
    xlabel('Easting [m]'); ylabel('Northing [m]'); zlabel('Depth [m]');
    grid on; axis('xy');
    colormap(colors); caxis([color_intensity(1), color_intensity(end)])
    c = colorbar; ylabel(c,'dv/v [%]');
end
set(gca,'ZDir','Reverse');



%% Test stretching and DTW on single receiver.

% source
i_src = 2;
% over trace
i_trc = 4;
 % test trace (i.e. day number)
test = 38;

% Compare all days
for i_day = 1 : 53
    rec_data = hydata{i_src,i_day};
    trc_data(:,i_day) = detrend(rec_data(:,i_trc)) .* taper_win; % apply detrend & taper
    tim_data(:,i_day) = t{i_src,i_day}; % get time data
end

% Filter and plot
trc_data = filtfilt(bf, af, trc_data);
plot(tim_data, trc_data); hold on;
plot(tim_data(:,1), taper_win, 'k');
xlabel('Time (s)');
close(clf)

% Do windowed strecthing on the filtered data


% stretching parameters
twin=0.005; % [s] width of moving window
tstep = 0.0001; % [s] distance to move window
dVmax = 0.025; % [%] maximum percentage change
dV = 0.001; % [%] sampling of velocity change percentage

[C,epsArray,tSamp] = movingWinStretch( trc_data(:,1), trc_data(:,test), ...
    dt, twin, tstep, dVmax, dV );

figure;
ax1 = subplot(3, 1, 1);
plot(tim_data(:,1), trc_data(:,1),'r'); hold on;
plot(tim_data(:,test), trc_data(:,test),'b');
xlabel('Time (s)');
ax2 = subplot(3, 1, 2);
plot( tSamp, C ); ylabel('Corr. Coeff.'); %ylim([0.9 1]);
xlim([0 tim_data(end)]);
ax3 = subplot(3, 1, 3);
plot( tSamp, epsArray ); ylabel('\epsilon'); ylim([-dVmax dVmax]);
xlabel('Time (s)'); xlim([0 tim_data(end)]);

% link axes for zooming
linkaxes([ax1 ax2 ax3],'x');

% warping

maxLag = 200; % max nuber of points to search forward and backward (can be npts, just takes longer and is unrealistic)
b = 4; % b-value to limit strain

err = computeErrorFunction( trc_data(:,test), trc_data(:,1), npts, maxLag ); % compute error function over lages

% accumuluate error in FORWARD direction
direction = 1; % direction to accumulate errors (1=forward, -1=backward)
dist  = accumulateErrorFunction( direction, err, npts, maxLag, b ); % forward accumulation to make distance function
stbar = backtrackDistanceFunction( -1*direction, dist, err, -maxLag, b ); % find shifts in backward direction
stbarTime = stbar .* dt; % convert from samples to time
tvec2 = tim_data(:,test) + stbarTime'; % make the warped time axis for the trace you are testing


t_start = 0.01; % [s] time to start the polynomial fitting
[~,n_start] = min(abs(tim_data(:,1) - t_start));

figure;
ax1 = subplot(3, 1, 1);
plot(tim_data(:,1), trc_data(:,1),'r'); hold on;
plot(tim_data(:,test), trc_data(:,test),'b');
legend('Reference','Day 2');
xlabel('Time (s)');
ax2 = subplot(3, 1, 2);
plot(tim_data(:,1), trc_data(:,1),'r'); hold on;
plot(tvec2', trc_data(:,test),'b');
xlabel('Warped time (s)');
ax3 = subplot(3, 1, 3);
plot( tim_data(:,1), stbarTime ); hold on;
plot( tim_data(n_start:end,1), stbarTime(n_start:end), 'r');
ylabel('dt'); ylim([-2e-3 2e-3]); grid on;
xlabel('Time (s)'); xlim([0 tim_data(end)]);

% link axes for zooming
linkaxes([ax1 ax2 ax3],'x');
linkaxes([ax1 ax2],'y');
xlim([0 tim_data(end,1)]);

% linear polynomial fitting
p = polyfit(tim_data(n_start:end,1), stbarTime(n_start:end), 1);
dv_v = -p(1);

