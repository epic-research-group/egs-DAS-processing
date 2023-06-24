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
set(gca,'ZDir','reverse');

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
OB_src_dep   = cassmGeom.Depth_m(    strcmp(src, cassmGeom.location));
OB_src_enc   = cassmGeom.CytekCh(    strcmp(src, cassmGeom.location)); % use cytek channel not encoder ch
OB_src_num   = cassmGeom.SourceNum(  strcmp(src, cassmGeom.location));
OB_src_east  = cassmGeom.Easting_m(  strcmp(src, cassmGeom.location));
OB_src_north = cassmGeom.Northing_m( strcmp(src, cassmGeom.location));
OB_src_elv   = cassmGeom.Elev_m(     strcmp(src, cassmGeom.location));


%% Search for file number to locate OT CASSM sources in the data folders

% Use the function findCollabSource here
files = dir(fullfile('/home/spri902/EGS_Collab/joint_inv/2018*'));% grab all data
datFolders = {files.name}'; % grab file numbers
%datFolders = datFolders(~(endsWith({files.name}','.txt')));
dataStartTimes = datetime(datFolders, 'InputFormat', 'yyyyMMddHHmmss');
srchDir = '/home/spri902/EGS_Collab/joint_inv/';% yours will be different

% srcNumstr  = ['00';'01';'02';'03';'04';'05'];
% thresh = 0.5;
% cytkCh = decodeCASSMChannel(t,data(:,94),thresh,0);
% srcChLoc = cassmGeom(find(cytkCh == cassmGeom.CytekCh),:);

%% Read in CASSM data and select only data from columns that are OT sensors

% OT hydrophone channels (13:23) in original data are 1:11 now
% OT accelerometer channels (70:78) in original data 12:20 now

% each cell contains the OT receiver traces from a single source near the
% middle of the OB bore hole directly beneath OT
for i = 1:length(datFolders)
    files2use = dir([srchDir,datFolders{i},'/*.dat']);
    for k = 1:length(files2use)
        
        [data{i,k} , t{i,k}] = ...
            readSeg2([srchDir,datFolders{i},'/',char(files2use(k).name)], 'noplot');
        
        src_Chan=decodeCASSMChannel(t{k},data{i,k}(:,94),0.5,0);
    end
end