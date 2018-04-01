
%% The first part of this script loads and pre-processes 3D time series of calcium imaging
% This part is specific for the properties of the time series and must be modified according to one's needs.
%% The second part (starting with "%% Select ROIs") uses the 3D time series to load a non-graphical user interface
% This part is based on a couple of helper scripts in the background that
% make drawing ROIs enjoyable.

clear all
global clut2b timetracesX_X ROI_map_X movie_AVG_X
addpath('non-GUI ROI analysis\')
load clut2b


%% Load list of files
% This code is optimized to load images acquired with Scanimage B
% (https://github.com/PTRRupprecht/Instrument-Control), but can be easily
% modified to load any tif-based time series

Filename = 'ca_imaging_example.tif';
clear meta
[A,result,meta.framerate,meta.zstep,meta.zoom,meta.motorpositions,meta.scalingfactors] = read_metadata_function(Filename);
L = imfinfo(Filename);
meta.height = L(1).Height;
meta.width = L(1).Width;
meta.numberframes = numel(L);

%% Read raw data from hard disk
meta.framerate = meta.framerate;
meta.framerate = meta.framerate;
nb_frames_total = sum(floor(meta.numberframes));
nb_planes = 1;
startingpoint = 1;
binning = 1;

movie = read_movie(Filename,meta.width,meta.height,nb_frames_total,startingpoint,binning,L,nb_planes);


%% Subtract baseline from movie
baseline = quantile(movie(:),0.03);
movie = movie - baseline;

%% Calculate average and maps of local correlations and responses

% average
AVG_movie = mean(movie,3);

% calculate activity maps
offset = 0;
f0_window = [1 300];
response_window = [418 min(550,size(movie,3))];
plot1 = 0; plot2 = 0; DF_movie_yesno = 0; % figure numbers
[DF_reponse,DF_master,DF_movie] = dFoverF(movie,offset,meta.framerate,plot1,plot2,DF_movie_yesno,f0_window,response_window);

% local correlation map (computational slightly expensive, but helpful)
tilesize = 16;
localCorrelations(:,:) = localCorrelationMap(movie,tilesize);

figure(88), 
subplot(1,2,1);
imagesc(DF_reponse(:,:),[-0.5 2]); axis off equal
subplot(1,2,2);
imagesc(localCorrelations(:,:),[-0.5 3]); axis off equal
akZoom('all_linked')


%% Select ROIs - this is the main part of the program, allowing to select ROIs
ROI_map_input = zeros(size(AVG_movie(:,:)));
trial_nb = 1; 
offset = 0;
df_scale = [-20 100];
AVG_Z = AVG_movie;
% AVG_Z(AVG_Z>80) = 80;

[ROI_mapX,timetracesX,timetracesX_raw] = timetraces_singleplane(movie,AVG_Z,offset,DF_reponse(:,:),DF_master(:,:),localCorrelations(:,:),df_scale,ROI_map_input,meta,1,AVG_Z);
% figure, imagesc(conv2(timetracesX,fspecial('gaussian',[25 1],23),'same'));
ROI_map_input = ROI_mapX;

plane{1}.ROI_map = ROI_mapX;
plane{1}.timetraces = timetracesX;
plane{1}.timetraces_raw = timetracesX_raw;
plane{1}.DF_reponse = DF_reponse;
plane{1}.meta = meta;
plane{1}.anatomy = AVG_movie;

save('Extracted_Data.mat','plane');

