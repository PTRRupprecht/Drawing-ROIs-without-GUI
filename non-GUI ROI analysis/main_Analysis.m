clear all
global clut2b timetracesX_X ROI_map_X movie_AVG_X
load clut2b
% load list of files
% cd distorted
FileList_0 = dir('Sag_Pos3_1-51_001*_.tif');
% cd ..
FileList = FileList_0;%dir('OB_4planes_paradigm1*_undistorted.tif');
FileListX = FileList;
% for k = 1:numel(FileList)
%     FileList(k) = FileListX(numel(FileList)+1-k);
% end
% cd distorted
% read metadata
clear meta
[A,result,meta.framerate,meta.zstep,meta.zoom,meta.motorpositions,meta.scalingfactors] = read_metadata_function(FileList_0(1).name);
% cd ..
for kkk = 1:numel(FileList)
    L{kkk} = imfinfo(FileList(kkk).name);
    meta.height = L{kkk}(1).Height;
    meta.width = L{kkk}(1).Width;
    meta.numberframes(kkk) = numel(L{kkk});
end
% initialize 3D matrix
binning = 1;

meta.framerate = meta.framerate/binning;
nb_frames_total = sum(floor(meta.numberframes/binning));
clear movie movie_p
useful_range_start = 1;

nb_planes = 1;
meta.framerate = meta.framerate/nb_planes;
nb_frames_perplane = floor(nb_frames_total/nb_planes);
for pp = 1:nb_planes
    counter_planes{pp} = 0;
    movie_p{pp} = zeros(meta.height,meta.width,nb_frames_perplane);
    for kkk = 1:numel(FileList)
        kkk
        % read raw data
        indicator = rem(sum(meta.numberframes(1:kkk-1)),nb_planes);
        if indicator < pp
            startingpoint = pp - indicator;
        else
            startingpoint = nb_planes + pp - indicator;
        end
        nb_frames_this_time = ceil((meta.numberframes(kkk)-startingpoint+1)/nb_planes);
        [movie_p{pp}(:,:,((counter_planes{pp}+1):(counter_planes{pp}+nb_frames_this_time))),movie_AVG_X{kkk,pp}] = read_movie(FileList(kkk).name,meta.width,meta.height,nb_frames_this_time,startingpoint,binning,L{kkk},nb_planes);
        counter_planes{pp} = counter_planes{pp} + nb_frames_this_time;
    end
    
    
    % discard einschwingvorgang
    movie_p{pp} = movie_p{pp}(:,:,useful_range_start:end);
    % align movie intra-trial
    reference{pp} = mean(movie_p{pp}(:,:,(-200:200) + round(nb_frames_total/2/nb_planes)),3);
    [movie_p{pp},offsety_resolved,offsetx_resolved] = alginWithinTrial(movie_p{pp},reference{pp});
    movie_AVG_X{pp} = mean(movie_p{pp},3);
    % check alignment visually
    % implay(movie_p{1}(:,:,1:20:end)/max(movie_p{1}(:)))
    % estimate offset of PMTs
    offset = quantile(movie_AVG_X{1}(:),0.02);
    % dF over F
end

% correct for DAQ board transient
% template = zeros(size(movie_p{pp}(:,:,1)));
% for pp = 1:nb_planes
%     dummy = mean(movie_p{pp}(:,:,1:end),3);
%     t1 =  nanmean(dummy(1:2:end,:),1);
%     t2 = nanmean(dummy(2:2:end,:),1);
%     template(1:2:end,:) = template(1:2:end,:) + repmat(  t1  , [size(template,2)/2 1]);
%     template(1:2:end,:) = template(1:2:end,:) + repmat( t2(end:-1:1)  , [size(template,2)/2 1]);
%     template(2:2:end,:) = template(2:2:end,:) +  repmat(t2 , [size(template,2)/2 1]);
%     template(2:2:end,:) = template(2:2:end,:) +  repmat( t1(end:-1:1), [size(template,2)/2 1]);
%     imagesc(template)
% end
% template = template/nb_planes/2;
% 
% for pp = 1:nb_planes
%     movie_pX{pp} = movie_p{pp} - repmat(template,[1 1 size(movie_p{pp},3)]);
% end

% LL = zeros(size(movie_AVG_X{pp}));
for pp = 1:nb_planes
    movie_AVG_X{pp} = mean(movie_p{pp}(:,:,200:800),3);
%     movie_AVG_X{pp} = movie_AVG_X{pp} - template;
% LL =  LL+ movie_AVG_X{pp};
end

for pp = 1:nb_planes
    plot1 = 34; plot2 = 0; DF_movie_yesno = 0; % figure number
    [DF_reponse{pp},DF_master{pp},DF_movie] = dFoverF(movie_p{pp},offset,meta.framerate,plot1,plot2,DF_movie_yesno);
    drawnow;
    % local correlation map (computational slightly expensive, but good
    % alternative to dF over F map; movie has to be 2^x for width and height
    tilesize = 16;
    localCorrelations{pp} = localCorrelationMap(movie_p{pp},tilesize);
    % localCorrelations = zeros(size(DF_reponse));
end

figure(88);
colormap(paruly)
for pp = 1:5
    anatomy{pp} = undistort_stack(mean(movie_p{pp}(:,:,1:end),3),meta.zoom);
    subplot(1,5,pp); imagesc(anatomy{pp},[8400 8700]); axis off equal; colormap(gray)
end

for pp = 1:5
   
    df_scale = [-10 200];
    last_ii = 1;
    % extract cellular time traces . semi-automated ROI-detection
    ROI_map_input = zeros(size(movie_p{pp}(:,:,1)));%   ROI_mapXX{pp};
    [ROI_mapX,timetracesX,timetracesX_raw] = timetraces_singleplane(movie_p{pp},movie_AVG_X{pp},offset,DF_reponse{pp},DF_master{pp},localCorrelations{pp},df_scale,ROI_map_input,meta,last_ii,movie_AVG_X{pp});
    % figure, imagesc(conv2(timetracesX,fspecial('gaussian',[25 1],23),'same'));
    ROI_map_input = ROI_mapX;

    ROI_mapXX{pp} = ROI_mapX;
    timetracesX_X{pp} = timetracesX;
    timetracesX_X_raw{pp} = timetracesX_raw;

    [size(timetracesX_X_raw{1},2), size(timetracesX_X_raw{2},2), size(timetracesX_X_raw{3},2),size(timetracesX_X_raw{4},2),size(timetracesX_X_raw{5},2)]
    sum([size(timetracesX_X_raw{1},2), size(timetracesX_X_raw{2},2), size(timetracesX_X_raw{3},2),size(timetracesX_X_raw{4},2),size(timetracesX_X_raw{5},2)])
end

figure(811);
colormap(lines)
for pp = 1:5
    ROI_mapXXZ{pp} = undistort_stack_integers(ROI_mapXX{pp},meta.zoom);
    subplot(5,1,pp); imagesc(ROI_mapXXZ{pp}); axis off equal; colormap(lines)
end

% fullX =     [(timetracesX_X_raw{1}), (timetracesX_X_raw{2}), (timetracesX_X_raw{3}),(timetracesX_X_raw{4}),(timetracesX_X_raw{5})]; 

DP_total.ROI_mapXX = ROI_mapXX;
DP_total.timetracesX_X = timetracesX_X;
DP_total.timetracesX_X_raw = timetracesX_X_raw;
DP_total.anatomy = anatomy;
DP_total.DF_reponse = DF_reponse;
DP_total.DF_master = DF_master;
DP_total.localCorrelations = localCorrelations;
save(strcat(FileList(1).name(1:end-4),'_timetrace.mat'),'DP_total','meta','offset');

