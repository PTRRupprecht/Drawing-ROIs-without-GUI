% controls the user interface with keyboard and mouse callback for
% analyzing a FOV as quickly as possible; check the
% switchImage()-subfunction for further details and a brief manual
%
% open window of FOV which allows to zoom and pan and select neurons
% manually ('f' for freehand) or automatically ('x' for
% crosscorrelation-based). Areasize can be chosen with the numbers 1-9 and
% 0; 'd' allows to delete a ROI; qwerty allow to switch between differents
% modes of presentation of the FOV.
function [ROI_mapX,timetracesX,timetracesX_raw] = timetraces_singleplane(movie,movie_AVG,offset,DF_reponse,DF_master,localCorrelations,df_scale,ROI_map_input,meta,last_ii,movie_AVG_X)
    global ROI_map neuron_index clut2b RGBimage timetraces timetracesX_raw display_state areaSize_userDefined replace_mode ROI_map_X avg_index
    width = meta.width; height = meta.height; framerate = meta.framerate;
    replace_mode = 0; 
    areaSize_userDefined = 0;
    nb_neurons = max(ROI_map_input(:)); % maixmally expected number of neurons
    timetraces = zeros(size(movie,3),nb_neurons);
    ROI_map = ROI_map_input;

    %% select neuron to be traced in all trials
    scaling = 1.0;
    mappe_X = ROI_map(:,:);
    mappe_X(mappe_X > 0) = 1;
    AVG_X = movie_AVG(:,:,1) - min(min(movie_AVG(:,:,1)));
    AVG_X = AVG_X/max(AVG_X(:));
%     AVG_X = imadjust(movie_AVG(:,:,1));
    AVG_X = adapthisteq(AVG_X,'NumTiles',[16 16]);
    RGBimage = zeros(size(movie_AVG,1),size(movie_AVG,2),3);
    RGBimage(:,:,1) = AVG_X*scaling; RGBimage(:,:,2) = AVG_X*scaling; RGBimage(:,:,3) = AVG_X*scaling;
    RGBimage(:,:,2) = RGBimage(:,:,2) + mappe_X/scaling; RGBimage(:,:,1) = RGBimage(:,:,1) + mappe_X/scaling;
    maxVal = max(max(RGBimage(:,:,1))); RGBimage(:,:,1) = RGBimage(:,:,1)/maxVal; RGBimage(:,:,2) = RGBimage(:,:,2)/max(max(RGBimage(:,:,2))); RGBimage(:,:,3) = RGBimage(:,:,3)/maxVal; 

    display_state = 'r';

    timetracesX_raw = zeros(size(movie,3),max(ROI_map(:)));
    for kk = 1:max(ROI_map(:))
        
        ROI_map_temp = ROI_map;
        ROI_map_temp(ROI_map_temp ~= kk) = 0;
        ROI_map_temp(ROI_map_temp>0) = 1;
        if sum(ROI_map_temp(:)) == 0
            timetraces(:,kk) = NaN;
            timetracesX_raw(:,kk) = NaN;
        else
            L = regionprops(ROI_map_temp,'BoundingBox');
            ROI_map_temp(ROI_map_temp == 0) = NaN;
            indizesx = ceil(L.BoundingBox(1)):floor((L.BoundingBox(1)+L.BoundingBox(3)));
            indizesy = ceil(L.BoundingBox(2)):floor((L.BoundingBox(2)+L.BoundingBox(3)));
            
            indizesx = unique(min(max(1,indizesx),size(movie,2)));
            indizesy = unique(min(max(1,indizesy),size(movie,1)));
             
            try
                movie_closeup = movie(indizesy,indizesx,:)  - offset ;
            catch
                keyboard
            end
            ROI_map_temp_closeup = ROI_map_temp(indizesy,indizesx);

            timetrace_X = squeeze(nanmean(nanmean(movie_closeup.*repmat(ROI_map_temp_closeup,[1 1 size(movie_closeup,3)]),1),2));
            timetracesX_rawX = timetrace_X;
            ftrace = smooth(timetrace_X,25); F0 = min(ftrace(25:end-25));
            timetrace_X = (timetrace_X - F0)/F0*100;
            timetraces(:,kk) = timetrace_X;
            timetracesX_raw(:,kk) = timetracesX_rawX;
        end
    end
    ROI_mapX = ROI_map;
    timetracesX = timetraces;
    
    try close 4121; end
    try close 4131; end
    try close 4233; end
    try close 5012; end
    try close 3331; end
    try close 777; end
    try close 8881; end
    
end







