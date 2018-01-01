% open window of FOV which allows to zoom and pan and select neurons
% manually ('f' for freehand) or automatically ('x' for
% crosscorrelation-based). Areasize can be chosen with the numbers 1-9 and
% 0; 'd' allows to delete a ROI; qwerty allow to switch between differents
% modes of presentation of the FOV.
function switchImage(~,event,AVG_X,DF_reponse,DF_master,localCorrelations,df_scale,movie,offset,handle1,framerate,zoom,last_ii,movie_AVG_X)
 
    global clut2b ROI_map neuron_index RGBimage timetraces display_state replace_mode areaSize_userDefined ROI_map_X timetracesX_X avg_index
    clear keyword
    keyword = event.Key;
    if length(keyword) == 1 && any(keyword == ['q' 'w' 'e' 'r' 't' 'y' 'u' 'n' 'm'])
        keyword = 'display';
    end
    if length(keyword) == 1 && any(keyword == ['1' '2' '3' '4' '5' '6' '7' '8' '9' '0'])
        keyword = 'changeAreaSize';
    end
    switch(keyword)
        % change user defined Area; if '0' is chosen, it will be set to the
        % default value, which is the average ROI area after n=10 ROIs have
        % been selected
        case 'changeAreaSize'
            areaSize_userDefined = str2double(event.Key)*25;
        % change view of FOV window
        case 'display'
            keyword2 = event.Key;
            change_display(keyword2,AVG_X,DF_reponse,df_scale,handle1,DF_master,localCorrelations,last_ii,movie_AVG_X)
        case 'd'
            x_x = ginput(1); x_x = round(x_x);
            delete_ROI(x_x);
            change_display(display_state,AVG_X,DF_reponse,df_scale,handle1,DF_master,localCorrelations,last_ii,movie_AVG_X)
        % draw ROI manually
        case 'f'
            h = imfreehand(gca);
            mappe = createMask(h); mappe = double(mappe);
            delete(h);
            mappe(mappe == 0) = nan;

            % determine selected position coordinates in zoom-in
            xy = find( mappe == 1 );
            [x, y] = ind2sub(size(mappe),xy);
            xx = round(mean(x)); yy = round(mean(y));

            replace_ROI([yy xx]);

            tsize = 35; 
            xindizes = min(max(xx-tsize:xx+tsize,1),size(movie,1));
            yindizes = min(max(yy-tsize:yy+tsize,1),size(movie,2));
            % zoom-in
            movie_closeup = movie(unique(xindizes),unique(yindizes),:)  - offset ;
            map_closeup = mappe(unique(xindizes),unique(yindizes));

            timetrace_X = squeeze(nanmean(nanmean(movie_closeup.*repmat(map_closeup,[1 1 size(movie,3)]),1),2));
            ftrace = smooth(timetrace_X,25); F0 = min(ftrace(25:end-25));
            timetrace_X = (timetrace_X - F0)/F0*100;

            if replace_mode == 0; index_cell = neuron_index; else index_cell = replace_mode; end
            ROI_map(xy) = index_cell;
            timetraces(:,index_cell) = timetrace_X;

            figure(4131); plot((1:length(timetrace_X))/framerate,timetrace_X); xlabel('time (seconds)');
            figure(4198);
            
            increase_neuroncounter();

            change_display(display_state,AVG_X,DF_reponse,df_scale,handle1,DF_master,localCorrelations,last_ii,movie_AVG_X)
            
        case 'g'
            [yyy,xxx]=ginput(1);
            xx = round(xxx); yy = round(yyy);

            colorindex = ROI_map(xx,yy);
            if colorindex > 0
                mappe = double((ROI_map == colorindex));
                timetraces(:,colorindex) = NaN;
                replace_mode = colorindex;
            end
            
            mappe(mappe == 0) = nan;

            % determine selected position coordinates in zoom-in
            xy = find( mappe == 1 );
            [x, y] = ind2sub(size(mappe),xy);
            xx2 = round(mean(x)); yy2 = round(mean(y));
            mappe2 = zeros(size(mappe));
            for j = 1:length(x)
                mappe2(min(size(mappe2,1),max(1,x(j) + xx - xx2)),min(size(mappe2,2),max(1,y(j) + yy - yy2))) = 1;
            end
            clear xy
            xy = find( mappe2 == 1 );

            tsize = 35; 
            xindizes = min(max(xx-tsize:xx+tsize,1),size(movie,1));
            yindizes = min(max(yy-tsize:yy+tsize,1),size(movie,2));
            % zoom-in
            movie_closeup = movie(unique(xindizes),unique(yindizes),:)  - offset ;
            map_closeup = mappe2(unique(xindizes),unique(yindizes));

            timetrace_X = squeeze(nanmean(nanmean(movie_closeup.*repmat(map_closeup,[1 1 size(movie,3)]),1),2));
            ftrace = smooth(timetrace_X,25); F0 = min(ftrace(25:end-25));
            timetrace_X = (timetrace_X - F0)/F0*100;

            if replace_mode == 0; index_cell = neuron_index; else index_cell = replace_mode; end
            ROI_map(ROI_map == index_cell) = 0;
            ROI_map(xy) = index_cell;
            timetraces(:,index_cell) = timetrace_X;

            figure(4131); plot((1:length(timetrace_X))/framerate,timetrace_X); xlabel('time (seconds)');
            figure(4198);
            
            increase_neuroncounter();

            change_display(display_state,AVG_X,DF_reponse,df_scale,handle1,DF_master,localCorrelations,last_ii,movie_AVG_X)
            
        % correlation-based point-seeded ROI
        case 'x'
            [yyy,xxx]=ginput(1);
            xx = round(xxx); yy = round(yyy);
            
            replace_ROI([yy xx]);
            
            tsize = 25; 
            [timetrace_X, reference_trace] = extract_timetraceAndROI(xx,yy,tsize,movie,offset);

            figure(4131); plot((1:length(reference_trace))/framerate,reference_trace); hold on; plot((1:length(timetrace_X))/framerate,timetrace_X); hold off; xlabel('time (seconds)');
            figure(4198);
            
            increase_neuroncounter();
            change_display(display_state,AVG_X,DF_reponse,df_scale,handle1,DF_master,localCorrelations,last_ii,movie_AVG_X)
        case 'v'
            x_x = ginput(1); x_x = round(x_x);
            colorindex = ROI_map(x_x(2),x_x(1));

            timetrace_X = timetraces(:,colorindex);
            
%             [yyy,xxx]=ginput(1);
%             xx = round(xxx); yy = round(yyy);
%             
%             tsize = 25; 
%             [timetrace_X, reference_trace] = extract_timetraceAndROI(xx,yy,tsize,movie,offset);
%             drawnow;
            
            xcorr_tot = zeros(size(movie,1),size(movie,2));
            for j = 1:size(movie,1)
                if mod(j,50) == 0; disp(j/size(movie,1)); end
                xcorr_tot(j,:) = corr(timetrace_X,squeeze(movie(j,:,:))');
            end
            figure(randi(2000,1)+99999); imagesc(xcorr_tot); colormap(gray);
            
            change_display(display_state,AVG_X,DF_reponse,df_scale,handle1,DF_master,localCorrelations,last_ii,movie_AVG_X)
        case 'z'
            [yyy,xxx]=ginput(1);
            xx = round(xxx); yy = round(yyy);
            
            replace_ROI([yy xx]);
            
            tsize = 35; 
            xindizes = min(max(xx-tsize:xx+tsize,1),size(movie,1));
            yindizes = min(max(yy-tsize:yy+tsize,1),size(movie,2));
            % zoom-in
            movie_closeup = movie(unique(xindizes),unique(yindizes),:)  - offset ;
            map_circle = zeros(size(movie_closeup,1),size(movie_closeup,2));
            centerx = tsize + 1 - max(sum(xindizes == 1),1) + 1; % accounts for border effects
            centery = tsize + 1 - max(sum(yindizes == 1),1) + 1;
            map_circle(centerx,centery) = 1; map_circle = bwdist(map_circle);
            [~,I] = sort(map_circle(:),'ascend'); map_circle(:) = 0;
            if areaSize_userDefined ~= 0
                AreaSize = areaSize_userDefined;         
            elseif length(unique(ROI_map)) > 10
                    ROI_map_thresh = ROI_map; ROI_map_thresh(ROI_map_thresh > 1) = 1;
                    AreaSize = 0.8*sum(ROI_map_thresh(:))/(length(unique(ROI_map))-1);
            else
                    AreaSize = 150;
            end
            map_circle(I(1:round(AreaSize))) = 1.0;
            ROI_map_closeup = ROI_map(unique(xindizes),unique(yindizes));
            ROI_map_closeup(ROI_map_closeup>1) = 1; ROI_map_closeup = abs(1-ROI_map_closeup);
            map_circle = map_circle.*ROI_map_closeup;
            timetrace_X = squeeze(nanmean(nanmean(movie_closeup.*repmat(map_circle,[1 1 size(movie,3)]),1),2));
            ftrace = smooth(timetrace_X,25); F0 = min(ftrace(25:end-25));
            timetrace_X = (timetrace_X - F0)/F0*100;

            if replace_mode == 0; index_cell = neuron_index; else index_cell = replace_mode; end
            ROI_map(unique(xindizes),unique(yindizes)) = ROI_map(unique(xindizes),unique(yindizes)) + map_circle*index_cell;
            timetraces(:,index_cell) = timetrace_X;
            
            figure(4131); plot((1:length(timetrace_X))/framerate,timetrace_X); hold off; xlabel('time (seconds)');
            figure(4198);
            
            increase_neuroncounter();
            change_display(display_state,AVG_X,DF_reponse,df_scale,handle1,DF_master,localCorrelations,last_ii,movie_AVG_X)
        
        case 'j'
            % without clustermap
            online_analysis_overview(ROI_map,timetraces,framerate,0,3,zoom);
        case 'k'
            % with clustermap
            online_analysis_overview(ROI_map,timetraces,framerate,1,3,zoom);
        case 'l'
            % with time-dependent correlation matrizes
            online_analysis_overview(ROI_map,timetraces,framerate,2,3,zoom);
        case 'h'
            % with saving data to disk
            online_analysis_overview(ROI_map,timetraces,framerate,3,3,zoom);
        case 'p'
            % single-cell trial-to-trial variations
            x_x = ginput(1); x_x = round(x_x);
            colorindex = ROI_map(x_x(2),x_x(1));
            figure(5012);
            plot((1:size(timetraces,1))/framerate,smooth(timetraces(:,colorindex),1),'k');
            hold on;
%             lll = length(timetracesX_X);
%             for ii = 1:lll
%                 try
%                     if nansum(timetracesX_X{ii}(:,colorindex)) ~= 0
%                         plot((1:size(timetracesX_X{ii},1))/framerate,smooth(timetracesX_X{ii}(:,colorindex),1),'Color',[ii/lll 1-ii/lll 1-ii/lll]);
%                     end
%                 end
%             end
            xlabel('time [sec]'); ylabel('dF/F [%]');
            hold off;
            figure(4198)
    end
end

% online analysis (i.e. during data analysis to get an overview or insight)
function online_analysis_overview(ROI_map,timetraces,framerate,operation_mode,nb_clusters,zoom)
    
    detect_traces = find(nansum(timetraces));
    dF = timetraces(:,detect_traces);            
    
    corrMatrix = corr(dF);
    corrMatrixOld = corrMatrix;
    % remove diagonal elements
    corrMatrix = corrMatrix - eye(size(corrMatrix));
    % and convert to a vector (as pdist)
    dissimilarity = 1 - corrMatrix';

    % perform complete linkage clustering
    Z = linkage(dissimilarity,'complete');
    % decide on a cutoff
    cutoff = 10;
    l_groups = 1;
    while l_groups < nb_clusters
        cutoff = cutoff - 0.1;
        groups = cluster(Z,'cutoff',cutoff,'criterion','distance');
        l_groups = length(unique(groups));
    end
    
    dF = timetraces(:,detect_traces);
    [~,IX] = sort(groups);
    dF_ordered = dF(:,IX);
    
    corrMatrix2 = corr(dF_ordered);
    
    [ROI_mapX_undistort,pixelsize] = undistort_stack(ROI_map,zoom);
%     ROI_mapX_undistort = round(ROI_mapX_undistort);
    ROI_mapX_undistort(mod(ROI_mapX_undistort,1)~=0)=0;
    
    
    center = zeros(2,length(detect_traces));
    for i = 1:length(detect_traces)
        [x,y] = find(ROI_mapX_undistort == detect_traces(i));
        m_x = mean(x); m_y = mean(y);
        center(:,i) = [m_x, m_y];
    end
    
    distMap = pixelsize*squareform(pdist(center(:,IX)'));
    distMapOld = pixelsize*squareform(pdist(center(:,:)'));

    corrMatrix_states = corr(dF');
    
    time_t = (1:size(timetraces,1))/framerate;
    
    figure(8881);
    subplot(2,2,1); imagesc(detect_traces,detect_traces,corrMatrix2); xlabel('neuron index'); ylabel('neuron index');
    subplot(2,2,2); imagesc(time_t,detect_traces,dF_ordered'); xlabel('time [sec]'); ylabel('neuron index');
    
    [values,axes] = hist3([distMap(:) corrMatrix2(:)],[100 100]);
    values = conv2(values,fspecial('gaussian',[5 5],3.5),'same');
    subplot(2,2,3);  imagesc(axes{2},axes{1},log(values))
    %scatter(distMap(:),corrMatrix2(:)); xlabel('distance [px]'); ylabel('correlation of timetraces'); box on;
    subplot(2,2,4); imagesc(time_t,time_t,corrMatrix_states); xlabel('time [sec]'); ylabel('time [sec]');
    
    if operation_mode == 2
        figure(777),
        for i = 1:6
            delta = floor(size(timetraces,1)/6);
            corrMatrix = corr(dF_ordered((i-1)*delta+1:i*delta,:));
            subplot(2,3,i); imagesc(detect_traces,detect_traces,corrMatrix,[-0.2 1]);
            xlabel('neuron index'); ylabel('neuron index');
        end
    elseif operation_mode == 1
        cluster_map = ROI_map;
        for k = 1:length(detect_traces)
            cluster_map(ROI_map == detect_traces(k)) = groups(k);
        end
        figure(3331); imagesc(cluster_map);
    elseif operation_mode == 3
        save('distance_X.mat','distMapOld','corrMatrixOld');
    end
    figure(4198)
end

% keep track of which cell number is replaced/added
function increase_neuroncounter()
    global neuron_index timetraces replace_mode
    if mod(neuron_index,100) == 0 && replace_mode == 0
        temp = timetraces;
        timetraces = zeros(size(timetraces,1),size(timetraces,2)+100);
        timetraces(:,1:end-100) = temp;
        clear temp
    end
    if replace_mode ~= 0
        replace_mode = 0;
    else
        neuron_index = neuron_index + 1;
    end
end

% delete ROI which contains the seeded point
function delete_ROI(x_x)
global timetraces ROI_map
    colorindex = ROI_map(x_x(2),x_x(1));
    if colorindex > 0
        ROI_map(ROI_map == colorindex) = 0;
        timetraces(:,colorindex) = NaN;
    end
    drawnow;
end

% replace ROI which contains the seeded point
function replace_ROI(x_x)
global timetraces ROI_map replace_mode
    colorindex = ROI_map(x_x(2),x_x(1));
    if colorindex > 0
        ROI_map(ROI_map == colorindex) = 0;
        timetraces(:,colorindex) = NaN;
        replace_mode = colorindex;
    end
    drawnow;
end

% update main display of FOV
function change_display(keyword2,AVG_X,DF_reponse,df_scale,handle1,DF_master,localCorrelations,last_ii,movie_AVG_X)
global display_state clut2b RGBimage ROI_map ROI_map_X avg_index
    switch keyword2
        case 'q'
             set(handle1,'CData',AVG_X); colormap(gray(256)); caxis auto
             display_state = 'q';
             avg_index = 0;
        case 'w'
              set(handle1,'CData',100*DF_reponse); caxis(df_scale); colormap(clut2b);
              display_state = 'w';
              avg_index = 0;
        case 'e'
             set(handle1,'CData',100*DF_master); caxis(df_scale); colormap(clut2b);
             display_state = 'e';
             avg_index = 0;
        case 'r'
            scaling = 2.5;
            ROI_map_thresh = ROI_map; ROI_map_thresh(ROI_map_thresh>1) = 1;
            RGBimage(:,:,1) = AVG_X*scaling; RGBimage(:,:,2) = AVG_X*scaling; RGBimage(:,:,3) = AVG_X*scaling;
            RGBimage(:,:,2) = RGBimage(:,:,2) + ROI_map_thresh/scaling; RGBimage(:,:,1) = RGBimage(:,:,1) + ROI_map_thresh/scaling;
            maxVal = max(max(RGBimage(:,:,1))); RGBimage(:,:,1) = RGBimage(:,:,1)/maxVal; RGBimage(:,:,2) = RGBimage(:,:,2)/max(max(RGBimage(:,:,2))); RGBimage(:,:,3) = RGBimage(:,:,3)/maxVal; 

            set(handle1,'CData',RGBimage); caxis auto
            display_state = 'r';
            avg_index = 0;
        case 't'
             cmap = lines(1000); cmap(1,:) = 1;
             set(handle1,'CData',ROI_map); colormap(cmap); caxis auto
             display_state = 't';
             avg_index = 0;
        case 'y'
            localCorrelationsX = conv2(localCorrelations,fspecial('gaussian',[2 2],1),'same');
            set(handle1,'CData',localCorrelationsX); colormap(jet(256)); caxis auto;
            display_state = 'y';
            avg_index = 0;
        case 'u'
            if isempty(ROI_map_X)
                change_display('r',AVG_X,DF_reponse,df_scale,handle1,DF_master,localCorrelations,last_ii,movie_AVG_X)
            else
                AVG_X2 = movie_AVG_X{last_ii}/max(max(movie_AVG_X{last_ii}));
                AVG_X2 = imadjust(AVG_X2);
                AVG_X2 = adapthisteq(AVG_X2,'NumTiles',[16 16]);
                scaling = 2.5;
                ROI_map_thresh = ROI_map_X{last_ii}; ROI_map_thresh(ROI_map_thresh>1) = 1;
                RGBimage(:,:,1) = AVG_X2*scaling; RGBimage(:,:,2) = AVG_X2*scaling; RGBimage(:,:,3) = AVG_X2*scaling;
                RGBimage(:,:,2) = RGBimage(:,:,2) + ROI_map_thresh/scaling; RGBimage(:,:,1) = RGBimage(:,:,1) + ROI_map_thresh/scaling;
                maxVal = max(max(RGBimage(:,:,1))); RGBimage(:,:,1) = RGBimage(:,:,1)/maxVal; RGBimage(:,:,2) = RGBimage(:,:,2)/max(max(RGBimage(:,:,2))); RGBimage(:,:,3) = RGBimage(:,:,3)/maxVal; 

                set(handle1,'CData',RGBimage); caxis auto
                display_state = 'u';
            end
            avg_index = 0;
        case 'n'
            avg_index = max(-last_ii,avg_index - 1);
            round_index = min(max(1,last_ii + avg_index),numel(movie_AVG_X));
            AVG_X2 = movie_AVG_X{round_index}/max(max(movie_AVG_X{round_index}));
            AVG_X2 = imadjust(AVG_X2);
            AVG_X2 = adapthisteq(AVG_X2,'NumTiles',[16 16]);
            scaling = 2.5;
            try
                ROI_map_thresh = ROI_map_X{round_index}; ROI_map_thresh(ROI_map_thresh>1) = 1;
                if isempty(ROI_map_X{round_index}); ROI_map_thresh = zeros(size(AVG_X2)); end
            catch
                ROI_map_thresh = ROI_map; ROI_map_thresh(ROI_map_thresh>1) = 1;
            end
            RGBimage(:,:,1) = AVG_X2*scaling; RGBimage(:,:,2) = AVG_X2*scaling; RGBimage(:,:,3) = AVG_X2*scaling;
            RGBimage(:,:,2) = RGBimage(:,:,2) + ROI_map_thresh/scaling; RGBimage(:,:,1) = RGBimage(:,:,1) + ROI_map_thresh/scaling;
            maxVal = max(max(RGBimage(:,:,1))); RGBimage(:,:,1) = RGBimage(:,:,1)/maxVal; RGBimage(:,:,2) = RGBimage(:,:,2)/max(max(RGBimage(:,:,2))); RGBimage(:,:,3) = RGBimage(:,:,3)/maxVal; 
            set(handle1,'CData',RGBimage); caxis auto
            display_state = 'u';
        case 'm'
            avg_index = min(avg_index + 1,numel(movie_AVG_X));
            round_index = min(max(1,last_ii + avg_index),numel(movie_AVG_X));
            AVG_X2 = movie_AVG_X{round_index}/max(max(movie_AVG_X{round_index}));
            AVG_X2 = imadjust(AVG_X2);
            AVG_X2 = adapthisteq(AVG_X2,'NumTiles',[16 16]);
            scaling = 2.5;
            try
                ROI_map_thresh = ROI_map_X{round_index}; ROI_map_thresh(ROI_map_thresh>1) = 1;
                if isempty(ROI_map_X{round_index}); ROI_map_thresh = zeros(size(AVG_X2)); end
            catch
                ROI_map_thresh = ROI_map; ROI_map_thresh(ROI_map_thresh>1) = 1;
            end
            RGBimage(:,:,1) = AVG_X2*scaling; RGBimage(:,:,2) = AVG_X2*scaling; RGBimage(:,:,3) = AVG_X2*scaling;
            RGBimage(:,:,2) = RGBimage(:,:,2) + ROI_map_thresh/scaling; RGBimage(:,:,1) = RGBimage(:,:,1) + ROI_map_thresh/scaling;
            maxVal = max(max(RGBimage(:,:,1))); RGBimage(:,:,1) = RGBimage(:,:,1)/maxVal; RGBimage(:,:,2) = RGBimage(:,:,2)/max(max(RGBimage(:,:,2))); RGBimage(:,:,3) = RGBimage(:,:,3)/maxVal; 
            set(handle1,'CData',RGBimage); caxis auto
            display_state = 'u';
    end
end

% calculates correlation-based ROI based on seed point and some additional spatial
% filtering; no overlap with other ROIs as another boundary condition
function [timetrace_X, reference_trace, mappe2, movie_closeup_AVGeq,centerx,centery] = extract_timetraceAndROI(xx,yy,tsize,movie,offset)

    global areaSize_userDefined ROI_map replace_mode timetraces neuron_index
    xindizes = min(max(xx-tsize:xx+tsize,1),size(movie,1));
    yindizes = min(max(yy-tsize:yy+tsize,1),size(movie,2));
    centerx = tsize + 1 - max(sum(xindizes == 1),1) + 1; % accounts for border effects
    centery = tsize + 1 - max(sum(yindizes == 1),1) + 1;

    % zoom-in
    movie_closeup = movie(unique(xindizes),unique(yindizes),:)  - offset ;
    movie_closeup_AVG = mean(movie_closeup,3);
    movie_closeup_AVGeq = adapthisteq(movie_closeup_AVG/max(movie_closeup_AVG(:)),'NumTiles',[4 4]);    
    sizetx = size(movie_closeup_AVG,1);
    sizety = size(movie_closeup_AVG,2);

    ROI_map_closeup = ROI_map(unique(xindizes),unique(yindizes));
    ROI_map_closeup(ROI_map_closeup>1) = 1; ROI_map_closeup = abs(1-ROI_map_closeup);

    reference_trace = squeeze(mean(mean(movie_closeup([-3:3]+centerx,[-3:3]+centery,:),1),2));
    tile_2D = reshape(movie_closeup,[sizetx*sizety size(movie_closeup,3)]);
    Mapp = corr(tile_2D',reference_trace);
    Mapp2D = reshape(Mapp,[sizetx sizety]).*ROI_map_closeup;

    [~,I] = sort(Mapp2D(:),'descend');

    if areaSize_userDefined ~= 0
        AreaSize = areaSize_userDefined;         
    elseif length(unique(ROI_map)) > 10
        ROI_map_thresh = ROI_map; ROI_map_thresh(ROI_map_thresh > 1) = 1;
        AreaSize = 0.8*sum(ROI_map_thresh(:))/(length(unique(ROI_map))-1);
    else
        AreaSize = 150;
    end
    mappe = zeros(sizetx,sizety);
    mappe(I(1:round(AreaSize))) = 1.0;
    mappe = imclose(mappe,strel('disk',1));
    cc = bwconncomp(mappe,4);
    S = regionprops(cc,'Area','PixelList','PixelIdxList');
    centerSpot = [centerx, centery];
    for i = 1:length(S)
        if ~(any(sum(abs(bsxfun(@minus,S(i).PixelList,centerSpot)),2) == 0))
           mappe(S(i).PixelIdxList) = 0;
        end
    end
    mappe = double(bwfill(mappe,'holes'));
    distmap = bwdist(mappe);
    weighted_environment = Mapp2D./distmap;
    [~,I] = sort(weighted_environment(:),'descend');
    mappe2 = zeros(sizetx,sizety);
    mappe2(I(1:(round(AreaSize*1.5)))) = 1.0;
    mappe2 = double(bwfill(mappe2,'holes'));
    cc = bwconncomp(mappe2,4);
    S = regionprops(cc,'Area','PixelList','PixelIdxList');
    centerSpot = [centerx, centery];
    for i = 1:length(S)
        if ~(any(sum(abs(bsxfun(@minus,S(i).PixelList,centerSpot)),2) == 0))
           mappe2(S(i).PixelIdxList) = 0;
        end
    end
    mappe2 = double(bwfill(mappe2,'holes')).*ROI_map_closeup;

    mappe2(mappe2 == 0) = nan;

    ftrace = smooth(reference_trace,25); F0 = min(ftrace(25:end-25));
    reference_trace = (reference_trace - F0)/F0*100;

    timetrace_X = squeeze(nanmean(nanmean(movie_closeup.*repmat(mappe2,[1 1 size(movie_closeup,3)]),1),2));
    ftrace = smooth(timetrace_X,25); F0 = min(ftrace(25:end-25));
    timetrace_X = (timetrace_X - F0)/F0*100;
    
    scaling = 1.5;
    mappeX = mappe2; mappeX(isnan(mappeX)) = 0;
    clear RGBim; RGBim(:,:,1) = movie_closeup_AVGeq*scaling; RGBim(:,:,2) = movie_closeup_AVGeq*scaling; RGBim(:,:,3) = movie_closeup_AVGeq*scaling;
    RGBim(:,:,2) = RGBim(:,:,2) + mappeX/scaling; RGBim(:,:,1) = RGBim(:,:,1) + mappeX/scaling;
    maxVal = max(max(RGBim(:,:,1))); RGBim(:,:,1) = RGBim(:,:,1)/maxVal; RGBim(:,:,2) = RGBim(:,:,2)/max(max(RGBim(:,:,2))); RGBim(:,:,3) = RGBim(:,:,3)/maxVal; 
    
    figure(4121); axes_h.a = subplot(1,2,1,'replace'); imagesc(RGBim); colormap(pink);
    hold on; plot(centery,centerx,'xy','MarkerSize',6); hold off; xlabel('Current trial');
    axes_h.b = subplot(1,2,2,'replace'); imagesc(conv2(Mapp2D,fspecial('gaussian',[2 2],1.5),'same'));
    set(gcf, 'units','normalized','Position',[0.1 0.1 0.6 0.4]);
    
    
    if replace_mode == 0; index_cell = neuron_index; else index_cell = replace_mode; end
    timetraces(:,index_cell) = timetrace_X;
    ROI_map(unique(xindizes),unique(yindizes)) = ROI_map(unique(xindizes),unique(yindizes)) + mappeX*index_cell;
end