function X = unwarp_precision(XX)

    pxperline = size(XX,2); % pixel per line
    xy_conversionfactor = 477/512; % empirically found through calibration for my setup
    pxperline_new = pxperline*xy_conversionfactor;
    f_res = 7.91e3; % resonant scanning frequency
    linetime = 1/f_res/2;

    linedwelltime = linetime*4096/(80e6/f_res/2); % in toto
    t_offset = (linetime - linedwelltime)/2;

    x = 1:pxperline;
    t = t_offset + x/pxperline*linedwelltime;
    y  = cos(t*2*pi*f_res);

    yy = y - min(y); yy = 1+(yy/max(yy)*(pxperline_new-1));
    yy = yy(end:-1:1);
    yy_floor = floor(yy);
    yy_remainder = yy - yy_floor; % which fraction goes to the left/right pixel?
    
    X = zeros(size(XX,1),pxperline_new,size(XX,3));
    histX = zeros(size(XX,1),pxperline_new);
    for k = 1:numel(yy_floor) % split each pixel into the two bins of the new grid
        X(:,yy_floor(k),:) = X(:,yy_floor(k),:) + XX(:,k,:)*(1-yy_remainder(k));
        histX(:,yy_floor(k)) = histX(:,yy_floor(k)) + (1-yy_remainder(k));
        if k < numel(yy_floor)
            X(:,yy_floor(k)+1,:) = X(:,yy_floor(k)+1,:) + XX(:,k+1,:)*yy_remainder(k);
            histX(:,yy_floor(k)+1) = histX(:,yy_floor(k)+1) + yy_remainder(k);
        end
    end
    X = X./repmat(histX,[1 1 size(X,3)]); % normalize
end

