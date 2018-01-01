function X = unwarp(XX)

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
    xpos_LUT = interp1(yy,x,1:pxperline_new);
    XX_reshaped = reshape(shiftdim(XX,1),size(XX,2),[]);
    X = interp1(1:pxperline,XX_reshaped,xpos_LUT,'linear');
    X = reshape(X,pxperline_new,size(XX,3),size(XX,1));
    X = shiftdim(X,2);
end

