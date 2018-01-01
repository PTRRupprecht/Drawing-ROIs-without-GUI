% Undistorting single images from resonant scanning
function Undistorted = undistort(Distorted,zoom)
    X = double(Distorted);

    imwidth = size(X,2);
    f_res = 7.91e3; % resonant scanning frequency

    dutyCycle = 12.5e-9*4096/(1/2/f_res); % percentage spent acquiring
    t_end = (1/2/f_res); 

    tx_start = (t_end - dutyCycle*t_end)/2; % start acquistion
    tx_end = t_end - (t_end - dutyCycle*t_end)/2; % end acquistion

    t = tx_start:1e-6:tx_end; % acquistion time window
    v = sin(t*2*pi*f_res);  % scanning velocity

    full_x = (cos(2*pi*f_res*tx_start) - cos(2*pi*f_res*tx_end))/(2*pi*f_res);  % spatial distance after scanning
    scalingFactor = ceil(imwidth/(max(v)/mean(v)))/full_x;  % spatial distance corresponds to xyz pixels
    v = scalingFactor*sin(t*2*pi*f_res); % scale velocity from a.u. to pixels

    % figure, plot(t,v); ylabel('velocity in px per sec');
    % xlabel('time in sec')

    % FF is the final image, FF_counts counts how often something has been
    % added to FF, allowing for subsequent normalisation;
    % interpolation is simply done by splitting the calculated pixel position
    % in the two neighoring pixels, with a ratio according to the tendence of
    % the un-rounded pixel value
    FF = zeros(size(X,1),1+ceil(imwidth/(max(v)/mean(v)))); FF_count = zeros(size(X,1),1+ceil(imwidth/(max(v)/mean(v))));
    for i = 1:imwidth
        tx_now = tx_start + (i)*12.5e-9*4096/imwidth;
        pixel_value = scalingFactor*(cos(2*pi*f_res*tx_start) - cos(2*pi*f_res*tx_now))/(2*pi*f_res);
        y_floor = max(1,floor(pixel_value));
        floor_perc =pixel_value - y_floor;
        FF(:,y_floor) = FF(:,y_floor) + X(:,i)*(1-floor_perc);
        FF(:,y_floor+1) = FF(:,y_floor+1) + X(:,i)*(floor_perc);
        FF_count(:,y_floor) = FF_count(:,y_floor) + 1-floor_perc;
        FF_count(:,y_floor+1) = FF_count(:,y_floor+1) + (floor_perc);
    end
    Undistorted = FF./FF_count; % final image

end