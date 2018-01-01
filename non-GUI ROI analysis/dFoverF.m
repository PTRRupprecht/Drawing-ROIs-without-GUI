% calculate dFoverF for a response window and for a moving window, taking
% the maximum response for each location
function [DF_reponse,DF_master,DF_movie] = dFoverF(movie,offset,framerate,plot1,plot2,DF_movie_yesno,f0_window,response_window)

% figure(66), plot(squeeze(mean(mean(movie,2),1)))
% xy = ginput(4);
xy = [ f0_window(1) 0; f0_window(2) 0; response_window(1) 0; response_window(2) 0];

F0_window = round(xy(1,1)):round(xy(2,1));
response_window = round(xy(3,1)):round(xy(4,1));
string_window_F0 = strcat('F0 window between frame',12,num2str(round(xy(1,1))),12,'(=',12,num2str(round(xy(1,1))/framerate),12,'sec) and frame',12,num2str(round(xy(2,1))),12,'(=',12,num2str(round(xy(2,1))/framerate),12,'sec).');
string_windows_response = strcat('Response window between frame',12,num2str(round(xy(3,1))),12,'(=',12,num2str(round(xy(3,1))/framerate),12,'sec) and frame',12,num2str(round(xy(4,1))),12,'(=',12,num2str(round(xy(4,1))/framerate),12,'sec).');
disp(string_window_F0);
disp(string_windows_response);
% close 66

DF = mean(movie(:,:,response_window),3);

% F0Z = zeros(size(movie,1),size(movie,2),round(size(movie,3)/200)-1);
% for i = 1:round(size(movie,3)/200)-1
%     F0Z(:,:,i) = mean(movie(:,:,((i-1):i)*200+1),3);
% end
% F0Z = smooth3(F0Z,'gaussian',[1 1 1],2);
% F0 = min(F0Z,[],3);
    
F0 = mean(movie(:,:,F0_window),3);
F0 = conv2(F0,fspecial('gaussian',[3 3], 1),'same');
DF = (DF - F0)./(F0-offset);
DF_reponse = conv2(DF,fspecial('gaussian',[3 3], 2),'same');
if plot1 ~=0; figure(plot1), imagesc(DF_reponse*100,[-5 130]); colormap(jet); axis equal off; end
drawnow;
DF_master = zeros(size(movie,1),size(movie,2));
sliding_window_size = 50;
for i = 1:floor(size(movie,3)/sliding_window_size)
    DF = mean(movie(:,:,((i-1)*sliding_window_size+1):(min(i*sliding_window_size+1,size(movie,3)))),3);
%     F0 = mean(movie(:,:,F0_window),3);
%     F0 = conv2(F0,fspecial('gaussian',[5 5], 4),'same');
    DF = (DF - F0)./(F0-offset*1.00);
    DF = conv2(DF,fspecial('gaussian',[3 3], 2),'same');
    DF_master = max(DF_master,DF);
%     figure(3); imagesc(DF_master,[0 200]);
end
if plot2 ~= 0; figure(plot2), imagesc(DF_master*100,[-5 130]); colormap(jet); axis equal off; end
 
if DF_movie_yesno
    F0 = mean(movie(:,:,F0_window),3);
    F0 = conv2(F0,fspecial('gaussian',[3 5], 3),'same');
    DF_movie = (movie - repmat(F0,[1 1 size(movie,3)]))./(repmat(F0,[1 1 size(movie,3)]) - offset);
else
    DF_movie = [];
end