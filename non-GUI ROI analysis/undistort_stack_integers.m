% undistort stack
function [movie_undistort,pixel_size] = undistort_stack(movie,zoom)
moviex = undistort(movie(:,:,1),zoom);
width_undistored = size(moviex,2);
width = size(movie,2);
[pixelsize_x, pixelsize_y] = pixelsize_xy(zoom,width,width_undistored);
x1 = 1:size(moviex,2);
x2 = linspace(1,size(moviex,2),size(moviex,2)*pixelsize_x/pixelsize_y);

movie_undistort = zeros(size(moviex,1),length(x2),size(movie,3));
for i = 1:size(movie,3)
    A = undistort_integers(movie(:,:,i),zoom);
    A2 = zeros(size(A,1),length(x2));
    for k = 1:size(A,1)
        A2(k,:) = interp1(x1,A(k,:),x2,'nearest');
    end
    movie_undistort(:,:,i) = A2;
end

pixel_size = pixelsize_y;
 