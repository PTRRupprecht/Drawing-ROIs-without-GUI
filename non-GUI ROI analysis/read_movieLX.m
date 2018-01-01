% read movie
function [movie,movie_AVG] = read_movieLX(filename,width,height,numberframes,startframe,binning,InfoImage,nb_planes)

% binning = 2;

movie = zeros(height,width,numberframes/binning);
for i = 1:numberframes/binning
    if mod(i,50) == 0; disp(strcat(num2str(i*binning),12,'frames read')); end
    if binning == 2
        movie(:,:,i) = (double(imread(filename,startframe-1 + binning*i-1,'Info',InfoImage)) + double(imread(filename,startframe-1 +binning*i,'Info',InfoImage)))/binning;
    else
        movie(:,:,i) = double(imread(filename,startframe+(i-1)*nb_planes,'Info',InfoImage));
    end
end
movie_AVG = mean(movie,3);

end