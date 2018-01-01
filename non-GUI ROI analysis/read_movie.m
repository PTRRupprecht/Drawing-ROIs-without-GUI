% read movie movie + AVG image
% implementation of the Tiff-read library by Anastasios Moressis (March 2016)
function [movie,movie_AVG] = read_movie(filename,width,height,numberframes,startframe,binning,InfoImage,nb_planes)
warning ('off','all');
% binning = 2;
TifLink = Tiff(filename, 'r');
movie = zeros(height,width,numberframes/binning,'uint16');
for i = 1:numberframes/binning
    if mod(i,50) == 0; disp(strcat(num2str(i*binning),12,'frames read')); end
    if binning == 2
        TifLink.setDirectory(startframe-1 + binning*i-1);
        movie(:,:,i) = TifLink.read();
        TifLink.setDirectory(startframe-1 +binning*i);
        movie(:,:,i) = (movie(:,:,i) + TifLink.read())/binning;
    else
        TifLink.setDirectory(startframe+(i-1)*nb_planes);
        movie(:,:,i) = TifLink.read();
    end
end 
movie = double(movie);
movie_AVG = mean(movie,3);
warning ('on','all');
end