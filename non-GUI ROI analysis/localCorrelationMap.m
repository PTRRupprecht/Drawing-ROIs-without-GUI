function Mppe = localCorrelationMap(movie,tilesize)

distMat = zeros(tilesize^2,tilesize^2);
for i = 1:tilesize^2
    for j = 1:tilesize^2
        [x1, y1] = ind2sub(size(distMat),i);
        [x2, y2] = ind2sub(size(distMat),j);
        distMat(i,j) = norm([x1 y1]-[x2 y2]);
    end
end

mask = 1./(distMat+1);
mask = mask - eye(size(mask)).*mask;

Mppe = zeros(size(movie,1),size(movie,2));

for i = 1:size(movie,1)/tilesize
    disp(strcat(num2str(i),12,'strip out of',12,num2str(size(movie,1)/tilesize)));
    for j = 1:size(movie,2)/tilesize
        junk = movie(((i-1)*tilesize:i*tilesize-1)+1,((j-1)*tilesize:j*tilesize-1)+1,:);
        junk = reshape(junk,[size(junk,1)*size(junk,2) size(junk,3)]);
        junk2 = corr(junk',junk');
        indizes = sum(junk2.*mask);
        yellow = reshape(indizes,[tilesize tilesize]);
        Mppe(((i-1)*tilesize:i*tilesize-1)+1,((j-1)*tilesize:j*tilesize-1)+1) =  Mppe(((i-1)*tilesize:i*tilesize-1)+1,((j-1)*tilesize:j*tilesize-1)+1) + yellow;

        if j < size(movie,2)/tilesize && i < size(movie,1)/tilesize
            junk = movie((((i-1)*tilesize+tilesize/2):(i*tilesize-1)+tilesize/2)+1,(((j-1)*tilesize+8):(j*tilesize-1)+8)+1,:);
            junk = reshape(junk,[size(junk,1)*size(junk,2) size(junk,3)]);
            junk2 = corr(junk',junk');
            indizes = sum(junk2.*mask);
            yellow = reshape(indizes,[tilesize tilesize]);
            Mppe((((i-1)*tilesize+tilesize/2):(i*tilesize-1)+tilesize/2)+1,(((j-1)*tilesize+tilesize/2):(j*tilesize-1)+tilesize/2)+1) =  Mppe((((i-1)*tilesize+tilesize/2):(i*tilesize-1)+tilesize/2)+1,(((j-1)*tilesize+tilesize/2):(j*tilesize-1)+tilesize/2)+1) + yellow;
        end
    
    end
end


end