% Bidirectionally aligns 3D movie
% Input: movie, an x*y*t 3D array
% Output: movieL, an x*y*t 3D array

function movieL = bidi_align(movie)

    % subpixel computation of the shift for each image of the stack
    shifts = zeros(size(movie,3),1);
    cross_template = zeros(size(movie,2),1)';
    cross_template2 = zeros(size(movie,2),1)';
    for j = 1:size(movie,3)
        j/size(movie,3)
        for o = 1:size(movie,1)/2
            dd1 = movie((o-1)*2+1,:,j);
            dd2 = movie((o-1)*2+2,:,j);
            cross_template = cross_template + ifftshift(ifft(conj(fft(dd1)).*fft(dd1)));
            cross_template2 = cross_template2 + ifftshift(ifft(conj(fft(dd1)).*fft(dd2)));
        end
        [~,C0max] = max(cross_template);
        subpixelCorrection = (log(cross_template(C0max-1))-log(cross_template(C0max+1))) / (log(cross_template(C0max-1)) + log(cross_template(C0max+1)) - 2*log(cross_template(C0max))) / 2;
        a = C0max + subpixelCorrection;
        [~,C0max] = max(cross_template2);
        subpixelCorrection = (log(cross_template2(C0max-1))-log(cross_template2(C0max+1))) / (log(cross_template2(C0max-1)) + log(cross_template2(C0max+1)) - 2*log(cross_template2(C0max))) / 2;
        b = C0max + subpixelCorrection;
        shifts(j) = b-a; 
    end

    % shift of every second line by the average value computed above 
    movieL = zeros(size(movie));
    for j = 1:size(movie,3)
        for o = 1:size(movie,1)/2
            movieL((o-1)*2+2,:,j) = circshift(movie((o-1)*2+2,:,j),[0 0]);
            movieL((o-1)*2+1,:,j) = circshift(movie((o-1)*2+1,:,j),[0 round(median(shifts))]);
        end
    end
    disp(strcat('Median change:',12,num2str(median(shifts)),12,'Std change:',12,num2str(std(shifts))));

end
           
