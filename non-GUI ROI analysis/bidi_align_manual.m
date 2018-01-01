% bidirectionally align 3D movie, manual shift
function movieL = bidi_align_manual(movie,shift)


movieL = zeros(size(movie));

for j = 1:size(movie,3)
    for o = 1:size(movie,1)/2
        movieL((o-1)*2+2,:,j) = circshift(movie((o-1)*2+2,:,j),[0 0]);
        movieL((o-1)*2+1,:,j) = circshift(movie((o-1)*2+1,:,j),[0 -shift]);
    end
end

end