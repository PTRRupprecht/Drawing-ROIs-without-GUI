function [movie,offsety_resolved,offsetx_resolved] = alginWithinTrial(movie,reference)

A = reference;

% save('reference_image.mat','A');

counter = 1;
for i = 1:1:(size(movie,3))

    B = mean(movie(:,:,i),3); % circshift(movie(:,:,1),[-3 +15]);

    result_conv =fftshift(real(ifft2(conj(fft2(A)).*fft2(B))));
    [y,x] = find(result_conv==max(result_conv(:))); %Find the 255 peak
    result_conv =fftshift(real(ifft2(conj(fft2(A)).*fft2(A))));
    [y0,x0] = find(result_conv==max(result_conv(:))); %Find the 255 peak
    offsety(counter) = y-y0;
    offsetx(counter) = x-x0;
    counter = counter + 1;
end

windowsize = 30;

tempx = round(medfilt1(offsetx,windowsize));
tempy = round(medfilt1(offsety,windowsize));

% offsety(30:end-30) = tempy(30:end-30);
% offsetx(30:end-30) = tempx(30:end-30);

offsety_resolved = offsety;% interp1(1:10:(size(movie,3)-9),offsety,1:size(movie,3),'nearest','extrap');
offsetx_resolved = offsetx;% interp1(1:10:(size(movie,3)-9),offsetx,1:size(movie,3),'nearest','extrap');

% figure(199), hold on, plot(offsety_resolved), plot(offsetx_resolved,'r'), 
if sum(offsety_resolved) ~=0 || sum(offsetx_resolved) ~=0 
    keyboard
    if 0
        offsetx_resolved(1:20) = offsetx_resolved(21);
        offsety_resolved(1:20) = offsety_resolved(21);
        offsetx_resolved(end-20:end) = offsetx_resolved(end-21);
        offsety_resolved(end-20:end) = offsety_resolved(end-21);
        offsety_resolved = round(smooth(offsety_resolved,80));
        offsetx_resolved = round(smooth(offsetx_resolved,80));
    end
end


% offsety_resolved(offsety_resolved ~= -2) = -2;
% offsetx_resolved(offsetx_resolved ~= 0) = 0;

for i = 1:size(movie,3)
    if mod(i,100) == 0; disp(i/size(movie,3)); end
    movie(:,:,i) = circshift(movie(:,:,i),[-offsety_resolved(i),-offsetx_resolved(i)]);
end

end