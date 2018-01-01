% calculate scale in x and y
function [pixelsize_x, pixelsize_y] = pixelsize_xy(zoom,width,width_undistored)

% Look-up-table, measured on the 14-01-2014 by moving the xy-stage using a fluorescent sample
% zsLUT = [1 870 934 ; 2 432 462; 3 283 302; 4 214 226; 8 111 120];
zsLUT = [1 913 944 ; 2 450 478; 3 307 319; 4 232 243; 8 119 122; 17 59.4 58.2];

% figure, plot(zsLUT(:,1),zsLUT(:,2),'k'); hold on; plot(zsLUT(:,1),zsLUT(:,3),'r');

% fit a function to data points for inter- and extrapolation
[xData, yData] = prepareCurveData( zsLUT(:,1), zsLUT(:,2) );
[xData2, yData2] = prepareCurveData( zsLUT(:,1), zsLUT(:,3) );
ft = fittype( 'a/x', 'independent', 'x', 'dependent', 'y');
opts = fitoptions( ft );
[fitresult, gof] = fit( xData, yData, ft, opts );
[fitresult2, gof] = fit( xData2, yData2, ft, opts);

pixelsize_y = fitresult2(zoom)/width;
pixelsize_x = fitresult(zoom)/width_undistored;
