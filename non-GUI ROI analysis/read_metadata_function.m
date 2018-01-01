% Peter Rupprecht 02-2015 using Scanimage B data
% reads metadata and returns image info file and framerate, zstep,
% framerate, zoom, motorpositions, scalingfactors
% optional input is the -full- filename/path
function [A,result,framerate,zstep,zoom,motorpositions,scalingfactors] = read_metadata_function(filename)

if strcmp(filename,'')
    %% read in path/file/folder of image
    PathName = '';
    [FileName,PathName,FilterIndex] = uigetfile(strcat(PathName,'*.tif')); file_name = strcat(PathName,FileName);
else
    file_name = filename;
end
A = imfinfo(file_name);


%% define metadata of interest 
snippet{1} = 'scanimage.SI4.scanFrameRate';
snippet{2} = 'scanimage.SI4.triggerOutDelay';
snippet{3} = 'scanimage.SI4.triggerOut';
snippet{4} = 'scanimage.SI4.stackZStepSize';
snippet{5} = 'scanimage.SI4.stackNumSlices';
snippet{6} = 'scanimage.SI4.scanPixelsPerLine';
snippet{7} = 'scanimage.SI4.scanZoomFactor';
snippet{8} = 'scanimage.SI4.scanLinesPerFrame';
snippet{9} = 'scanimage.SI4.savedBitdepthX';
snippet{10} = 'scanimage.SI4.triggerTime';
snippet{11} = 'scanimage.SI4.framerate_user';
snippet{12} = 'scanimage.SI4.framerate_user_check';
snippet{13} = 'scanimage.SI4.beamPowers';
snippet{14} = 'scanimage.SI4.beamLengthConstants';
snippet{15} = 'scanimage.SI4.autoscaleSavedImages';
snippet{16} = 'scanimage.SI4.acqNumFrames';
snippet{17} = 'scanimage.SI4.acqNumAveragedFrames';

snippet{18} = 'scanimage.SI4.motorPosition';
snippet{19} = 'scalingFactorAndOffset';
snippet{20} = 'framerate_precise';

%% read out metadata
clear result
for jj = 1:17
    k = strfind( A(1).ImageDescription(1:end),snippet{jj});
    check = 0;
    stringlength = 10;
    while check == 0
    result{jj} = str2num( A(1).ImageDescription(k + length(snippet{jj})+3:k + length(snippet{jj}) + 3 + stringlength) );
    if isempty(result{jj}); stringlength = stringlength - 1; else check = 1; end
    end
end

% exceptions

k = strfind( A(1).ImageDescription(1:end),snippet{18});
result{18} = A(1).ImageDescription(k + length(snippet{18})+3:k + length(snippet{18}) + 3 + 24);
k = strfind( A(1).ImageDescription(1:end),snippet{19});
try
    result{19} = A(1).ImageDescription(k + length(snippet{19})+3:k + length(snippet{19}) + 3 + 19);
catch
    result{19} = A(1).ImageDescription(k + length(snippet{19})+3:end);
end
k = strfind( A(1).ImageDescription(1:end),snippet{20});
try
    result{20} = A(1).ImageDescription(k + length(snippet{20})+3:k + length(snippet{20}) + 3 + 5);
catch
    result{20} = [];
end
%% make metadata it easier to process

if result{12}; framerate = min(result{11},result{1}); else framerate = result{1}; end
if ~isempty(result{20}); framerate = str2double(result{20}); end
zstep = result{4};
zoom = result{7};
motorpositions = result{18};
scalingfactors = result{19};

end