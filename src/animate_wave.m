close all
clear all
%path = '~/saguaroASU/cse598-HPC/p2/src/';
path = '';
gifname = 'wave.gif';
plottype = 'mesh';
prefix = 'output';
dirlist = dir([prefix,'*.txt']);
numframes = length(dirlist) - 2;
load([prefix,'1.txt']);
domSize = sqrt(length(output1));
images = zeros(domSize,domSize,numframes); 
imgseq = 1:2:numframes;

for(j=imgseq)
    fname = ['output',num2str(j)];
    fullfile = [path,fname,'.txt'];
    load(fullfile)
    eval(['output = ',fname,';'])
    eval(['clear ',fname])
    X=1:domSize;
    Y=1:domSize;
    Z_orig = output(:,end);
    Z=Z_orig;
    %avgZ = mean(Z);
    for(i=1:length(Z))
        if(Z(i) >1000000)
            Z(i)=4; 
        end
    end
    Z=reshape(Z,domSize,domSize);
    Z = Z - min(Z(:));
    images(:,:,j) = Z;
%end

%for(j = imgseq)
    amin = min(images(:,:,j));
    amax = max(images(:,:,j));
    im = mat2gray(images(:,:,j), [amin, amax]);
    [imind,cm] = gray2ind(im,256);
    if(j==1)
        imwrite(imind, cm, gifname, 'gif', 'Loopcount',inf);
    else
        imwrite(imind, cm, gifname, 'gif', 'WriteMode', 'append');
    end
end
disp([gifname,' written @ ',datestr(now)]);
