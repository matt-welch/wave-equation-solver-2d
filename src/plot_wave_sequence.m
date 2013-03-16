close all
clear all
%path = '~/saguaroASU/cse598-HPC/p2/src/';
path = '';
plottype = 'mesh';
prefix = 'output';
dirlist = dir([prefix,'*.txt']);
numframes = length(dirlist) - 2;
fname = ['output1'];
fullfile = [path,fname,'.txt'];
load(fullfile)
eval(['output = ',fname,';'])
domSize = sqrt(length(output));
X=1:domSize;
Y=1:domSize;
Z_orig = output(:,end);
Z=Z_orig;
Z=reshape(Z,domSize,domSize);
h=figure;
if(plottype == 'mesh')
    plothandle = meshz(X,Y,Z);
else
    plothandle = surf(X, Y, Z); % ends up with lots of black
end

azimuth = -7;%0;%-7;
elevation = 43;%90;%43;
view(azimuth,elevation);
title(fullfile);

for(i=2:5:numframes)
    fname = ['output',num2str(i)];
    fullfile = [path,fname,'.txt'];
    load(fullfile)
    eval(['output = ',fname,';'])
    eval(['clear ',fname,';'])
    Z = reshape(output(:,end), domSize, domSize);
     temp = get(plothandle, 'ZData');
    temp(1:domSize,1:domSize) = Z;
    set(plothandle,'ZData',temp);
    title(fullfile);
    drawnow();
    pause(0.25);
end

disp(['min = ',num2str(min(output))])
disp(['max = ',num2str(max(output))])
disp(['median = ',num2str(median(output))])
disp(['std = ',num2str(std(output))])
disp(['mean = ',num2str(mean(output))])

Z_orig=reshape(Z,domSize, domSize);

run('show_high_vals')

