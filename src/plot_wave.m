%*******************************************************************************
% * FILENAME:    plot_wave.m
% * DESCRIPTION: Plots a single mesh of a 2-D wave domain
% * AUTHOR:      James Matthew Welch [JMW]
% * SCHOOL:      Arizona State University
% * CLASS:       CSE598: High Performance Computing
% * INSTRUCTOR:  Dr. Gil Speyer
% * SECTION:     20520
% * TERM:        Spring 2013
% *******************************************************************************/

close all
clear all
%path = '~/saguaroASU/cse598-HPC/p2/src/';
path = '';
defPath = '/home/mwelch/Dropbox/classes/cse598-HPC/projects/p2/src/';
defFileBase = 'output';
defFileExt = 'txt';
defFile = [defFileBase '.' defFileExt];
defFullFile = '';
uiTitle = 'Select a Wave Solver output file';
uiFilterSpec = {'*.txt','TEXT files'};
defFullFile = [defPath, defFile];
[FILENAME, PATHNAME] = uigetfile(uiFilterSpec, uiTitle, defFullFile);
if(FILENAME ~= 0)
    fullfile = [PATHNAME FILENAME]; 
    [pathstr,fname,ext] = fileparts(FILENAME);
else
    disp('No File Selected, aborting. ');
%     disp(defFullFile);
%     fullfile = defFullFile;
%     fname = 'output';
    return
end
load(fullfile)
eval(['output = ',fname,';'])
plottype = 'mesh';
domSize = sqrt(length(output));
X=1:domSize;
Y=1:domSize;
Z_orig = output(:,end);
Z=Z_orig;
%avgZ = mean(Z);
for i=1:length(Z)
    if(Z(i) >1000000)
        Z(i)=4; 
    end
end
Z=reshape(Z,domSize,domSize);
h=figure;
if(strcmp(plottype,'mesh'))
    meshz(X,Y,Z);
else
    surf(X, Y, Z); % ends up with lots of black
end

azimuth = 0;%-7;
elevation = 90;%43;
view(azimuth,elevation);
title(fullfile);

%title(['Rot=',num2str(rotation),', Elev=',num2str(elevation)]);
%interval = 50;
%axis([0 domSize+4 0 domSize+4 0 max(Z(:))]);
disp(['min = ',num2str(min(output))])
disp(['max = ',num2str(max(output))])
disp(['median = ',num2str(median(output))])
disp(['std = ',num2str(std(output))])
disp(['mean = ',num2str(mean(output))])

Z_orig=reshape(Z,domSize, domSize);

run('show_high_vals')

