close all
clear all
%path = '~/saguaroASU/cse598-HPC/p2/src/';
path = '';
plottype = 'mesh';
for(i=1:100)
    fname = ['output',num2str(i)];
    fullfile = [path,fname,'.txt'];
    load(fullfile)
    eval(['output = ',fname,';'])
    domSize = sqrt(length(output));
    X=1:domSize;
    Y=1:domSize;
    Z_orig = output(:,end);
    Z=Z_orig;
    %avgZ = mean(Z);    
    run('plot_wave_scr')
    title(fullfile);
    pause(0.5);
    close(h)
end

disp(['min = ',num2str(min(output))])
disp(['max = ',num2str(max(output))])
disp(['median = ',num2str(median(output))])
disp(['std = ',num2str(std(output))])
disp(['mean = ',num2str(mean(output))])

Z_orig=reshape(Z,domSize, domSize);

run('show_high_vals')

