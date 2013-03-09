close all
clear all
load output.txt
domSize = sqrt(length(output));
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
figure;
meshz(X,Y,Z);
azimuth = -7;
elevation = 43;
view(azimuth,elevation);
%title(['Rot=',num2str(rotation),', Elev=',num2str(elevation)]);
interval = 50;
%axis([0 domSize+4 0 domSize+4 0 max(Z(:))]);
disp(['min = ',num2str(min(output))])
disp(['max = ',num2str(max(output))])
disp(['median = ',num2str(median(output))])
disp(['std = ',num2str(std(output))])
disp(['mean = ',num2str(mean(output))])

Z_orig=reshape(Z,domSize, domSize);

run('show_high_vals')

