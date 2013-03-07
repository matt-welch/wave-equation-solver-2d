close all
clear all
load output.txt
domSize = sqrt(length(output));
X=1:domSize;
Y=1:domSize;
Z = output(:,end);
%avgZ = mean(Z);
%for(i=1:length(Z))
%	if(Z(i) >avgZ)
%		Z(i)=avgZ; 
%	endif
%endfor
Z=reshape(Z,domSize,domSize);
figure;
meshz(X,Y,Z);
rotation = 330;
elevation = 15;
view(rotation,elevation);
title(['Rot=',num2str(rotation),', Elev=',num2str(elevation)]);
interval = 50;
%axis([241-interval 241+interval 0 2*interval 0 max(max(Z(:,:)))]);
disp(['Min = ',num2str(min(output))])
disp(['Max = ',num2str(max(output))])


