close all
clear all
load output.txt
X=1:480;
Y=1:480;
Z = output(:,end);
%avgZ = mean(Z);
%for(i=1:length(Z))
%	if(Z(i) >avgZ)
%		Z(i)=avgZ; 
%	endif
%endfor
Z=reshape(Z,480,480);
figure;
meshz(X,Y,Z);
view(100,15);
interval = 50;
%axis([241-interval 241+interval 0 2*interval 0 max(max(Z(:,:)))]);


