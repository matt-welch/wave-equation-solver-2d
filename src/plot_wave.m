close all
clear all
load output.txt
X=1:480;
Y=1:480;
Z=reshape(output(:,3),480,480);
figure;
meshz(X,Y,Z);
axis([241-100 241+100 0 200 0 max(Z(:,:))]);
view(90,0);

