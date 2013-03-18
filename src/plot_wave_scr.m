%*******************************************************************************
% * FILENAME:    plot_wave_scr.m
% * DESCRIPTION: subscript used for animation of a 3D meshplot 
% * AUTHOR:      James Matthew Welch [JMW]
% * SCHOOL:      Arizona State University
% * CLASS:       CSE598: High Performance Computing
% * INSTRUCTOR:  Dr. Gil Speyer
% * SECTION:     20520
% * TERM:        Spring 2013
% *******************************************************************************/
for(i=1:length(Z))
    if(Z(i) >1000000)
        Z(i)=4; 
    end
end
Z=reshape(Z,domSize,domSize);
h=figure;
if(plottype == 'mesh')
    plothandle = meshz(X,Y,Z);
else
    plothandle = surf(X, Y, Z); % ends up with lots of black
end

azimuth = 0;%-7;
elevation = 90;%43;
view(azimuth,elevation);