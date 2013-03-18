%*******************************************************************************
% * FILENAME:    show_high_vals.m
% * DESCRIPTION: displays matrix maxima
% * AUTHOR:      James Matthew Welch [JMW]
% * SCHOOL:      Arizona State University
% * CLASS:       CSE598: High Performance Computing
% * INSTRUCTOR:  Dr. Gil Speyer
% * SECTION:     20520
% * TERM:        Spring 2013
% *******************************************************************************/
indices=[];
init=1;
for(i=init:domSize)
    dim='col';
    if(dim=='row')
        mz = max(Z_orig(i,:));
    else
        mz = max(Z_orig(:,i));
    end
    if(mz > 10)
        disp([dim,'[',num2str(i-init),']']);
        indices(end+1)=i;
    end
end
            
