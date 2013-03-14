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
            
