function parallelNewAlgebraic(minDt, maxDt, tFinal, regFac, N, bending)
minDt=str2num(minDt); 
maxDt=str2num(maxDt); 
tFinal=str2num(tFinal); 
regFac=str2num(regFac);
N = str2num(N); 
bending = str2num(bending);

poolobj = parpool;
fprintf('Number of workers: %g\n', poolobj.NumWorkers);
searching = 1; 
while(searching) 
    curDt = (minDt + maxDt)/2; 
    fprintf('trying candidate time step dt = %0.16e \n', curDt)
    tic
    isstable=parElasticNonlinearNewAlgebraic(curDt, tFinal, regFac,N,bending); 
    toc
    if isstable 
        minDt = curDt;
    else
        maxDt = curDt; 
    end 
    % disp(['next candidate critical time step = ' num2str(curDt)])
end
end
