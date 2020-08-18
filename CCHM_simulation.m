%% The version will perform do-calculus and use causal effect score when the BIC score is equivalent


construct_verbose= false;
tol             = 10e-3;
maxCondSet      = 4;
pdSep           = 4;
heuristic       = 3;
nnSamples       = [10000];
nnLatent        = [0.1];
maxinDegree     = [3];
nnVarsL         = [10];
startingIter    = 1;
nIters          = 10;
threshold       = 0.01;
conservative = true;



for inSamples = 1:length(nnSamples)
    nSamples=nnSamples(inSamples);
    fprintf('--------------------------------nSamples: %d------------------------------------------------------------------\n', nSamples);


    
    for inLatent = 1:length(nnLatent)
        fprintf('--------------------------------nLatent: %.4f----------------------------------------\n', nnLatent(inLatent));
        pLatent = 100*nnLatent(inLatent);
     
        

        
        % Control random number generator
        rng(0,'combRecursive');
        
        for inmaxinDegree= 1:length(maxParents)
            nmaxinDegree= maxParents(inMaxParents);
            fprintf('--------------------------------nMaxinDegree: %d-----------------------------------\n', nMaxParents);
       
            for inVarsL = 1:length(nnVarsL)
                nVarsL = nnVarsL(inVarsL);
                fprintf('--------------------------------nVars: %d--------------------\n', nVarsL);
                nLatent = ceil(nnLatent(inLatent)*nVarsL);
               % nLatent=1;
                fprintf('--------------------------------totalLatent: %d--------------------\n', nLatent);
             
                

                for iter = startingIter:nIters
                    fprintf('Iter %d:\n', iter);
                    
                    % Control rng behavior
                    stream = RandStream.getGlobalStream();
                    % Unique Substream for each combination of iter-nVarsL-nMaxParents
                    stream.Substream = (iter + 60*nVarsL + 2100*nMaxParents);
                 
                    % Generate new data
                    dag = randomdag(nVarsL, nMaxParents);
                
                    % Choose latent variables.
               
                    isLatent = false(1, nVarsL);
                    isLatent(randsample(1:nVarsL, nLatent)) = true;
                     Lat = find(isLatent)';

                 

 
                    % Simulate data
                    bn = dag2randBN(dag, 'gaussian');
                    ds = simulatedata(bn, nSamples, 'gaussian', 'isLatent', isLatent);
                    
                    % Create the true MAG/PAG
                    magT = convertDagToMag(dag, Lat, []);
                    pagT = mag2pag(magT);
                    nVars = sum(~isLatent);
                    
                  
                    % Compute the covariance/correlation matrix
                  
                      covMat = cov(ds.data(:,~isLatent));
                      rawdata=ds.data(:, ~isLatent);

                    % Run CCHM
                    tic;
                    ds.data= ds.data(:,~isLatent);            
                    [cchmMag,time, chcmIters,pairs] = CCHM(ds, nSamples, tol, false,construct_verbose,heuristic, threshold, maxCondSet, pdSep);

                    metric_runtime(iter)=toc;
  
                 %    convert to MAG
                    cchmMag = ag2mag(cchmMag);
                       
%                    convert to PAG
                    cchmPag = mag2pag(cchmMag);
                    
    
%    find numbers edge

                     nEdges(iter) = 0.5*nnz((cchmPag)~=0); 
                     nEdgesGT(iter) = 0.5*nnz((pagT)~=0);
                    
    
  
%% Run measurment metrics

                metric_SHDcchm(iter)=structuralHammingDistancePAG(cchmPag, pagT);
                [metric_precision(iter), metric_recall(iter),metric_bsf(iter)] = precisionRecall(cchmPag, pagT);

              


               
                    
                end
            end
        end
    end
end

