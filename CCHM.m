function [mag, iter,gs, pairs] = CCHM(ds, nSamples, tol, verbose,construct_verbose,heuristic, threshold, maxCondSet, pdSep)

covMat=  cov(ds.data);
nVarsL = size(ds.data,2);



%% Phase 1 find the skeleton and blacklist, whitelist



%               if strcmp(ds.type, 'continuouse')

[FCI_skeleton,sepSets,pvalue,maxSepSet,~] = fciskeleton(ds, heuristic, threshold, maxCondSet, pdSep, false);





%              elseif strcmp(ds.type ,'discrete')
%              [FCI_skeleton,sepSets,pvalue,maxSepSet,~]= fciskeleton_all_type(ds,'g2test_2', heuristic, threshold, maxCondSet,  pdSep, false);
% 
% 
%              end


[pag, dnc, cols, unfaithful, triangles,pvalue2] = R0_conservative(FCI_skeleton, sepSets, @fisher, ds, threshold, maxCondSet, false);

ass_score=pvalue;

for i = 1:nVarsL
    skeleton_final{1,i} = find(pag(:,i)>=1);
end



%Change o- to -

nvars = size(pag,1);
for i = 1:nvars
    for j=1:nvars
        if(pag(i,j)==1)
            pag(i,j) = 3;
            
        end
    end
end

nVars=size(covMat,1);
% begin from the empty mag
isAncestor = false(nVars);
stepAncestor = zeros(nVars);
isParent= false(nVars);
isBidir = false(nVars);
whichPair=nan;
pairs_withpvalue = zeros(nVars,3);

% initially each node is a component changed to --> FCI-M3HC jump to R0
mag =zeros(nVars);
omega_mat =zeros(nVars);
beta_mat =zeros(nVars);
omega_mat1 =zeros(nVars);
beta_mat1 =zeros(nVars);
omega_mat2 =zeros(nVars);
beta_mat2 =zeros(nVars);
omega_mat3 =zeros(nVars);
beta_mat3 =zeros(nVars);
[nComps, sizes, comps, inComponent]= concomp(mag);
nsf = -nSamples/2;
scores = zeros(1,nComps);
nEdges=0;

for iComp =1:nComps
    [compMag, district] = componentMag(comps{iComp}, nVars, mag, isParent);
    [scores(iComp),beta,omega]  = score_contrib(compMag, comps{iComp}, district, sizes(iComp), covMat,nSamples, tol);
    
end



tmpSll = nsf*sum(scores);
CurScore = -2*tmpSll+bicPenalty(nSamples,nVars, nEdges);
if construct_verbose  
    fprintf('initiate CurScore= %.4f\n', CurScore);
end




IterMag=cell(100,1);
iter=1;
gs = struct();
gs(iter).score=CurScore;


pairs = [];
for j = 1:length(skeleton_final)
    for i = 1:length(skeleton_final{j})
        if j<skeleton_final{j}(i)
            pairs = [pairs; [j, skeleton_final{j}(i)]];
            
        else
        end
    end
end




nPairs = size(pairs,1);
%% from , to, visit? , mem_Action, minScore, prohibit action -> for dnc, prohibit action <- for dnc,prohibit action <-> for dnc, p-value, collider p-value
whichPair = pairs;

visit = false(nPairs,1);
whichPair = [whichPair visit];
which_operation = zeros(nPairs,1);
whichPair = [whichPair which_operation];
candidate_score = Inf(nPairs,1);
whichPair = [whichPair candidate_score];
prohi_action = zeros(nPairs,3);
whichPair = [whichPair prohi_action];
%bidirect and directed PAG

% use p-value for collider to sort



for iPair =1:nPairs
    pairs(iPair, 3)= pvalue2(pairs(iPair, 1),pairs(iPair, 2));
    
    
    
    if size(cols) ~=0
        
        
        % find a-b-*
        [row1,~] = find(cols(:,1)==pairs(iPair, 1) & cols(:,2)==pairs(iPair, 2));
        
        if size(row1,1)>0
            
            
            for i_row=1:size(row1,1)
                
                a=cols(row1(i_row),1);
                b=cols(row1(i_row),2);
                c=cols(row1(i_row),3);
                
                
                pairs(iPair, 4)= pvalue(a,b)+pvalue(b,c);
                
            end
        end
        % find b-a-*
        [row1,~] = find(cols(:,1)==pairs(iPair, 2) & cols(:,2)==pairs(iPair, 1));
        
        if size(row1,1)>0
            
            for i_row=1:size(row1,1)
                
                a=cols(row1(i_row),1);
                b=cols(row1(i_row),2);
                c=cols(row1(i_row),3);
                
                pairs(iPair, 4)= pvalue(a,b)+pvalue(b,c);
                
            end
            
        end
        
        % find *-a-b
        [row1,~] = find(cols(:,2)==pairs(iPair, 1) & cols(:,3)==pairs(iPair, 2));
        
        if size(row1,1)>0
            
            for i_row=1:size(row1,1)
                
                a=cols(row1(i_row),1);
                b=cols(row1(i_row),2);
                c=cols(row1(i_row),3);
                
                pairs(iPair, 4)= pvalue(a,b)+pvalue(b,c);
                
            end
        end
        % find *-b-a
        [row1,~] = find(cols(:,2)==pairs(iPair, 2) & cols(:,3)==pairs(iPair, 1));
        
        if size(row1,1)>0
            
            for i_row=1:size(row1,1)
                
                a=cols(row1(i_row),1);
                b=cols(row1(i_row),2);
                c=cols(row1(i_row),3);
                
                pairs(iPair, 4)= pvalue(a,b)+pvalue(b,c);
                
            end
            
        end
    end
    
    
end
if size(cols) ==0
    
    whichPair = [whichPair  zeros(nPairs)];
    whichPair = [whichPair zeros(nPairs)];
else
    whichPair = [whichPair pairs(:,3)];
    whichPair = [whichPair pairs(:,4)];
end



% heuristic collider P-value min to max sortrows(whichPair,9) 
% sortrows(whichPair,10) pvalue for collider

whichPair = sortrows(whichPair,9);




for iPair =1:nPairs
    
    from = pairs(iPair, 1); to = pairs(iPair, 2);
    
    
    % add bidirected
    if pag(from, to)==2 && pag(to, from)==2
        
        if  iter==1
            
            % fprintf('iter %d\n', iter)
            %  fprintf('%d<->%d already in the graph\n', from, to)
            nEdges=nEdges+1;
            [CurScore, scoreContribs, newNComps, newSizes, newComps, newInComponent, IterMag{iter},beta,omega,district] = addBidirectedEdge(from, to, nEdges, nComps, sizes, comps, inComponent, mag, covMat, isParent, scores, nSamples, nVars,  tol);
            gs(iter).score=CurScore;
            
            iter=iter+1;
            
            
        elseif iter>1
            
            mag_check=IterMag{iter-1};
            mag_check(from,to)=2;
            mag_check(to,from)=2;
            if isAcyclic(mag_check)
                
                %fprintf('iter %d\n', iter)
                %fprintf('%d<->%d already in the graph\n', from, to)
                nEdges=nEdges+1;
                [CurScore, scoreContribs, newNComps, newSizes, newComps, newInComponent, IterMag{iter},beta,omega,district] = addBidirectedEdge(from, to, nEdges,newNComps, newSizes, newComps, newInComponent,IterMag{iter-1}, covMat, isParent, scoreContribs, nSamples, nVars,  tol);
                gs(iter).score=CurScore;
                
                
                iter=iter+1;
            else
                fprintf('There is an error');
            end
        end
        
    end
    
    
    % add directed ->
    if pag(from, to)==2 && pag(to, from)==3
        
        if  iter==1
            
            %fprintf('iter %d\n', iter)
            % fprintf('%d->%d already in the graph\n', from, to)
            isParent(from,to)=true;
            nEdges=nEdges+1;
            [CurScore, scoreContribs, newNComps, newSizes, newComps, newInComponent, IterMag{iter},beta,omega,district] = addDirectedEdge(from, to, nEdges, nComps, sizes, comps, inComponent, mag, covMat, isParent, scores, nSamples, nVars,  tol);
            
            gs(iter).score=CurScore;
            iter=iter+1;
            
        elseif iter>1
            mag_check=IterMag{iter-1};
            mag_check(from,to)=2;
            mag_check(to,from)=3;
            if isAcyclic(mag_check)
                % fprintf('iter %d\n', iter)
                %fprintf('%d->%d already in the graph\n', from, to)
                isParent(from,to)=true;
                nEdges=nEdges+1;
                [CurScore, scoreContribs, newNComps, newSizes, newComps, newInComponent, IterMag{iter},beta,omega,district] = addDirectedEdge(from, to, nEdges,newNComps, newSizes, newComps, newInComponent,IterMag{iter-1}, covMat, isParent, scoreContribs, nSamples, nVars,  tol);
                
                gs(iter).score=CurScore;
                iter=iter+1;
            else
                % change colliders to undirected graph as it creates cyclic
                
                fprintf('There is an error from %d to %d \n',from,to);
                [row1,~] = find(cols(:,1)==from & cols(:,2)==to);
                if size(row1,1)==1
                    
                    a=cols(row1,1);
                    b=cols(row1,2);
                    c=cols(row1,3);
                    
                    pag(a,b)=3;
                    pag(b,a)=3;
                    pag(b,c)=3;
                    pag(c,b)=3;
                    %delete collider *-b-c
                    if IterMag{iter-1}(c,b)==2 &  IterMag{iter-1}(b,c)==3
                        isParent(c,b)=false;
                        nEdges=nEdges-1;
                        [CurScore, scoreContribs, newNComps, newSizes, newComps, newInComponent, IterMag{iter},beta,omega,district] = removeDirectedEdge(c, b, nEdges,newNComps, newSizes, newComps, newInComponent,IterMag{iter-1}, covMat, isParent, scoreContribs, nSamples, nVars,  tol);
                        
                        gs(iter).score=CurScore;
                        iter=iter+1;
                        
                        
                    end
                    
                    
                end
                
            end
        end
        
    end
    
    % add directed <-
    if pag(from, to)==3 && pag(to, from)==2
        
        if  iter==1
            
            %fprintf('iter %d\n', iter)
            % fprintf('%d->%d already in the graph\n', to, from)
            isParent(to,from)=true;
            nEdges=nEdges+1;
            [CurScore, scoreContribs, newNComps, newSizes, newComps, newInComponent, IterMag{iter},beta,omega,district] = addDirectedEdge(to, from, nEdges, nComps, sizes, comps, inComponent, mag, covMat, isParent, scores, nSamples, nVars,  tol);
            
            gs(iter).score=CurScore;
            iter=iter+1;
            
        elseif iter>1
            
            mag_check=IterMag{iter-1};
            mag_check(from,to)=3;
            mag_check(to,from)=2;
            
            if isAcyclic(mag_check)
                
                
                % fprintf('iter %d\n', iter)
                % fprintf('%d->%d already in the graph\n', to, from)
                isParent(to,from)=true;
                nEdges=nEdges+1;
                [CurScore, scoreContribs, newNComps, newSizes, newComps, newInComponent, IterMag{iter},beta,omega,district] = addDirectedEdge(to, from, nEdges,newNComps, newSizes, newComps, newInComponent,IterMag{iter-1}, covMat, isParent, scoreContribs, nSamples, nVars,  tol);
                
                gs(iter).score=CurScore;
                iter=iter+1;
            else
                % change colliders to undirected graph as it creates cyclic
                
                fprintf('There is an error from %d to %d ',to,from);
                [row1,~] = find(cols(:,1)==to & cols(:,2)==from);
                if size(row1,1)==1
                    
                    a=cols(row1,1);
                    b=cols(row1,2);
                    c=cols(row1,3);
                    
                    pag(a,b)=3;
                    pag(b,a)=3;
                    pag(b,c)=3;
                    pag(c,b)=3;
                    if IterMag{iter-1}(c,b)==2 &  IterMag{iter-1}(b,c)==3
                        isParent(c,b)=false;
                        nEdges=nEdges-1;
                        [CurScore, scoreContribs, newNComps, newSizes, newComps, newInComponent, IterMag{iter},beta,omega,district] = removeDirectedEdge(c, b, nEdges,newNComps, newSizes, newComps, newInComponent,IterMag{iter-1}, covMat, isParent, scoreContribs, nSamples, nVars,  tol);
                        
                        gs(iter).score=CurScore;
                        
                        iter=iter+1;
                    end
                    
                end
                
            end
        end
        
    end
    
    
    
    
    
end
% fprintf('direct phase %d  \n',iter);




%undirected PAG

for iPair =1:nPairs
    
    from = whichPair(iPair, 1); to = whichPair(iPair, 2);
    
    if pag(from, to)==3 && pag(to, from)==3
        whichPair(iPair,3)= 1;
        
        if size(dnc) ~=0
            
            [row1,~] = find((dnc(:,1)==from & dnc(:,2)==to));
            [row2,~] = find((dnc(:,1)==to & dnc(:,2)==from));
            [row3,~] = find((dnc(:,2)==from & dnc(:,3)==to));
            [row4,~] = find((dnc(:,2)==to & dnc(:,3)==from));
            
            if  iter==1
                
                
                
                
                %% Node Operation A->B checked Acyclic
                
                nEdges=nEdges+1;
                isParent(from,to)=true;
                [score1,~,~,~,~,~,Mag_A,beta1,omega1,district] = addDirectedEdge(from, to, nEdges, nComps, sizes, comps, inComponent, mag, covMat, isParent, scores, nSamples, nVars,  tol);
                isParent(from,to)=false;
                
                
                omega_mat1(district,district)= omega1;
                beta_mat1(district,district)= beta1;
                
                causal_effect1=abs(beta_mat1(to,from));
                
                if ~isAcyclic(Mag_A)
                    score1 =Inf;
                end
                
                
                
                
                
                
                %% Node Operation A <- B checked Acyclic
                isParent(to,from)=true;
                [score2,~,~,~,~,~,Mag_B,beta2,omega2,district] = addDirectedEdge(to, from, nEdges, nComps, sizes, comps, inComponent, mag, covMat, isParent, scores, nSamples, nVars,  tol);
                isParent(to,from)=false;
                
                
                
                omega_mat2(district,district)= omega2;
                beta_mat2(district,district)= beta2;
                
                causal_effect2=abs(beta_mat2(from,to));
                if ~isAcyclic(Mag_B)
                    score2 =Inf;
                end
                
                
                [score3,~,~,~,~,~,~,~,~,~] = addBidirectedEdge(from, to, nEdges, nComps, sizes, comps, inComponent, mag, covMat, isParent, scores, nSamples, nVars,  tol);
                nEdges=nEdges-1;
                
                
                
                
                
            elseif iter>1  & IterMag{iter-1}(dnc(row1,3),dnc(row1,2))==2
                whichPair(iPair,6)=1;
                whichPair(iPair,8)=1;
                score1=Inf;
                nEdges=nEdges+1;
                isParent(to,from)=true;
                [score2,~,~,~,~,~,Mag_B,beta1,omega1,district] = addDirectedEdge(to, from, nEdges, newNComps, newSizes, newComps, newInComponent, IterMag{iter-1}, covMat, isParent, scoreContribs, nSamples, nVars,  tol);
                nEdges=nEdges-1;
                isParent(to,from)=false;
                score3=Inf;
                
                
                beta_mat1(district,district)= beta1;
                omega_mat1(district,district)= omega1;
                causal_effect1=abs(beta_mat1(to,from));
                
                
                if ~isAcyclic(Mag_B)
                    score2 =Inf;
                end
                
                %fprintf('else 2 have dnc= %d<--%d ',dnc(row1,2),dnc(row1,3))
                %fprintf('avoid dnc= %d %d %d',dnc(row1,1),dnc(row1,2),dnc(row1,3))
                % %
                
            elseif iter>1 & IterMag{iter-1}(dnc(row2,3),dnc(row2,2))==2
                whichPair(iPair,7)=1;
                whichPair(iPair,8)=1;
                
                nEdges=nEdges+1;
                isParent(from,to)=true;
                [score1,~,~,~,~,~,Mag_A,beta1,omega1,district] = addDirectedEdge(from, to, nEdges, newNComps, newSizes, newComps, newInComponent, IterMag{iter-1}, covMat, isParent, scoreContribs, nSamples, nVars,  tol);
                isParent(from,to)=false;
                nEdges=nEdges-1;
                
                omega_mat1(district,district)= omega1;
                beta_mat1(district,district)= beta1;
                
                causal_effect1=abs(beta_mat1(from,to));
                
                if ~isAcyclic(Mag_A)
                    score1 =Inf;
                end
                score2=Inf;
                score3=Inf;
                
                %                 fprintf('else 3 have dnc= %d<--%d ',dnc(row2,2),dnc(row2,3))
                %                 fprintf('avoid dnc= %d %d %d',dnc(row2,1),dnc(row2,2),dnc(row2,3))
                %
                % %
                
            elseif iter>1 & IterMag{iter-1}(dnc(row3,1),dnc(row3,2))==2
                whichPair(iPair,7)=1;
                whichPair(iPair,8)=1;
                nEdges=nEdges+1;
                isParent(from,to)=true;
                [score1,~,~,~,~,~,Mag_A,beta1,omega1,district] = addDirectedEdge(from, to, nEdges, newNComps, newSizes, newComps, newInComponent, IterMag{iter-1}, covMat, isParent, scoreContribs, nSamples, nVars,  tol);
                nEdges=nEdges-1;
                isParent(from,to)=false;
                score2=Inf;
                score3=Inf;
                omega_mat1(district,district)= omega1;
                beta_mat1(district,district)= beta1;
                
                causal_effect1=abs(beta_mat1(to,from));
                if ~isAcyclic(Mag_A)
                    score1 =Inf;
                end
                %
                %                 fprintf('else 3 have dnc= %d<--%d ',dnc(row3,2),dnc(row3,3))
                %                 fprintf('avoid dnc= %d %d %d',dnc(row3,1),dnc(row3,2),dnc(row3,3))
                %
                
            elseif iter>1 & IterMag{iter-1}(dnc(row4,1),dnc(row4,2))==2
                whichPair(iPair,7)=1;
                whichPair(iPair,8)=1;
                
                score1=Inf;
                
                nEdges=nEdges+1;
                isParent(to,from)=true;
                [score2,~,~,~,~,~,Mag_B,beta2,omega2,district] = addDirectedEdge(to, from, nEdges, newNComps, newSizes, newComps, newInComponent, IterMag{iter-1}, covMat, isParent, scoreContribs, nSamples, nVars,  tol);
                isParent(to,from)=false;
                nEdges=nEdges-1;
                if ~isAcyclic(Mag_B)
                    score2 =Inf;
                end
                
                omega_mat2(district,district)= omega2;
                beta_mat2(district,district)= beta2;
                
                causal_effect2=abs(beta_mat2(from,to));
                
                score3=Inf;
                %
                %                 fprintf('else 3 have dnc= %d<--%d ',dnc(row4,2),dnc(row4,3))
                %                 fprintf('avoid dnc= %d %d %d',dnc(row4,1),dnc(row1,2),dnc(row4,3))
                %
                % %
            else
                nEdges=nEdges+1;
                isParent(from,to)=true;
                [score1,~,~,~,~,~,Mag_A,beta1,omega1,district] = addDirectedEdge(from, to, nEdges, newNComps, newSizes, newComps, newInComponent, IterMag{iter-1}, covMat, isParent, scoreContribs, nSamples, nVars,  tol);
                isParent(from,to)=false;
                omega_mat1(district,district)= omega1;
                beta_mat1(district,district)= beta1;
                
                causal_effect1=abs(beta_mat1(to,from));
                %                          causal_effect1
                
                if ~isAcyclic(Mag_A)
                    score1 =Inf;
                end
                
                
                isParent(to,from)=true;
                [score2,~,~,~,~,~,Mag_B,beta2,omega2,district] = addDirectedEdge(to, from, nEdges, newNComps, newSizes, newComps, newInComponent, IterMag{iter-1}, covMat, isParent, scoreContribs, nSamples, nVars,  tol);
                isParent(to,from)=false;
                
                omega_mat2(district,district)= omega2;
                beta_mat2(district,district)= beta2;
                
                causal_effect2=abs(beta_mat2(from,to));
                %                          causal_effect2
                if ~isAcyclic(Mag_B)
                    score2 =Inf;
                end
                
                
                
                
                [score3,~,~,~,~,~,~,~,~,~] = addBidirectedEdge(from, to, nEdges, newNComps, newSizes, newComps, newInComponent, IterMag{iter-1}, covMat, isParent, scoreContribs, nSamples, nVars,  tol);
                
                nEdges=nEdges-1;
            end
            
        elseif  iter==1
            
            nEdges=nEdges+1;
            isParent(from,to)=true;
            [score1,~,~,~,~,~,Mag_A,beta1,omega1,district] = addDirectedEdge(from, to, nEdges, nComps, sizes, comps, inComponent, mag, covMat, isParent, scores, nSamples, nVars,  tol);
            isParent(from,to)=false;
            
            %                 district
            
            omega_mat1(district,district)= omega1;
            beta_mat1(district,district)= beta1;
            
            causal_effect1=abs(beta_mat1(to,from));
            %                 causal_effect1
            %
            %                 beta_mat1
            %                 omega_mat1
            
            if ~isAcyclic(Mag_A)
                score1 =Inf;
            end
            
            
            isParent(to,from)=true;
            [score2,~,~,~,~,~,Mag_B,beta2,omega2,district] = addDirectedEdge(to, from, nEdges, nComps, sizes, comps, inComponent, mag, covMat, isParent, scores, nSamples, nVars,  tol);
            isParent(to,from)=false;
            
            %                 district
            
            omega_mat2(district,district)= omega2;
            beta_mat2(district,district)= beta2;
            
            causal_effect2=abs(beta_mat2(from,to));
            %                 causal_effect2
            %
            %                 beta_mat2
            %                 omega_mat2
            
            if ~isAcyclic(Mag_B)
                score2 =Inf;
            end
            
            
            [score3,~,~,~,~,~,~,beta3,omega3,district] = addBidirectedEdge(from, to, nEdges, nComps, sizes, comps, inComponent, mag, covMat, isParent, scores, nSamples, nVars,  tol);
            %                 district
            
            omega_mat3(district,district)= omega3;
            beta_mat3(district,district)= beta3;
            error_term=(omega_mat3(from,to));
            %
            %                 beta_mat3
            %                 omega_mat3
            
            
            nEdges=nEdges-1;
            
            
            
            if construct_verbose
                
                fprintf('from  %d to %d \n',from,to);
                fprintf('score1: %0.10f causal effect1=%0.10f \n',score1,causal_effect1);
                fprintf('score2: %0.10f causal effect2=%0.10f \n',score2,causal_effect2);
                fprintf('score3: %0.10f error_term3=%0.10f\n',score3,error_term);
            end
            
        else
            nEdges=nEdges+1;
            isParent(from,to)=true;
            [score1,~,~,~,~,~,Mag_A,omega1,beta1,district] = addDirectedEdge(from, to, nEdges, newNComps, newSizes, newComps, newInComponent, IterMag{iter-1}, covMat, isParent, scoreContribs, nSamples, nVars,  tol);
            isParent(from,to)=false;
            
            
            
            omega_mat1(district,district)= omega1;
            beta_mat1(district,district)= beta1;
            causal_effect1=abs(beta_mat1(to,from));
            %                 causal_effect1
            %                 beta_mat1
            %                 omega_mat1
            
            if ~isAcyclic(Mag_A)
                score1 =Inf;
            end
            
            
            isParent(to,from)=true;
            [score2,~,~,~,~,~,Mag_B,omega2,beta2,district] = addDirectedEdge(to, from, nEdges, newNComps, newSizes, newComps, newInComponent, IterMag{iter-1}, covMat, isParent, scoreContribs, nSamples, nVars,  tol);
            isParent(to,from)=false;
            
            
            omega_mat2(district,district)= omega2;
            beta_mat2(district,district)= beta2;
            causal_effect2=abs(beta_mat2(from,to));
            %                 causal_effect2
            %                 beta_mat1
            %                 omega_mat1
            if ~isAcyclic(Mag_B)
                score2 =Inf;
            end
            
            
            [score3,~,~,~,~,~,~,omega3,beta3,district] = addBidirectedEdge(from, to, nEdges, newNComps, newSizes, newComps, newInComponent, IterMag{iter-1}, covMat, isParent, scoreContribs, nSamples, nVars,  tol);
            omega_mat3(district,district)= omega3;
            beta_mat3(district,district)= beta3;
            
            omega_mat3(district,district)= omega3;
            beta_mat3(district,district)= beta3;
            
            nEdges=nEdges-1;
            fprintf('iter:%d from  %d to %d \n',iter,from,to);
            
            fprintf('score1: %0.10f causal effect1=%0.10f \n',score1,causal_effect1);
            fprintf('score2: %0.10f causal effect2=%0.10f \n',score2,causal_effect2);
            fprintf('score3: %0.10f \n',score3);
        end
        
        
        if round(score1,10) < round(score2,10) && round(score1,10) < round(score3,10)
            if construct_verbose
                fprintf('consider   %d ->  %d \n',from,to);
                fprintf('score1: %0.10f causal effect1=%0.10f \n',score1);
                fprintf('score2: %0.10f causal effect2=%0.10f \n',score2);
                fprintf('score3: %0.10f \n',score3);
                
            end
            whichPair(iPair,4)= 1;
            whichPair(iPair,5)= score1;
            
            
            
        elseif round(score2,10) < round(score1,10) && round(score2,10) < round(score3,10)
            if construct_verbose
                fprintf('consider   %d ->  %d \n',to,from);
                fprintf('score1: %0.10f causal effect1=%0.10f \n',score1);
                fprintf('score2: %0.10f causal effect2=%0.10f \n',score2);
                fprintf('score3: %0.10f \n',score3);
            end
            whichPair(iPair,4)= 2;
            whichPair(iPair,5)= score2;
            
            
        elseif round(score3,10) < round(score2,10) && round(score3,10) < round(score1,10)
            if construct_verbose
                fprintf('consider   %d <->  %d \n',to,from);
                fprintf('score1: %0.15f causal effect1=%0.10f \n',score1,causal_effect1);
                fprintf('score2: %0.15f causal effect2=%0.10f \n',score2,causal_effect2);
                fprintf('score3: %0.15f \n',score3);
            end
            
            whichPair(iPair,4)= 3;
            whichPair(iPair,5)= score3;
            
            
            
        else
            %if BIC scores for ->, <-, <-> are equal use direct causal effect
            %to identify direction
            
            
            %if BIC is not decreased, return undirect graph
            if score1==Inf && score2==Inf && score3==Inf
                whichPair(iPair,3)= 0;
                
            else
                
                if causal_effect2 < causal_effect1
                    if construct_verbose
                        fprintf('consider   %d ->  %d as score1=%d score2=%d score3=%d\n',from,to,score1,score2,score3);
                        fprintf('causal effect =%d  %d ->  %d \n',causal_effect1,from,to);
                        fprintf('causal effect =%d  %d <-  %d \n',causal_effect2,from,to);
                    end
                    whichPair(iPair,5)= score1;
                    
                    whichPair(iPair,4)= 1;
                    
                elseif causal_effect1 < causal_effect2
                    if construct_verbose
                        fprintf('consider   %d ->  %d score1=%d score2=%d score3=%d\n',to,from,score1,score2,score3);
                        fprintf('causal effect =%d  %d ->  %d \n',causal_effect2,to,from);
                        fprintf('causal effect =%d  %d <-  %d \n',causal_effect1,to,from);
                    end
                    whichPair(iPair,5)= score2;
                    whichPair(iPair,4)= 2;
                    
                else
                    %if BIC scores are equal and direct causal effect are equal, return undirect graph
                    whichPair(iPair,3)=0;
                    if construct_verbose
                        fprintf('undirected from %d to  %d \n',from,to);
                        fprintf('score1: %0.10f causal effect1=%0.10f \n',score1,causal_effect1);
                        fprintf('score2: %0.10f causal effect2=%0.10f \n',score2,causal_effect2);
                        fprintf('score3: %0.10f \n',score3);
                    end
                    %
                end
            end
        end
    end
    
end


[all,~] = max(whichPair(:,3));

while (all ==1)
    
    [M,I] = min(whichPair(:,5));
    
    %if there is undired edges but non-collider block the orientation, all
    %pair score will be infinity
    
    if M>=CurScore
        
        %find *A-from-to if no prohibit action from -> to, try A<-from->to
        
        pair= find(whichPair(:,1)==from & whichPair(:,2)==to);
        
        
        if (pair>0 & whichPair(pair,6)==0 & whichPair(pair,7)==1 & whichPair(pair,8)==1)
           
            nEdges=nEdges+1;
            isParent(from,to)=true;
            [score1, scoreContribs1, newNComps1, newSizes1, newComps1, newInComponent1, mag1] = addDirectedEdge(from, to, nEdges, newNComps, newSizes, newComps, newInComponent, IterMag{iter-1}, covMat, isParent, scoreContribs, nSamples, nVars,  tol);
            
            
            %if from-> to creates cyclic try to reverse A->from->to to
            %A<-from->to
            if ~isAcyclic(mag1)
                CurScore=score1;
                scoreContribs=scoreContribs1;
                newNComps=newNComps1;
                newSizes=newSizes1;
                newComps=newComps1;
                newInComponent=newInComponent1;
                
                % find the matrix [1xN] for the list of direct path
                % of cyclic
                [cyclic,ncyclic]  = findCyclic(mag1);
                mat_cyclic=cyclic{1,1};
                [i_index,~]= find(mag1(:,from)==2 & mag1(from,:)==3);
                a=intersect(unique(i_index)',mat_cyclic);
              
                %more than one cyclic?
                
                if (ncyclic>1)
                    
                   
                    nEdges=nEdges+1;
                    isParent(to,from)=true;
                    [score1, scoreContribs1, newNComps1, newSizes1, newComps1, newInComponent1, mag1] = addDirectedEdge(to, from, nEdges, newNComps, newSizes, newComps, newInComponent, IterMag{iter-1}, covMat, isParent, scoreContribs, nSamples, nVars,  tol);
                    
                    
                    
                    
                    CurScore=score1;
                    scoreContribs=scoreContribs1;
                    newNComps=newNComps1;
                    newSizes=newSizes1;
                    newComps=newComps1;
                    newInComponent=newInComponent1;
                    IterMag{iter}=mag1;
                    iter=iter+1;
                    whichPair(pair,3)=0;
                    whichPair(:,4)=0;
                    whichPair(:,5)=Inf;
                    
                    
                else
                    
                    %reverse a->from
                    isParent(a,from)=false;
                    nEdges=nEdges-1;
                    [CurScore2, scoreContribs2, newNComps2, newSizes2, newComps2, newInComponent2, mag2,beta2,omega2,district2] = removeDirectedEdge(a, from, nEdges,newNComps, newSizes, newComps, newInComponent,mag1, covMat, isParent, scoreContribs, nSamples, nVars,  tol);
                    
                    isParent(from,a)=true;
                    nEdges=nEdges+1;
                    [CurScoreA, scoreContribsA, newNCompsA, newSizesA, newCompsA, newInComponentA, magA,betaA,omegaA,districtA] = addDirectedEdge(from, a, nEdges,newNComps2, newSizes2, newComps2, newInComponent2,mag2, covMat, isParent, scoreContribs2, nSamples, nVars,  tol);
                    
                    isParent(a,from)=true;
                    isParent(from,a)=false;
                    
                    %end revers a-> from section
                    
                    
                    %reverse to->c
                    [i_index,~]= find(mag1(to,:)==2 & mag1(:,to)==3);
                    c=intersect(unique(i_index)',mat_cyclic);
                    isParent(to,c)=false;
                    nEdges=nEdges-1;
                    [CurScore3, scoreContribs3, newNComps3, newSizes3, newComps3, newInComponent3,  mag3,beta3,omega3,district3] = removeDirectedEdge(to, c, nEdges,newNComps, newSizes, newComps, newInComponent,mag1, covMat, isParent, scoreContribs, nSamples, nVars,  tol);
                    
                    isParent(c,to)=true;
                    nEdges=nEdges+1;
                    [CurScoreB, scoreContribsB, newNCompsB, newSizesB, newCompsB, newInComponentB, magB,betaB,omegaB,districtB] = addDirectedEdge(c, to, nEdges,newNComps3, newSizes3, newComps3, newInComponent3,mag3, covMat, isParent, scoreContribs3, nSamples, nVars,  tol);
                    %causal_effectB=abs(beta_matB(from,a));
                    isParent(to,c)=true;
                    isParent(c,to)=false;
                    %end reverse to->c section
                    
                    if (CurScoreA<CurScoreB & isAcyclic(magA))
                       
                        isParent(a,from)=false;
                        isParent(from,a)=true;
                        
                        CurScore=CurScoreA;
                        scoreContribs=scoreContribsA;
                        newNComps=newNCompsA;
                        newSizes=newSizesA;
                        newComps=newCompsA;
                        newInComponent=newInComponentA;
                        beta=betaA;
                        omega=omegaA;
                        district=districtA;
                        IterMag{iter}=magA;
                        iter=iter+1;
                        whichPair(pair,3)=0;
                        whichPair(:,4)=0;
                        whichPair(:,5)=Inf;
                        
                    elseif (CurScoreB<CurScoreA & isAcyclic(magB))
                       
                        isParent(to,c)=false;
                        isParent(c,to)=true;
                        %assignin('base','xmagB',magB);
                        CurScore=CurScoreB;
                        scoreContribs=scoreContribsB;
                        newNComps=newNCompsB;
                        newSizes=newSizesB;
                        newComps=newCompsB;
                        newInComponent=newInComponentB;
                        beta=betaB;
                        omega=omegaB;
                        district=districtB;
                        IterMag{iter}=magB;
                        iter=iter+1;
                        whichPair(pair,3)=0;
                        whichPair(:,4)=0;
                        whichPair(:,5)=Inf;
                        
                    else
                      
                        
                        %                           CurScore=score1;
                        %                           scoreContribs=scoreContribs1;
                        %                           newNComps=newNComps1;
                        %                           newSizes=newSizes1;
                        %                           newComps=newComps1;
                        %                           newInComponent=newInComponent1;
                        %                           IterMag{iter}=mag1;
                        %                           iter=iter+1;
                        %                           whichPair(pair,3)=0;
                        %                           whichPair(:,4)=0;
                        %                           whichPair(:,5)=Inf;
                    end
                end
            else
                
                CurScore=score1;
                scoreContribs=scoreContribs1;
                newNComps=newNComps1;
                newSizes=newSizes1;
                newComps=newComps1;
                newInComponent=newInComponent1;
                IterMag{iter}=mag1;
                iter=iter+1;
                whichPair(pair,3)=0;
                whichPair(:,4)=0;
                whichPair(:,5)=Inf;
                
            end
            
        end
        
        
        %find *A-from<-to if no prohibit action from<-to but prohibit from
        %-> and <->
        
        if (pair>0 & whichPair(pair,6)==1 & whichPair(pair,7)==0 & whichPair(pair,8)==1)
          
            nEdges=nEdges+1;
            isParent(to,from)=true;
            [score1, scoreContribs1, newNComps1, newSizes1, newComps1, newInComponent1, mag1] = addDirectedEdge(to, from, nEdges, newNComps, newSizes, newComps, newInComponent, IterMag{iter-1}, covMat, isParent, scoreContribs, nSamples, nVars,  tol);
            
            %if from<- to creates cyclic try to reverse A<-from<-to to
            %A->from<-to
            if ~isAcyclic(mag1)
                CurScore=score1;
                scoreContribs=scoreContribs1;
                newNComps=newNComps1;
                newSizes=newSizes1;
                newComps=newComps1;
                newInComponent=newInComponent1;
                
                % find the matrix [1xN] for the list of direct path
                % of cyclic
                [cyclic,ncyclic] = findCyclic(mag1);
                
                %assignin('base', 'xcyclic', cyclic);
                
                mat_cyclic=cyclic{1,1};
                
                %find a<-from
                
                [i_index,~]= find(mag1(from,:)==2 & mag1(:,from)==3);
                
                a=intersect(unique(i_index)',mat_cyclic);
                
                %more than one cyclic?
                %assignin('base', 'xa', a);
                if (ncyclic>1)
                    
                  
                    nEdges=nEdges+1;
                    isParent(from,to)=true;
                    [score1, scoreContribs1, newNComps1, newSizes1, newComps1, newInComponent1, mag1] = addDirectedEdge(from, to, nEdges, newNComps, newSizes, newComps, newInComponent, IterMag{iter-1}, covMat, isParent, scoreContribs, nSamples, nVars,  tol);
                    
                    CurScore=score1;
                    scoreContribs=scoreContribs1;
                    newNComps=newNComps1;
                    newSizes=newSizes1;
                    newComps=newComps1;
                    newInComponent=newInComponent1;
                    IterMag{iter}=mag1;
                    iter=iter+1;
                    whichPair(pair,3)=0;
                    whichPair(:,4)=0;
                    whichPair(pair,5)=Inf;
                    
                else
                    %%reverse a<-from
                    isParent(from,a)=false;
                    nEdges=nEdges-1;
                    [CurScore2, scoreContribs2, newNComps2, newSizes2, newComps2, newInComponent2, mag2,beta2,omega2,district2] = removeDirectedEdge(a, from, nEdges,newNComps, newSizes, newComps, newInComponent,mag1, covMat, isParent, scoreContribs, nSamples, nVars,  tol);
                    
                    isParent(a,from)=true;
                    nEdges=nEdges+1;
                    [CurScoreA, scoreContribsA, newNCompsA, newSizesA, newCompsA, newInComponentA, magA,betaA,omegaA,districtA] = addDirectedEdge(a, from, nEdges,newNComps2, newSizes2, newComps2, newInComponent2,mag2, covMat, isParent, scoreContribs2, nSamples, nVars,  tol);
                    
                    %restore isParent
                    isParent(from,a)=true;
                    isParent(a,from)=false;
                    
                    %%end revers a-> from section
                    
                    
                    %find the score to reverse to<-c
                    [i_index,~]= find(mag1(:,to)==2 & mag1(to,:)==3);
                    c=intersect(unique(i_index)',mat_cyclic);
                    if isempty(a)
                        
                       
                        nEdges=nEdges+1;
                        isParent(to,from)=true;
                        [score1, scoreContribs1, newNComps1, newSizes1, newComps1, newInComponent1, mag1] = addDirectedEdge(to, from, nEdges, newNComps, newSizes, newComps, newInComponent, IterMag{iter-1}, covMat, isParent, scoreContribs, nSamples, nVars,  tol);
                        CurScore=score1;
                        scoreContribs=scoreContribs1;
                        newNComps=newNComps1;
                        newSizes=newSizes1;
                        newComps=newComps1;
                        newInComponent=newInComponent1;
                        IterMag{iter}=mag1;
                        iter=iter+1;
                        whichPair(pair,3)=0;
                        whichPair(:,4)=0;
                        whichPair(:,4)=score1;
                    else
                        
                        isParent(c,to)=false;
                        nEdges=nEdges-1;
                        [CurScore3, scoreContribs3, newNComps3, newSizes3, newComps3, newInComponent3,  mag3,beta3,omega3,district3] = removeDirectedEdge(c, to, nEdges,newNComps, newSizes, newComps, newInComponent,mag1, covMat, isParent, scoreContribs, nSamples, nVars,  tol);
                        
                        isParent(to,c)=true;
                        nEdges=nEdges+1;
                        [CurScoreB, scoreContribsB, newNCompsB, newSizesB, newCompsB, newInComponentB, magB,betaB,omegaB,districtB] = addDirectedEdge(to, c, nEdges,newNComps3, newSizes3, newComps3, newInComponent3,mag3, covMat, isParent, scoreContribs3, nSamples, nVars,  tol);
                        %causal_effectB=abs(beta_matB(from,a));
                        isParent(c,to)=true;
                        isParent(to,c)=false;
                        %end reverse to->c section
                        
                        if (CurScoreA<CurScoreB & isAcyclic(magA))
                          
                            isParent(from,a)=false;
                            isParent(a,from)=true;
                            
                            CurScore=CurScoreA;
                            scoreContribs=scoreContribsA;
                            newNComps=newNCompsA;
                            newSizes=newSizesA;
                            newComps=newCompsA;
                            newInComponent=newInComponentA;
                            beta=betaA;
                            omega=omegaA;
                            district=districtA;
                            IterMag{iter}=magA;
                            iter=iter+1;
                            whichPair(pair,3)=0;
                            whichPair(:,4)=0;
                            whichPair(:,4)=CurScoreA;
                            
                        elseif (CurScoreB<CurScoreA & isAcyclic(magB))
                           
                            isParent(c,to)=false;
                            isParent(to,c)=true;
                            %assignin('base','xmagB',magB);
                            CurScore=CurScoreB;
                            scoreContribs=scoreContribsB;
                            newNComps=newNCompsB;
                            newSizes=newSizesB;
                            newComps=newCompsB;
                            newInComponent=newInComponentB;
                            beta=betaB;
                            omega=omegaB;
                            district=districtB;
                            IterMag{iter}=magB;
                            iter=iter+1;
                            whichPair(pair,3)=0;
                            whichPair(:,4)=0;
                            whichPair(:,4)=CurScoreB;
                        else
                            %fprintf('score equivalent  \n');
  
                          
                        end
                    end
                    
                end
            elseif ~isAcyclic(mag1) & score1 <= CurScore
                CurScore=score1;
                scoreContribs=scoreContribs1;
                newNComps=newNComps1;
                newSizes=newSizes1;
                newComps=newComps1;
                newInComponent=newInComponent1;
                IterMag{iter}=mag1;
                iter=iter+1;
                whichPair(pair,3)=0;
                whichPair(:,4)=0;
                whichPair(pair,5)=Inf;
                
            else
                nEdges=nEdges-1;
                isParent(to,from)=false;
               
                
            end
        end
        
        
        %find -> , <- or  <->
        if (pair>0 & ((whichPair(pair,6)==0 & whichPair(pair,7)==0 & whichPair(pair,8)==0) | (whichPair(pair,6)==1 & whichPair(pair,7)==1 & whichPair(pair,8)==1) ))
           
            
            nEdges=nEdges+1;
            isParent(from,to)=true;
            [score1, scoreContribs1, newNComps1, newSizes1, newComps1, newInComponent1, mag1,beta1,omega1,district1] = addDirectedEdge(from, to, nEdges, newNComps, newSizes, newComps, newInComponent, IterMag{iter-1}, covMat, isParent, scoreContribs, nSamples, nVars,  tol);
            isParent(from,to)=false;
            
            omega_mat1(district1,district1)= omega1;
            beta_mat1(district1,district1)= beta1;
            causal_effect1=abs(beta_mat1(to,from));
            
            if ~isAcyclic(mag1)
                score1=Inf;
                causal_effect1 =-Inf;
                
            end
            
            isParent(to,from)=true;
            [score2, scoreContribs2, newNComps2, newSizes2, newComps2, newInComponent2, mag2,beta2,omega2,district2] = addDirectedEdge(to, from, nEdges, newNComps, newSizes, newComps, newInComponent, IterMag{iter-1}, covMat, isParent, scoreContribs, nSamples, nVars,  tol);
            isParent(to,from)=false;
            omega_mat2(district2,district2)= omega2;
            beta_mat2(district2,district2)= beta2;
            causal_effect2=abs(beta_mat2(from,to));
            
            if ~isAcyclic(mag2)
                score2=Inf;
                causal_effect2=-Inf;
            end
            
            [score3, scoreContribs3, newNComps3, newSizes3, newComps3, newInComponent3, mag3,~,~,~] = addBidirectedEdge(to, from, nEdges, newNComps, newSizes, newComps, newInComponent, IterMag{iter-1}, covMat, isParent, scoreContribs, nSamples, nVars,  tol);
            nEdges=nEdges-1;
            
            
            
            if round(score1,10) < round(score2,10) && round(score1,10) < round(score3,10)
              
                isParent(from,to)=true;
                nEdges=nEdges+1;
                CurScore=score1;
                scoreContribs=scoreContribs1;
                newNComps=newNComps1;
                newSizes=newSizes1;
                newComps=newComps1;
                newInComponent=newInComponent1;
                IterMag{iter}=mag1;
                iter=iter+1;
                whichPair(pair,3)=0;
                whichPair(:,4)=0;
                whichPair(pair,5)=Inf;
            elseif round(score2,10) < round(score1,10) && round(score2,10) < round(score3,10)
               
                nEdges=nEdges+1;
                isParent(to,from)=true;
                CurScore=score2;
                scoreContribs=scoreContribs2;
                newNComps=newNComps2;
                newSizes=newSizes2;
                newComps=newComps2;
                newInComponent=newInComponent2;
                IterMag{iter}=mag2;
                iter=iter+1;
                whichPair(pair,3)=0;
                whichPair(:,4)=0;
                whichPair(pair,5)=Inf;
            elseif round(score3,10) < round(score2,10) && round(score3,10) < round(score1,10)
               
                nEdges=nEdges+1;
                CurScore=score3;
                scoreContribs=scoreContribs3;
                newNComps=newNComps3;
                newSizes=newSizes3;
                newComps=newComps3;
                newInComponent=newInComponent3;
                IterMag{iter}=mag3;
                iter=iter+1;
                whichPair(pair,3)=0;
                whichPair(:,4)=0;
                whichPair(pair,3)=score3;
            elseif (score3<=CurScore)
                
               
                
                if causal_effect2 < causal_effect1
                    if construct_verbose
                        fprintf('consider %d ->%d \n',from,to);
                    end
                  
                    isParent(from,to)=true;
                    CurScore=score1;
                    scoreContribs=scoreContribs1;
                    newNComps=newNComps1;
                    newSizes=newSizes1;
                    newComps=newComps1;
                    newInComponent=newInComponent1;
                    IterMag{iter}=mag1;
                    iter=iter+1;
                    
                    whichPair(pair,3)=0;
                    whichPair(:,4)=0;
                    whichPair(pair,5)=Inf;
                  
                elseif causal_effect1 < causal_effect2
                    if construct_verbose
                        fprintf('consider  %d->%d \n',to,from);
                    end
                  
                    
                    isParent(to,from)=true;
                    CurScore=score2;
                    scoreContribs=scoreContribs2;
                    newNComps=newNComps2;
                    newSizes=newSizes2;
                    newComps=newComps2;
                    newInComponent=newInComponent2;
                    IterMag{iter}=mag2;
                    
                    iter=iter+1;
                    whichPair(pair,3)=0;
                    whichPair(:,4)=0;
                    whichPair(pair,5)=Inf;
                   
                else
                    
                    pag(to,from)=0;
                    pag(from,to)=0;
                    whichPair(pair,3)=0;
                    whichPair(:,4)=0;
                    whichPair(pair,5)=Inf;
                end
            else
                
              
                %remove undirected edge from the skeleton
                pag(to,from)=0;
                pag(from,to)=0;
                whichPair(pair,3)=0;
                whichPair(:,4)=0;
                whichPair(pair,5)=Inf;
                
            end
            
        else
           
            %fprintf('pair=%d whichPair(6)=%d whichPair(7)=%d whichPair(8)=%d  \n',pair,whichPair(pair,6),whichPair(pair,7),whichPair(pair,8));
            pag(to,from)=0;
            pag(from,to)=0;
            whichPair(pair,3)=0;
            whichPair(:,4)=0;
            whichPair(pair,5)=Inf;
            
        end
        
    end
    
    
    from = whichPair(I,1); to = whichPair(I,2);
    
    if (whichPair(I,4)==1)
        if construct_verbose
            fprintf('Added %d -> %d \n',from,to);
        end
        
        nEdges=nEdges+1;
        isParent(from,to)=true;
        if iter==1
            
            [CurScore, scoreContribs, newNComps, newSizes, newComps, newInComponent, IterMag{iter}] = addDirectedEdge(from, to, nEdges, nComps, sizes, comps, inComponent, mag, covMat, isParent, scores, nSamples, nVars,  tol);
            
        else
            
            [CurScore, scoreContribs, newNComps, newSizes, newComps, newInComponent, IterMag{iter}] = addDirectedEdge(from, to, nEdges, newNComps, newSizes, newComps, newInComponent, IterMag{iter-1}, covMat, isParent, scoreContribs, nSamples, nVars,  tol);
            
        end
        
        gs(iter).score=CurScore;
        
        
        iter=iter+1;
        
        whichPair(I,3)=0;
        whichPair(:,4)=0;
        whichPair(:,5)=Inf;
        if size(dnc) ~=0
            
            [row1,~] = find(dnc(:,1)==from & dnc(:,2)==to);
            if size(row1,1)>0
                dnc_index =size(row1,1);
                for i =1:dnc_index
                    [row2,~]= find(whichPair(:,1)==dnc(row1(dnc_index,1),2) & whichPair(:,2)==dnc(row1(dnc_index,1),3));
                    
                    if size(row2,1)>0
                        whichPair(row2,7)=1;
                        whichPair(row2,8)=1;
                        %                                   fprintf('avoid dnc= %d->%d-*%d',whichPair(row2,1),whichPair(row2,2));
                        %
                    end
                    
                    
                    [row2,~]= find(whichPair(:,2)==dnc(row1(dnc_index,1),2) & whichPair(:,1)==dnc(row1(dnc_index,1),3));
                    
                    if size(row2,1)>0
                        whichPair(row2,6)=1;
                        whichPair(row2,8)=1;
                        %                                    fprintf('avoid dnc= %d->%d-*%d',whichPair(row2,1),whichPair(row2,2));
                        %
                    end
                end
                
                
            end
            
            [row1,~] = find(dnc(:,3)==from & dnc(:,2)==to);
            if size(row1,1)>0
                dnc_index =size(row1,1);
                for i =1:dnc_index
                    [row2,~]= find(whichPair(:,1)==dnc(row1(dnc_index,1),1) & whichPair(:,2)==dnc(row1(dnc_index,1),2));
                    
                    if size(row2,1)>0
                        whichPair(row2,6)=1;
                        whichPair(row2,8)=1;
                        %                                  fprintf('avoid dnc= %d->%d-*%d',whichPair(row2,1),whichPair(row2,2));
                        %
                    end
                    
                    [row2,~]= find(whichPair(:,2)==dnc(row1(dnc_index,1),1) & whichPair(:,1)==dnc(row1(dnc_index,1),2));
                    
                    if size(row2,1)>0
                        whichPair(row2,7)=1;
                        whichPair(row2,8)=1;
                        %                                    fprintf('avoid dnc= %d->%d-*%d',whichPair(row2,1),whichPair(row2,2));
                        %
                    end
                end
                
                
            end
        end
        
    elseif (whichPair(I,4)==2)
        if construct_verbose
            fprintf('Added %d -> %d \n',to,from);
        end
        nEdges=nEdges+1;
        isParent(to,from)=true;
        if iter==1
            
            [CurScore, scoreContribs, newNComps, newSizes, newComps, newInComponent, IterMag{iter}] = addDirectedEdge(to, from, nEdges, nComps, sizes, comps, inComponent, mag, covMat, isParent, scores, nSamples, nVars,  tol);
            
        else
            
            [CurScore, scoreContribs, newNComps, newSizes, newComps, newInComponent, IterMag{iter}] = addDirectedEdge(to, from, nEdges, newNComps, newSizes, newComps, newInComponent, IterMag{iter-1}, covMat, isParent, scoreContribs, nSamples, nVars,  tol);
            %                      fprintf('%d -> %d \n',to,from);
            
            
        end
        
        gs(iter).score=CurScore;
        
        iter=iter+1;
        
        
        whichPair(I,3)=0;
        whichPair(:,4)=0;
        whichPair(:,5)=Inf;
        if size(dnc) ~=0
            
            [row1,~] = find(dnc(:,1)==to & dnc(:,2)==from);
            if size(row1,1)>0
                dnc_index =size(row1,1);
                for i =1:dnc_index
                    [row2,~]= find(whichPair(:,1)==dnc(row1(dnc_index,1),2) & whichPair(:,2)==dnc(row1(dnc_index,1),3));
                    
                    if size(row2,1)>0
                        whichPair(row2,7)=1;
                        whichPair(row2,8)=1;
                        %                                   fprintf('avoid dnc= %d->%d-*%d',whichPair(row2,1),whichPair(row2,2));
                        %
                    end
                    
                    
                    
                    [row2,~]= find(whichPair(:,2)==dnc(row1(dnc_index,1),2) & whichPair(:,1)==dnc(row1(dnc_index,1),3));
                    
                    if size(row2,1)>0
                        whichPair(row2,6)=1;
                        whichPair(row2,8)=1;
                        %                                    fprintf('avoid dnc= %d->%d-*%d',whichPair(row2,1),whichPair(row2,2));
                        %
                    end
                    
                end
                
                
            end
            
            [row1,~] = find(dnc(:,3)==to & dnc(:,2)==from);
            if size(row1,1)>0
                dnc_index =size(row1,1);
                for i =1:dnc_index
                    [row2,~]= find((whichPair(:,1)==dnc(row1(dnc_index,1),1) & whichPair(:,2)==dnc(row1(dnc_index,1),2)));
                    
                    if size(row2,1)>0
                        whichPair(row2,6)=1;
                        whichPair(row2,8)=1;
                        %                                    fprintf('avoid dnc= %d->%d-*%d',whichPair(row2,1),whichPair(row2,2));
                        %
                    end
                    
                    [row2,~]= find((whichPair(:,2)==dnc(row1(dnc_index,1),1) & whichPair(:,1)==dnc(row1(dnc_index,1),2)));
                    
                    if size(row2,1)>0
                        whichPair(row2,7)=1;
                        whichPair(row2,8)=1;
                        %                                    fprintf('avoid dnc= %d->%d-*%d',whichPair(row2,1),whichPair(row2,2));
                        %
                    end
                end
            end
            
        end
        
        
        
        
        
    elseif (whichPair(I,4)==3)
        
        previous_score=CurScore;
        nEdges=nEdges+1;
        if iter==1
            [CurScore, scoreContribs, newNComps, newSizes, newComps, newInComponent, IterMag{iter},beta,omega,district] = addBidirectedEdge(from, to, nEdges, nComps, sizes, comps, inComponent, mag, covMat, isParent, scores, nSamples, nVars,  tol);
            if construct_verbose
                fprintf('added    %d <->  %d \n',from,to);
            end
        else
            [CurScore, scoreContribs, newNComps, newSizes, newComps, newInComponent, IterMag{iter},beta,omega,district] = addBidirectedEdge(from, to, nEdges, newNComps, newSizes, newComps, newInComponent, IterMag{iter-1}, covMat, isParent, scoreContribs, nSamples, nVars,  tol);
            
        end
        if construct_verbose
            fprintf('Added %d <-> %d score=%.4f \n',from ,to,(previous_score-CurScore));
        end
        gs(iter).score=CurScore;
        
        iter=iter+1;
        
        omega_mat(district,district)=omega;
        beta_mat(district,district)=beta;
        %                     omega_mat
        %                     beta_mat
        whichPair(I,3)=0;
        whichPair(:,4)=0;
        whichPair(:,5)=Inf;
        if size(dnc) ~=0
            [row1,~] = find(dnc(:,1)==from & dnc(:,2)==to);
            if size(row1,1)>0
                dnc_index =size(row1,1);
                for i =1:dnc_index
                    [row2,~]= find(whichPair(:,1)==dnc(row1(dnc_index,1),2) & whichPair(:,2)==dnc(row1(dnc_index,1),3));
                    
                    if size(row2,1)>0
                        whichPair(row2,7)=1;
                        whichPair(row2,8)=1;
                        %                                    fprintf('avoid dnc= %d->%d-*%d',whichPair(row2,1),whichPair(row2,2));
                        %
                    end
                    
                    
                    [row2,~]= find(whichPair(:,2)==dnc(row1(dnc_index,1),2) & whichPair(:,1)==dnc(row1(dnc_index,1),3));
                    
                    if size(row2,1)>0
                        whichPair(row2,6)=1;
                        whichPair(row2,8)=1;
                        %                                    fprintf('avoid dnc= %d->%d-*%d',whichPair(row2,1),whichPair(row2,2));
                        %
                    end
                end
                
                
            end
            
            [row1,~] = find(dnc(:,3)==from & dnc(:,2)==to);
            if size(row1,1)>0
                dnc_index =size(row1,1);
                for i =1:dnc_index
                    [row2,~]= find(whichPair(:,1)==dnc(row1(dnc_index,1),1) & whichPair(:,2)==dnc(row1(dnc_index,1),2));
                    
                    if size(row2,1)>0
                        whichPair(row2,6)=1;
                        whichPair(row2,8)=1;
                        %                                    fprintf('avoid dnc= %d->%d-*%d',whichPair(row2,1),whichPair(row2,2));
                        %
                    end
                    
                    [row2,~]= find(whichPair(:,2)==dnc(row1(dnc_index,1),1) & whichPair(:,1)==dnc(row1(dnc_index,1),2));
                    
                    if size(row2,1)>0
                        whichPair(row2,7)=1;
                        whichPair(row2,8)=1;
                        %                                    fprintf('avoid dnc= %d->%d-*%d',whichPair(row2,1),whichPair(row2,2));
                        %
                    end
                end
                
                
            end
            
        end
        
        
    end
    
    %find undirected edges
    for iPair =1:nPairs
        
        from = whichPair(iPair, 1); to = whichPair(iPair, 2);
        
        if pag(from, to)==3 && pag(to, from)==3 && IterMag{iter-1}(from, to)==0 &&IterMag{iter-1}(to, from)==0
            whichPair(iPair,3)= 1;
            %               fprintf('undirected %d-%d  \n' ,whichPair(iPair,1) ,whichPair(iPair,2));
            %               assignin('base','whichPair2',whichPair);
            
        end
    end
    
    
    [row,~] = find(whichPair(:,3)==1);
    
    for v=row'
        
        from = whichPair(v,1); to = whichPair(v,2);
        if whichPair(v,7)==1 && whichPair(v,8)==1
            nEdges=nEdges+1;
            isParent(from,to)=true;
            [score1,~,~,~,~,~,Mag_A,beta1,omega1,district] = addDirectedEdge(from, to, nEdges, newNComps, newSizes, newComps, newInComponent, IterMag{iter-1}, covMat, isParent, scoreContribs, nSamples, nVars,  tol);
            isParent(from,to)=false;
            nEdges=nEdges-1;
            score2=Inf;
            score3=Inf;
            omega_mat1(district,district)= omega1;
            beta_mat1(district,district)= beta1;
            
            causal_effect1=abs(beta_mat1(to,from));
            
            if ~isAcyclic(Mag_A)
                score1 =Inf;
            end
        elseif whichPair(v,6)==1 && whichPair(v,8)==1
            
            score1=Inf;
            nEdges=nEdges+1;
            isParent(to,from)=true;
            [score2,~,~,~,~,~,Mag_B,beta2,omega2,district] = addDirectedEdge(to, from, nEdges, newNComps, newSizes, newComps, newInComponent, IterMag{iter-1}, covMat, isParent, scoreContribs, nSamples, nVars,  tol);
            nEdges=nEdges-1;
            isParent(to,from)=false;
            score3=Inf;
            
            omega_mat2(district,district)= omega2;
            beta_mat2(district,district)= beta2;
            
            causal_effect2=abs(beta_mat2(from,to));
            
            if ~isAcyclic(Mag_B)
                score2 =Inf;
            end
            
            
        else
            nEdges=nEdges+1;
            
            isParent(from,to)=true;
            [score1,~,~,~,~,~,Mag_A,beta1,omega1,district] = addDirectedEdge(from, to, nEdges, newNComps, newSizes, newComps, newInComponent, IterMag{iter-1}, covMat, isParent, scoreContribs, nSamples, nVars,  tol);
            isParent(from,to)=false;
            omega_mat1(district,district)= omega1;
            beta_mat1(district,district)= beta1;
            causal_effect1=abs(beta_mat1(to,from));
          
            if ~isAcyclic(Mag_A)
                score1 =Inf;
            end
            
            
            
            isParent(to,from)=true;
            [score2,~,~,~,~,~,Mag_B,beta2,omega2,district] = addDirectedEdge(to, from, nEdges, newNComps, newSizes, newComps, newInComponent, IterMag{iter-1}, covMat, isParent, scoreContribs, nSamples, nVars,  tol);
            isParent(to,from)=false;
            omega_mat2(district,district)= omega2;
            beta_mat2(district,district)=beta2;
            causal_effect2=abs(beta_mat2(from,to));
        
            
            if ~isAcyclic(Mag_B)
                score2 =Inf;
            end
            
            
            
            [score3,~,~,~,~,~,~,beta3,omega3,district] = addBidirectedEdge(from, to, nEdges, newNComps, newSizes, newComps, newInComponent, IterMag{iter-1}, covMat, isParent, scoreContribs, nSamples, nVars,  tol);
            omega_mat3(district,district)= omega3;
            beta_mat3(district,district)=beta3;
          
            nEdges=nEdges-1;
            if construct_verbose
                fprintf('iter:%d from  %d to %d \n',iter,from,to);
                fprintf('Current Score: %0.10f \n',CurScore);
                fprintf('score1: %0.10f causal effect1=%0.10f \n',score1,causal_effect1);
                fprintf('score2: %0.10f causal effect2=%0.10f \n',score2,causal_effect2);
                fprintf('score3: %0.10f \n',score3);
            end
            
        end
        
        
        if round(score1,10) < round(score2,10) && round(score1,10) < round(score3,10)
            whichPair(v,4)= 1;
            whichPair(v,5)= score1;
        elseif round(score2,10) < round(score1,10) && round(score2,10) < round(score3,10)
            whichPair(v,4)= 2;
            whichPair(v,5)= score2;
        elseif round(score3,10) < round(score2,10) && round(score3,10) < round(score1,10)
            whichPair(v,4)= 3;
            whichPair(v,5)= score3;
        else
            if causal_effect2 < causal_effect1
                if construct_verbose
                    fprintf('consider   %d ->  %d \n',from,to);
                end
                whichPair(v,4)= 1;
                whichPair(v,5)= score1;
                
            elseif causal_effect1 < causal_effect2
                if construct_verbose
                    fprintf('consider   %d ->  %d \n',to,from);
                end
                whichPair(v,4)= 2;
                whichPair(v,5)= score2;
                
            else
                whichPair(v,3)=0;
            end
            
            
        end
        
    end
    
    [all,~] = max(whichPair(:,3));
    
    
    if iter >1
        
        gs(iter).score=CurScore;
        
    end
    
    
   
end




if iter>1
    mag=IterMag{iter-1};
    
end



%Change Mag o- to -
nvars = size(FCI_skeleton,1); % KC- Return number of variable
for i = 1:nvars
    for j=1:nvars
        if(FCI_skeleton(i,j)==3 && mag(i,j)==0 )
            mag(i,j) = 1;
        end
    end
end




end


function [newScore, scoreContribs, newNComps, newSizes, newComps, newInComponent, newMag,beta,omega,district] = addDirectedEdge(from, to, nEdges, nComps, sizes, comps, inComponent, mag, covMat, isParent, scores, nSamples, nVars,  tol)

% if the edge is bidirected you first have to update the components.
if mag(to, from)==2 & mag(from, to)==2
    [~, scores, newNComps, newSizes, newComps, newInComponent, newMag] = removeBidirectedEdge(from, to, nEdges, nComps, sizes, comps, inComponent, mag, covMat, isParent, scores, nSamples, nVars,  tol);
else
    [newNComps, newSizes, newComps, newInComponent, newMag] = deal(nComps, sizes, comps, inComponent, mag);
end
% change mag
newMag(from, to)=2;newMag(to, from)=3;
tmpIsParent = isParent;
tmpIsParent(from, to)=true;
% if the edge was directed in reverse, you must also update
% inComponent(from);
scoreContribs = scores;

if mag(from, to)==3;
    tmpIsParent(to, from)=false;
    component= comps{inComponent(from)};
    [compMag, district] = componentMag(component, nVars, newMag, tmpIsParent);
    [scoreContribs(inComponent(from)),beta, omega] = score_contrib(compMag, component, district, sizes(inComponent(from)), covMat, nSamples, tol);
end
% update district
component = comps{inComponent(to)};
[compMag, district] = componentMag(component, nVars, newMag, tmpIsParent);
% get new score
[scoreContribs(inComponent(to)),beta, omega] = score_contrib(compMag, component, district, sizes(inComponent(to)), covMat, nSamples, tol);

tmpSll = (-nSamples/2)*sum(scoreContribs);

if mag(from, to)==0
    nEdges = nEdges+1;
end
% newScore =-2*tmpSll+log(nSamples)*(nVars+nEdges);
newScore = -2*tmpSll+ bicPenalty(nSamples, nVars, nEdges);


end

function [newScore, scoreContribs, newNComps, newSizes, newComps, newInComponent, newMag,beta,omega,district] = addBidirectedEdge(from, to, nEdges, nComps, sizes, comps, inComponent, mag, covMat, isParent, scores, nSamples, nVars,  tol)

% if edge was directed, update isParent
if mag(to, from)==3
    isParent(from, to)=false;
elseif mag(from, to)==3
    isParent(to, from)=false;
end
[newNComps, newSizes, newComps, newInComponent, k, m] = updateConcomp(from, to, nComps, sizes, comps, inComponent);
% add edge
newMag = mag;
newMag(to, from)=2;newMag(from, to)=2;
component= newComps{k};
[compMag, district] = componentMag(component, nVars, newMag, isParent);
% keep old scores
keepScores =[1:m-1, m+1:nComps];
if k<m
    newScores = scores(keepScores);
else
    newScores= scores;
end
[newScores(k),beta, omega]   = score_contrib(compMag, component, district, newSizes(k), covMat, nSamples,  tol);
scoreContribs = newScores;

tmpSll = (-nSamples/2)*sum(newScores);
tmp_nEdges = nEdges+1;
%newScore =-2*tmpSll+log(nSamples)*(nVars+tmp_nEdges);
newScore = -2*tmpSll+ bicPenalty(nSamples, nVars, tmp_nEdges);



end


function [newScore, scoreContribs, nComps, sizes, comps, inComponent, newMag,beta,omega,district] = removeDirectedEdge(from, to, nEdges, nComps, sizes, comps, inComponent, mag, covMat, isParent, scores, nSamples, nVars,  tol)
newMag = mag;
newMag(from, to)=0;newMag(to, from)=0;
% update district
tmpIsParent = isParent;
tmpIsParent(from, to)=false;
component = comps{inComponent(to)};
[compMag, district] = componentMag(component, nVars, newMag, tmpIsParent);
% get new score
scoreContribs = scores;
[scoreContribs(inComponent(to)),beta,omega] = score_contrib(compMag, component, district, sizes(inComponent(to)), covMat, nSamples, tol);

tmpSll = (-nSamples/2)*sum(scoreContribs);
tmp_nEdges = nEdges+1;
%newScore =-2*tmpSll+log(nSamples)*(nVars+tmp_nEdges);
newScore = -2*tmpSll+ bicPenalty(nSamples, nVars, tmp_nEdges);
end


function [newScore, scoreContribs, newNComps, newNsizes, newComps, newInComponent, newMag] = removeBidirectedEdge(from, to, nEdges, nComps, sizes, comps, inComponent, mag, covMat, isParent, scores, nSamples, nVars,  tol)
% remove edge
newMag = mag;
newMag(to, from)=0;newMag(from, to)=0;
% update components
[newNComps, newNsizes, newComps, newInComponent, k, m] = updateConcompRem(from, to, nComps, sizes, comps, inComponent, mag);
% keep old scores
if isnan(m) % if you have not split the component.
    tmp_scores=scores;
else
    keepScores =[1:k-1, k+1:m-1, m+1:newNComps];
    tmp_scores = nan(1, newNComps);
    tmp_scores(keepScores) = scores([1:k-1,k+1:nComps]);
    % update m-score
    component= newComps{m};
    [compMag, district] = componentMag(component, nVars, newMag, isParent);
    [tmp_scores(m),beta, omega]   = score_contrib(compMag, component, district, newNsizes(m), covMat, nSamples,  tol);
end
%update k-score
component = newComps{k};
[compMag, district] = componentMag(component, nVars, newMag, isParent);
[tmp_scores(k),beta, omega]   = score_contrib(compMag, component, district, newNsizes(k), covMat, nSamples,  tol);

scoreContribs = tmp_scores;

% calculate new score
tmpSll = (-nSamples/2)*sum(tmp_scores);
tmp_nEdges = nEdges-1;

newScore = -2*tmpSll+ bicPenalty(nSamples, nVars, tmp_nEdges);

end

function [isParent, isAncestor, stepAncestor] = updateAncestors(from, to, mag, isParent, isAncestor, stepAncestor, flag)


if flag
    isParent(to, from)=false;
    isParent(from, to)= true;
    [isAncestor, stepAncestor]= warshall2(isParent);
else
    isParent(from, to)= true;
    Anc_from = isAncestor(:, from);
    Desc_to = isAncestor(to, :);
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Anc_from = logical(Anc_from);
    Desc_to = logical(Desc_to);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    stepAncestor(Anc_from, to) = 2;
    stepAncestor(Anc_from, Desc_to) = 2;
    stepAncestor(from, Desc_to)= 2;
    
    Anc_from(from)= true;
    Desc_to(to)= true;
    isAncestor(Anc_from, Desc_to)= true;
end

end

function bp = bicPenalty(nSamples, nVars, nEdges)

bp = log(nSamples)*(2*nVars+nEdges);
end


function bool = isTabu(from, to, iAct, mag, tabu)

bool=false;
switch iAct
    case 1
        mag(from, to)=2; mag(to, from)=3;
    case 2
        mag(from, to)=3; mag(to, from)=2;
    case 3
        mag(from, to)=2; mag(to, from)=2;
    case 4
        mag(from, to)=0; mag(to, from)=0;
end



for m=1:size(tabu,1)
    if all(all(mag==tabu{m}))==1
        bool = true;
        break;
    end
end
end


