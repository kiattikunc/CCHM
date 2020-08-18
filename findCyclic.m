function [cyclic,n_cycle] = findCyclic(MAG_graph)


 nvars = size(MAG_graph,1); % KC- Return number of variable
A=zeros(nvars,nvars);

for i = 1:nvars
    for j=1:nvars
        if(MAG_graph(i,j)==2) && (MAG_graph(j,i)==3) 
          A(i,j) = 1; 
          A(j,i) = 0; 
     
%  

        elseif MAG_graph(i,j)==2 && MAG_graph(j,i)==2
           A(i,j)= 0; 
           A(i,j) = 0; 
 
        end
        
    
    end
end
% A = [0 1 0 1 0 0 ; 
%      0 0 1 0 0 0 ; 
%      0 0 0 1 0 1 ;
%      0 0 0 0 1 1 ; 
%      1 0 0 0 0 0 ; 
%      0 0 0 0 0 0] ; 
   
G = digraph(A);
G.Nodes.Names=cumsum(ones(1,nvars))';
%G.Edges.Weight = [1 1 1 1 1 1 1 1]';
% figure(2) 
%  p = plot(G);
% % 
 e = dfsearch(G, 1, 'edgetodiscovered', 'Restart', true);
 n = size(e,1);
 n_cycle =n;
 dummy = numnodes(G) + 1;
% 
%  fprintf('n=%d \n',n)
%   fprintf('dummy=%d \n',dummy);
%  fprintf('e(:,2)=%d',e(1,2));
  G = addedge(G, repmat(dummy, n, 1) , e(1,2),[0]);
  cyclic = shortestpathtree(G, dummy , e(1,1), 'OutputForm', 'cell', 'Method', 'unweighted');
 

 end
