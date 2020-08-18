function [isAcyclic] = isAcyclic(MAG_graph)

 nvars = size(MAG_graph,1); % KC- Return number of variable

for i = 1:nvars
    for j=1:nvars
        if(MAG_graph(i,j)==2) && (MAG_graph(j,i)==3) 
           MAG_graph(i,j) = 1; 
           MAG_graph(j,i) = 0; 
      
        elseif MAG_graph(i,j)==3 && (MAG_graph(j,i)==2) 
            MAG_graph(i,j) = 0; 
            MAG_graph(j,i) = 1; 
%  

        elseif MAG_graph(i,j)==2 && MAG_graph(j,i)==2
            MAG_graph(i,j) = 0; 
            MAG_graph(j,i) = 0; 
 
        end
        
    
    end
end

   G = digraph(MAG_graph);
%   figure(2)   
%  h = plot(G);
  isAcyclic = isdag(G);
end

                    
                    