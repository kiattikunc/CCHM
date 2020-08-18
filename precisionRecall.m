function [precision, recall,bsf,tp,tn,fp,fn,a,i] = precisionRecall(pag, gtpag)
% [PRECISION, RECALL] = PRECISIONRECALL(PCG, GTPCG) Returns precision and
% recall for comparing a  pag to the ground truth pag(as described in Tillman et al, learning
% equivalence classes...)
% precision = # of edges in pag that are in the gtpag (with correct orientations)/# of edges in pag
% recall = # of edges in pag that are in the gtpag (with correct orientation)/# of edges in gtpag
% NOTE: if i understood correctly, they do not take under consideration
% whether the edge is dashed or solid.

nEdgesPag = nnz(pag);
nEdgesGtPag = nnz(gtpag);
nCorrectEdges = nnz(((pag-gtpag)==0)& ((pag'-gtpag')==0) & ~~pag);
%nCorrectEdges = nnz(~triu(pcg-gtpcg & pcg'-gtpcg')&~~pcg);
precision = nCorrectEdges/nEdgesPag;
recall = nCorrectEdges/nEdgesGtPag;
nvar = size(gtpag,2);

tps = nnz(((pag-gtpag)==0)& ((pag'-gtpag')==0) & ~~pag);

fps = nnz(((pag ~= gtpag)|(pag'~=gtpag')) & ~~pag);

fns = nnz(((pag-gtpag)~=0) | ((pag'-gtpag')~=0) & ~~gtpag);

tns = nnz((gtpag==0) & (pag==0));

 tp = 0.5*tps;
 fp =0.5*fps;
 fn = 0.5*fns;
 tn = 0.5*(tns -nvar);
a = 0.5*nEdgesGtPag;
i = (nvar*(nvar-1)/2)-a;

p1=0.5*(tp/a)+(tn/i);
p2=0.5*(fp/i)+(fn/a);

bsf = 0.5*((tp/a)+(tn/i)-(fp/i)-(fn/a));


end

