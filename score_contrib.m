function [sc,beta,omega] = score_contrib(compMag, component, district, compSize, covMat, nSamples,  tol)
if sum(district)==1
    sc = logdet(2*pi*covMat(district, district))+(nSamples-1)/nSamples;
    beta=0;
    omega=0;
 
else
    curCovMat = covMat(district, district);
    compMag = compMag(district, district);
    [beta, omega, curHatCovMat, ~] = RICF_fit(compMag, curCovMat, tol);
  
    
   
    
%    Rcpp::Named("SigmaHat", graph.getSigma()),
%                               Rcpp::Named("OmegaHat", graph.getOmegaInit()),
%                               Rcpp::Named("BHat", graph.getBInit()),
%                               Rcpp::Named("Iter", counter),
%                               Rcpp::Named("Converged", (convCrit < tol)));
%[beta, omega, hatCovMat, ricf] = RICF_fit(smm, covMat, tol)

    remParents=district;
    remParents(component)=false;
    parInds = remParents(district);

    if any(remParents)
        l1 = compSize*log(2*pi); % part 1
        l2 = logdet(curHatCovMat) - log(prod(diag(curHatCovMat(parInds, parInds)))); % part 2
        l3 = (nSamples-1)/nSamples*(trace(curHatCovMat\curCovMat)-sum(remParents)); % part 3
        sc = l1+l2+l3;
    else
        l1 = compSize*log(2*pi);
        l2 = logdet(curHatCovMat);
        l3 = (nSamples-1)/nSamples*trace(curHatCovMat\curCovMat);
        sc = l1+l2+l3;
    end
    


end
