function NLV = lint_score(x,Gamma) 
    nnlv = size(x,1);                         
    nobs = size(Gamma,1);
    NLV = zeros(nobs,nnlv);                     
    for j = 1:nnlv
        NLV(:,j) = Gamma(:,x(j,1)).*Gamma(:,x(j,2));  
    end