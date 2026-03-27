function[H,V,plus,minus,dim] = Generate_States_test(pn)
%pn: the number of photon per signal
% This function will automatically generate 4 states with pn photons used in lossy SARG04
% protocol.
if pn<=3
    dim = pn+1;
    H = zket(dim,1);
    V = zket(dim,2);
% we order the entries so that the first two entries correspond to H and V
% |pn,0>_HV and |0,pn>_HV
    plus = 0;
    minus = 0;
    for i = 1:pn-1
        plus = plus + nchoosek(pn,i)*sqrt(factorial(pn-i)*factorial(i))*zket(dim,i+2);
        minus = minus + ((-1)^i)*nchoosek(pn,i)*sqrt(factorial(pn-i)*factorial(i))*zket(dim,i+2);
    end
    % add |pn,0>_HV and |0,pn>_HV component
    plus = plus + sqrt(factorial(pn))*(zket(dim,1)+zket(dim,2));
    minus = minus + sqrt(factorial(pn))*(zket(dim,1)+(-1)^(pn)*zket(dim,2));
    
    %normalize the states
    plus = plus/sqrt(plus'*plus);
    minus = minus/sqrt(minus'*minus);
   
else
    [H,V,plus,minus]=photons_state(pn);
    dim = 4;
end
end