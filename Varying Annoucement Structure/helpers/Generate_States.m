function[H,V,plus,minus,dim] = Generate_States(pn)
%pn: the number of photon per signal
% This function will automatically generate 4 states used in lossy SARG04
% protocol.
% dim = pn+2;
dim = pn+1;
H = zket(dim,1);
V = zket(dim,2);

plus = 0;
minus = 0;
for i = 1:pn-1
    % plus = plus + nchoosek(pn,i)*sqrt(factorial(pn-i)*factorial(i))*zket(dim,i+3);
    % minus = minus + ((-1)^i)*nchoosek(pn,i)*sqrt(factorial(pn-i)*factorial(i))*zket(dim,i+3);
    plus = plus + nchoosek(pn,i)*sqrt(factorial(pn-i)*factorial(i))*zket(dim,i+2);
    minus = minus + ((-1)^i)*nchoosek(pn,i)*sqrt(factorial(pn-i)*factorial(i))*zket(dim,i+2);
end

%normalize the states
plus = plus + sqrt(factorial(pn))*(zket(dim,1)+zket(dim,2));
minus = minus + sqrt(factorial(pn))*(zket(dim,1)+(-1)^(pn)*zket(dim,2));


plus = plus/sqrt(plus'*plus);
minus = minus/sqrt(minus'*minus);

end