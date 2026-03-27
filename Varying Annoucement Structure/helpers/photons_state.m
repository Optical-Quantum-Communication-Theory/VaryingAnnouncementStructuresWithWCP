function [newH,newV,newP,newM]=photons_state(n)
% we use the same 4 state with n excitations. But we simplify it by finding
% a new basis using gram schimdt procedure such that the entries in rest of the dimensions is 0.
[H,V,P,M] = Generate_States(n);

orth3 = P - (H'*P)*H - (V'*P)*V;
orth3 = orth3/sqrt(orth3'*orth3);
orth4 = M - H*M'*H - V*M'*V - orth3'*M*orth3;
orth4 = orth4/sqrt(orth4'*orth4);

newP = [P(1);P(2);P'*orth3;0];
newM = [M(1);M(2);M'*orth3;M'*orth4];

newH = zket(4,1);
newV = zket(4,2);
end