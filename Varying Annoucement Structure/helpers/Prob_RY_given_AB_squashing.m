function [p_ab,rho_RY_AB] = Prob_RY_given_AB_squashing(a,b,rhoXY)
%calculate the Prob(RY|AB)
% R:keybit
% X:Alice's measurement
% Y:Bob's measurement
% A:Alice's announcement
% B:Bob's announcement
dimX = 4;
dimY = 5;
dimY_prime = size(rhoXY,2);
dimA = 4;
dimB = 2;
dimR = 3;
% isometries to calculate the distributions of Prob(RBAXY)
% the key, Bob's annoucements and Alice's annoucements are based on X and Y
% registers. These are conditional probabilities 
NA = Produce_NA_squashing();  %generate p(A|XY)
NB = Produce_NB();  %generate p(B|AXY) 
G = Produce_G();  %generate p(R|ABXY)
rho = reshape(rhoXY',[dimX*dimY_prime,1]);

rho_RBAXY = G*NB*NA*rho; % get the joint prob vec Prob(RBAXY)

%obtain the marginals we want
% Prob(RBAY) and % Prob(BA)
rho_RBAY = kron(kron(eye(dimR*dimB*dimA),[1;1;1;1]'),eye(dimY))*rho_RBAXY; 
rho_BA = kron(kron([1;1;1]',eye(dimB*dimA)),[1;1;1;1;1]')*rho_RBAY;

% compute Prob(RY|AB)
p_ab = kron(zket(dimB,b),zket(dimA,a))'*rho_BA;
rho_RY_AB_vec = kron(eye(dimR),kron(kron(zket(dimB,b),zket(dimA,a))',eye(dimY)))*rho_RBAY/p_ab;
rho_RY_AB = reshape(rho_RY_AB_vec,[dimY,dimR])';

end