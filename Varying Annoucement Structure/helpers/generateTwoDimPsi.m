function [psi, rho, rhoA, POVMsA, sqrtPOVMsA] = generateTwoDimPsi(nCut, lambda, sigProb, opt)
%Takes photon number cutoff, photon source intensity and signal
%probabilities. Generates psi in tensor product of Alice, shield and
%channel input spaces. Channel input space is structured as vacuum, H, V,
%H^2, V^2, HV, ... where HV is an excitation in both H and V. This function
%makes a call to Michael's Generate States function. 
arguments
    nCut (1,1) double {mustBeInteger(nCut), mustBeGreaterThan(nCut, 0)}
    lambda (1,1) double {mustBeGreaterThanOrEqual(lambda, 0)}
    sigProb (1, 4) double {mustBeProbDist(sigProb)}
    opt.dimFlag (1,1) double {mustBeOneOrDimA(opt.dimFlag)} = 4;
    opt.catchHermitian (1,1) double {mustBeNumericOrLogical} = 1;
end
pois = @(lamb, n) exp(-lamb)*lamb^n/factorial(n); 
triang = @(n) nchoosek(n+1, 2); 
dimA = 4;
dimFlag = opt.dimFlag; 
dimSig = triang(nCut); %will concatonate flags after
A = 0; %shield A and A prime
POVMsA = cell(1,4); 
POVMsB = cell(1, 2*dimA/2 + 1); 
for a = 1:dimA
    gammaA = 0; 
    gammaB = 0; 
    Aa = 0; 
    for n = 1:nCut-1
        %generate all dimSig dimensional states associated with a
        [H, V, P, M] = generateStates(n); 
        nStates = {H, V, P, M}; % m photon version of these states
        Tn = triang(n);  %"subspaces before" first term of n photon subspace
        probn = pois(lambda, n); 
        %probn = dirac(nDirac, n);
        anState = [zeros(Tn, 1); nStates{a}; zeros(dimSig - triang(n+1),1)]; %position in dimSig dimensional state
        anState = [anState; zeros(dimFlag, 1)]; %flag states 
        ketna = kron( zket(nCut + 1, n + 1), zket(dimA, a) );  % |n>|a>
        Aa = Aa + sqrt(probn)*(ketna*anState');        
        gammaA = gammaA + ketna*ketna'; 
        % gammaB = gammaB + anState*anState'; 
    end
    %vacuum contribution. Can prob be moved to end
    vac = kron(zket((nCut+1), 1), zket(dimA, a)); 
    gammaA = gammaA + vac*vac';
    vac = (vac*zket(dimSig + dimFlag, 1)'); 
    Aa = Aa + sqrt(pois(lambda, 0)) * vac;
    % n > nCut flag state contribution
    aFlagState = kron(zket(nCut + 1, nCut + 1), zket(dimA, a)); % n>nCut flag state for a
    gammaA = gammaA + aFlagState*aFlagState';
    if dimFlag == 4
        fa = a;
    else
        fa = 1;
    end
    aFlagState = (aFlagState * zket(dimSig + dimFlag, dimSig + fa)'); 
    pleqn = sum(arrayfun(@(m) pois(lambda, m), 0:nCut-1)); 
    Aa = Aa + sqrt(1-pleqn)*aFlagState; 

    %add psia to psitot
    A = A + sqrt(sigProb(a))*Aa; 
    POVMsA{a} = gammaA;
    POVMsB{a} = gammaB;
end
%just make sure nothing has gone horribly wrong
%return psi
[U, S, V] = svd(A); 
PiA = 0;
psi = 0;
dimNew = size(V, 1); 
%projOps = arrayfun(@(i) zket(dimNew, i) * U(:, i)', 1:dimNew, UniformOutput = false); 
%veci = arrayfun(@(i) S(i, i) * kron(zket(dimNew, i), V(:, i)), uniformOutput = false);  
for i = 1:dimNew
    PiA = PiA + zket(dimNew, i) * U(:, i)';
    psi =  psi + S(i, i) * kron(zket(dimNew, i), V(:, i));
end
try
    isUnitVector(psi); 
catch error
    rethrow(error); 
end
sqrtPOVMsA = POVMsA; 

for a = 1:numel(POVMsA)
    gammaA = POVMsA{a};
    tol = 1e-10;
    gammaA = PiA * gammaA * PiA' ; %map to new Schmidt basis
    gammaA((gammaA > -tol) & (gammaA < tol)) = 0;

    if opt.catchHermitian == 1
        try
            isHermPSD(gammaA); 
        catch ME
            if strcmp(ME.identifier, 'generateTwoDimPsi:NonHermPSD')
                warning("Correcting Non Hermitian/PSD POVM Element %d of %d.",a,numel(POVMsA));
                gammaA = (gammaA + gammaA')/2;
            end
        end
    end

    [P, D] = eig(gammaA);  %gammaA hermitian psd, find diagonalization
    sqrtEigenvals = real(sqrt(D)); %all eigenvals real, so remove numerical errors
    sqrtEigenvals(sqrtEigenvals < 0) = 0;  %''
    
    %P((P > -tol) & (P < tol)) = 0;

    sqrtGammaA = P*sqrtEigenvals*P'; 
    
    POVMsA{a} = gammaA; 
    sqrtPOVMsA{a} = sqrtGammaA; 
end
rho = psi*psi'; 
rhoA = PartialTrace(rho, 2, [dimNew, dimNew]); 
try
    assert(issymmetric(rhoA));
catch
    warning('Computed rhoA is assymetric. Cannot account for imaginary parts.');
end
rhoA((rhoA > -eps) & (rhoA < eps)) = 0;
end

function[H,V,plus,minus] = generateStates(n)
%pn: the number of photon per signal
% This function will automatically generate 4 states used in lossy SARG04
% protocol. Written by Michael. Added for convenience. 
% dim = pn+2;
arguments
    n (1,1) double {mustBeInteger(n)}
end
dim = n+1;
H = zket(dim,1);
V = zket(dim,2);
plus = 0;
minus = 0;
for i = 1:n-1
    % plus = plus + nchoosek(pn,i)*sqrt(factorial(pn-i)*factorial(i))*zket(dim,i+3);
    % minus = minus + ((-1)^i)*nchoosek(pn,i)*sqrt(factorial(pn-i)*factorial(i))*zket(dim,i+3);
    plus = plus + nchoosek(n,i)*sqrt(factorial(n-i)*factorial(i))*zket(dim,i+2);
    minus = minus + ((-1)^i)*nchoosek(n,i)*sqrt(factorial(n-i)*factorial(i))*zket(dim,i+2);
end
%normalize the states
plus = plus + sqrt(factorial(n))*(zket(dim,1)+zket(dim,2));
minus = minus + sqrt(factorial(n))*(zket(dim,1)+(-1)^(n)*zket(dim,2));
plus = plus/sqrt(plus'*plus);
minus = minus/sqrt(minus'*minus);
end

function isUnitVector(vec)
if ~ismembertol(vec'*vec, 1)
    throw(MException('generateTwoDimPsi:NotUnitVector', ...
        'Computed vector is not a unit vector.'));
end
end
function isHermPSD(F)
if ishermitian(F) == 0 ||  IsPSD(F) == 0
    %warning("on"); 
    %warning("generateTwoDimPsi:nonHermPSD",...
    %    "Matlab believes POVM element %d of %d is Non Hermitian or PSD.",a,numel(POVM));

end
end
function mustBeOneOrDimA(x)
if x ~= 1 && x ~= 4
        throw(MException('generateTwoDimPsi:NotValidFlag', ...
        'Flag Hilbert Space must have dimension 1 or 4 (dimA).'));
end
end