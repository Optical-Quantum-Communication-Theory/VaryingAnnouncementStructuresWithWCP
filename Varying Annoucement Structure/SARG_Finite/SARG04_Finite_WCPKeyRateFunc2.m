function [keyRate, modParser] = SARG04_Finite_WCPKeyRateFunc2(params,options,mathSolverFunc,debugInfo)
% Input parameters:
% * dimA: dimension of Alice's system.
% * dimB: dimension of Bob's system.
% * announcementsA: Array of announcements made for each measurement Alice
%   made ordered in the same way as the columns of expectationsJoint.
% * announcementsB: Array of announcements made for each measurement Bob
%   made ordered in the same way as the rows of expectationsJoint.
% * keyMap: An array of KeyMapElement objects that contain pairs of accepted
%   announcements and an array dictating the mapping of Alice's measurement
%   outcome to key bits (May be written with Strings).
% * krausOps: A cell array of matrices. The Kraus operators that form the G
%   map on Alice and Bob's joint system. These should form a completely
%   positive trace non-increasing linear map. Each Kraus operator must be
%   the same size.
% * keyProj:  A cell array of projection operators that extract the key from
%   G(\rho). These projection operators should sum to identity.
% * f: error correction effiency. It 1 means we are correcting at the
%   Shannon limit. A more practical value would be around 1.16.
% * observablesJoint: The joint observables from Alice and Bob's
%   measurements which they perform on the (ideally) max entangled state. The
%   observables must be Hermitian and each must be the size dimA*dimB by
%   dimA*dimB. The observables assume the spaces are ordered A \otimes B.
%   They also should be positive semi-definite and should sum to identity,
%   but this is hard to check because of machine precision issues.
% * expectationsJoint: The joint expectations (as an array) from Alice and
%   Bob's measurements that line up with it's corresponding observable in
%   observablesJoint. These values should be between 0 and 1.
% * rhoA (nan): The density matrix of Alice, if known to be unchanged (for
%   example, in the source replacement scheme). Optional, as entanglement
%   based protocols do not use rhoA, but including rhoA when possible
%   can significantly improve key rate.
% Outputs:
% * keyrate: Key rate of the QKD protocol.
% Options:
% * verboseLevel: (global option) See makeGlobalOptionsParser for details.
% DebugInfo:
% * krausSum: sum_i K^\dagger_i*K_i. For a completely positive trace
%   non-increasing map, this sum should be <=I. 
%
% See also QKDKeyRateModule, BasicBB84_4DAliceDescriptionFunc, makeGlobalOptionsParser
arguments
    params (1,1) struct
    options (1,1) struct
    mathSolverFunc (1,1) function_handle
    debugInfo (1,1) DebugInfo
end

%% options parser
optionsParser = makeGlobalOptionsParser(mfilename);
optionsParser.parse(options);
options = optionsParser.Results;


%% modParser
modParser = moduleParser(mfilename);

modParser.addRequiredParam("observablesJoint",@(x) allCells(x,@ishermitian));
modParser.addRequiredParam("expectationsConditional");
modParser.addRequiredParam("pz");

modParser.addOptionalParam("eta", 1, @(x) mustBeInRange(x, 0, 1));

modParser.addRequiredParam("krausOps", @isCPTNIKrausOps);
modParser.addRequiredParam("keyProj", @(x) mustBeAKeyProj(x));

modParser.addRequiredParam("dimA",@mustBeInteger);
modParser.addRequiredParam("dimB", @mustBeInteger);
modParser.addAdditionalConstraint(@mustBePositive,"dimA")
modParser.addAdditionalConstraint(@mustBePositive,"dimB")
modParser.addAdditionalConstraint(@observablesAndDimensionsMustBeTheSame,["observablesJoint","dimA","dimB"])

modParser.addRequiredParam("f", @(x) mustBeGreaterThanOrEqual(x,1));

modParser.addOptionalParam("rhoA", nan, @isDensityOperator);

modParser.addOptionalParam("blockDimsA", nan); % figure out check on these later (both are given or not, dimensions sum to total dimensions etc.)
modParser.addOptionalParam("blockDimsB", nan);
modParser.addAdditionalConstraint(@(blockDimsA,blockDimsB) ~xor(all(isnan(blockDimsA),"all"),all(isnan(blockDimsB),"all")),["blockDimsA","blockDimsB"]);

%finite size parameters
modParser.addRequiredParam("pGen",@(x) mustBeInRange(x,0,1));
modParser.addRequiredParam("N", @(N) mustBeInteger(N));

modParser.addRequiredParam("epsilon");

modParser.parse(params);

params = modParser.Results;

%% simple setup
debugMathSolver = debugInfo.addLeaves("mathSolver");
mathSolverInput = struct();

%% Turn the expectationCons to squashed expectationjoints.
pz = params.pz;
squashingMap = Produce_Squashing(0);
expectationsJoint = pagemtimes(diag([pz/2,pz/2, (1-pz)/2, (1-pz)/2])*params.expectationsConditional,squashingMap.');

true_expect = [1/2*expectationsJoint(1,:);1/2*expectationsJoint(3,:);...
    1/2*expectationsJoint(2,:);1/2*expectationsJoint(3,:);...
    1/2*expectationsJoint(2,:);1/2*expectationsJoint(4,:);...
    1/2*expectationsJoint(1,:);1/2*expectationsJoint(4,:)];

params.expectationsJoint = true_expect;
debugInfo.storeInfo("true_expect",true_expect);

%% Error correction
% The error correction for SARG04 protocol cannot be automatically handled
% by the software due to the B92 like announcement structure. Here, we use
% the original statistic (click patterns for Alice and Bob rho_XY) to
% obtain probability distributions Prob(R(key)Y(Bob's outcome)|AB(Alice and
% Bob's announcements)). We use this to compute the error correction cost.
rhoXY = diag([pz/2,pz/2,(1-pz)/2,(1-pz)/2])*params.expectationsConditional(:,:,1);
[gains,conExps] = Produce_gains_joints_squashing(rhoXY);
% we need the error part.
gains = gains(5:8);
conExps = conExps(5:8);
deltaLeak = generalECFunc(gains',conExps,params.f);
debugInfo.storeInfo("deltaLeak",deltaLeak);


%% translate for the math solver
mathSolverInput.krausOps = params.krausOps;
mathSolverInput.keyProj = params.keyProj;

if ~isnan(params.rhoA)
    mathSolverInput.rhoA = params.rhoA;
end

numObs = numel(params.observablesJoint);

%% setting up Finite size constrains

pGen = params.pGen; 
epsilon = params.epsilon;
epsAT = epsilon.AT;

%% ADAPTIVE LENGTH NEW CONSTRAINTS %%
constraintOps = params.observablesJoint; 
sigmaAlphabet = numel(constraintOps) + 2; 
N = params.N; 
fObs = true_expect; 

nTest = floor((1 - pGen) * N); 
testStatistics = floor(nTest * fObs); 
testStatistics = testStatistics(:); 

[~, bounds] = binofit(testStatistics, N, epsAT / (2 * sigmaAlphabet)); 

kappaL = params.expectationsJoint - reshape(bounds(:, 1), size(fObs)); 
kappaU = - params.expectationsJoint + reshape(bounds(:, 2), size(fObs));

testConstraintsL = (params.expectationsJoint - kappaL) / (1 - pGen);
testConstraintsU = (params.expectationsJoint + kappaU) / (1 - pGen); 


testRoundConstraints = arrayfun(@(index) InequalityConstraint(...
    constraintOps{index}, testConstraintsL(index), testConstraintsU(index)), 1:numObs);
mathSolverInput.inequalityConstraints = testRoundConstraints; 


%%
% if block diag information was give, then pass it to the solver.
if ~all(isnan(params.blockDimsA),"all")
    mathSolverInput.blockDimsA = params.blockDimsA;
    mathSolverInput.blockDimsB = params.blockDimsB;
end

% now we call the math solver function on the formulated inputs, with
% it's options.
[relEnt,~] = mathSolverFunc(mathSolverInput,debugMathSolver);

keyRate = renyiEntropyAdaptiveKeyRateFunc(relEnt, deltaLeak, gains,sigmaAlphabet,params, debugInfo);

if options.verboseLevel>=1
    %ensure that we cut off at 0 when we display this for the user.
    fprintf("Key rate: %e\n",max(keyRate,0));
end

%set the upper bound as well for debugging
if isfield(debugMathSolver.info,"relEntUpperBound")
    keyRateUpperBound = debugMathSolver.info.relEntUpperBound - deltaLeak;
    debugInfo.storeInfo("keyRateUpperBound",keyRateUpperBound)

    if options.verboseLevel>=2
        fprintf("Key rate upper bound: %e",max(keyRateUpperBound,0))
    end
end
end
function [keyRate, debugInfo] = renyiEntropyAdaptiveKeyRateFunc(relEnt, deltaLeak, gains, sigmaAlphabet, params, debugInfo)
    %Compute Finite Size Key Rate
    
    dimA = 2; % Alice key register dimension
    pGen = params.pGen;  %Probability of Generation Round
    epsPA = params.epsilon.PA;  %Privacy Amplification Security
    epsEC = params.epsilon.EC;
    epsAT = params.epsilon.AT; 
    N = params.N;  %Total signal count
    nSift = floor(sum(gains)* N); %Generation signals surviving sifting
    c = sqrt(log2(1/epsPA)/nSift)/log2(dimA+1);
    alpha = 1 + c; %Renyi entropy order
    
    pSiftGen = sum(gains) * pGen; 

    [~, bounds] = binofit(floor(pSiftGen * N), N, epsAT / (2 * sigmaAlphabet)); 
    prefactor = pSiftGen / bounds(2); 

    ECCost = pGen * deltaLeak + log2(2/epsEC)/N; 
    PATerm = alpha / c * (log2(1/2/epsPA) + 2/alpha)/N; 
    renyiTerm = sqrt(nSift) * c * (log2(dimA+1))^2/N;

    debugInfo.storeInfo("kappaPrefactor", prefactor); 
    debugInfo.storeInfo("ErrorCorrectionCost", ECCost);
    debugInfo.storeInfo("PrivacyAmplificationTerm", PATerm); 
    debugInfo.storeInfo("RenyiEntropyCorrection", renyiTerm); 

    keyRate = prefactor * relEnt - ECCost - PATerm - renyiTerm;
end

function observablesAndDimensionsMustBeTheSame(observables,dimA,dimB)
if ~allCells(observables,@(x) size(x,1) == dimA*dimB)
    throwAsCaller(MException("BasicBB84_4DAliceKeyRateFunc:ObservablesAndDimensionsMustBeTheSame","The Observables must have the same dimensions as Alice and Bob multiplied together."));
end
end
