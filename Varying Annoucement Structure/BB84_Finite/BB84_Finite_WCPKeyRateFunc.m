function [keyRate, modParser] = BB84_Finite_WCPKeyRateFunc(params,options,mathSolverFunc,debugInfo)
% BasicBB84_4DAliceKeyRateFunc A simple key rate function for a qubit BB84 
% protocol with no loss.
%
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
%   postive trace non-increasing linear map. Each Kraus operator must be
%   the same size.
% * keyProj:  A cell array of projection operators that extract the key from
%   G(\rho). These projection operators should sum to identity.
% * f: error correction effiency. It 1 means we are correcting at the
%   Shannon limit. A more practicle value would be around 1.16.
% * observablesJoint: The joint observables from Alice and Bob's
%   measurments which they perform on the (idealy) max entangled state. The
%   observables must be hermitian and each must be the size dimA*dimB by
%   dimA*dimB. The observables assume the spaces are ordered A \otimes B.
%   They also should be positive semi-definite and should sum to identity,
%   but this is hard to check because of machine precision issues.
% * expectationsConditional: The expectations (as an array) from Alice and
%   Bob's measurements that line up with it's corresponding observable in
%   observablesJoint. These values should be betwen 0 and 1.
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
modParser.addOptionalParam("eta", 1, @(x) mustBeInRange(x, 0, 1));
modParser.addRequiredParam("pz");


% modParser.addRequiredParam("rhoXY");
% modParser.addAdditionalConstraint(@isEqualSize,["observablesJoint","expectationsJoint"]);

modParser.addRequiredParam("krausOps", @isCPTNIKrausOps);

modParser.addRequiredParam("keyProj", @(x) mustBeAKeyProj(x));

modParser.addRequiredParam("dimA",@mustBeInteger);
modParser.addRequiredParam("dimB", @mustBeInteger);
modParser.addAdditionalConstraint(@mustBePositive,"dimA")
modParser.addAdditionalConstraint(@mustBePositive,"dimB")
modParser.addAdditionalConstraint(@observablesAndDimensionsMustBeTheSame,["observablesJoint","dimA","dimB"])

modParser.addRequiredParam("announcementsA")
modParser.addRequiredParam("announcementsB")
modParser.addRequiredParam("keyMap",@(x)mustBeA(x,"KeyMapElement"))
% modParser.addAdditionalConstraint(@mustBeSizedLikeAnnouncements,["expectationsJoint","announcementsA","announcementsB"])

modParser.addRequiredParam("f", @(x) mustBeGreaterThanOrEqual(x,1));

modParser.addOptionalParam("rhoA", nan, @isDensityOperator);

modParser.addOptionalParam("blockDimsA", nan); % figure out check on these later (both are given or not, dimensions sum to total dimensions etc.)
modParser.addOptionalParam("blockDimsB", nan);
modParser.addAdditionalConstraint(@(blockDimsA,blockDimsB) ~xor(all(isnan(blockDimsA),"all"),all(isnan(blockDimsB),"all")),["blockDimsA","blockDimsB"]);

%finite size parameters

modParser.addRequiredParam("pGen",@(x) mustBeInRange(x,0,1));
modParser.addRequiredParam("N", @(N) mustBeInteger(N));
modParser.addRequiredParam("epsilon");
% modParser.addRequiredParam("t", @(t) t == 0); 


modParser.parse(params);

params = modParser.Results;

%% simple setup
debugMathSolver = debugInfo.addLeaves("mathSolver");

mathSolverInput = struct();
pz=params.pz;
%% Turn the expectationCons to squashed expectationjoints.
squashingMap = Produce_Squashing(0);
expectationsJoint = pagemtimes(diag([pz/2,pz/2, (1-pz)/2, (1-pz)/2])*params.expectationsConditional,squashingMap.');
%% Error correction
params.expectationsJoint = expectationsJoint;
[deltaLeak,gains] = errorCorrectionCost(params.announcementsA,params.announcementsB,...
    params.expectationsJoint,params.keyMap,params.f);
debugInfo.storeInfo("deltaLeak",deltaLeak);
debugInfo.storeInfo("expectationsJoint",expectationsJoint);

%% translate for the math solver
mathSolverInput.krausOps = params.krausOps;
mathSolverInput.keyProj = params.keyProj;
% also include rhoA from the description if it was given
if ~isnan(params.rhoA)
    mathSolverInput.rhoA = params.rhoA;
end

numObs = numel(params.observablesJoint);

pGen = params.pGen; 
epsAT = params.epsilon.AT;
%% ADAPTIVE LENGTH NEW CONSTRAINTS %%
constraintOps = params.observablesJoint; 
sigmaAlphabet = numel(constraintOps) + 2; 
N = params.N; 
fObs = params.expectationsJoint; 

%Eq A17 bounds
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

% if block diag information was give, then pass it to the solver.
if ~all(isnan(params.blockDimsA),"all")
    mathSolverInput.blockDimsA = params.blockDimsA;
    mathSolverInput.blockDimsB = params.blockDimsB;
end

[relEnt,~] = mathSolverFunc(mathSolverInput,debugMathSolver);
%store the key rate (even if negative)
keyRate = renyiEntropyAdaptiveKeyRateFunc(relEnt, deltaLeak, gains, sigmaAlphabet, params, debugInfo); 

if options.verboseLevel>=1
    %ensure that we cut off at 0 when we display this for the user.
    fprintf("Key rate: %e\n",max(keyRate,0));
end

%set the upper bound as well for debuging
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
