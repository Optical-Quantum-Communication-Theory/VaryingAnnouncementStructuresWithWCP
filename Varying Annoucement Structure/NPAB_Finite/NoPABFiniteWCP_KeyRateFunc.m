function [keyRate, modParser] = NoPABFiniteWCP_KeyRateFunc(params,options,mathSolverFunc,mathSolverOptions,debugInfo)
% Key rate function for NPAB BB84 with WCPs.
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
% * expectationsJoint: The joint expectations (as an array) from Alice and
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
    mathSolverOptions (1,1) struct
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

modParser.addRequiredParam("f", @(x) mustBeGreaterThanOrEqual(x,1));

modParser.addOptionalParam("rhoA", nan, @isDensityOperator);
modParser.addRequiredParam("pz"); 
modParser.addOptionalParam("blockDimsA", nan); 
modParser.addOptionalParam("blockDimsB", nan);
modParser.addAdditionalConstraint(@(blockDimsA,blockDimsB) ~xor(all(isnan(blockDimsA),"all"),all(isnan(blockDimsB),"all")),["blockDimsA","blockDimsB"]);

modParser.addRequiredParam("pGen",@(x) mustBeInRange(x,0,1));
modParser.addRequiredParam("N", @(N) mustBeInteger(N));
modParser.addRequiredParam("epsilon");
modParser.addRequiredParam("t", @(t) t == 0); 

modParser.parse(params);

params = modParser.Results;

%% simple setup
debugMathSolver = debugInfo.addLeaves("mathSolver");

mathSolverInput = struct();

%% Turn the expectationCons to squashed expectationjoints.
pz = params.pz; 
%squashingMap = passiveSquashingMap();
randDC = activeSquashingMap(); %randomly assign double clicks
flagDC = activeSquashingMapDCFlag(); %flag double click events (used for EC)

squashedExpectationsJoint = params.expectationsConditional*transpose(randDC); 
expectationsJoint = diag([pz/2,pz/2, (1-pz)/2, (1-pz)/2])*squashedExpectationsJoint; 

params.expectationsJoint = expectationsJoint;


squashedExpectationsJointDC = params.expectationsConditional*transpose(flagDC); 
expectationsJointDC = diag([pz/2,pz/2, (1-pz)/2, (1-pz)/2])*squashedExpectationsJointDC;

announcementsBobDC = [params.announcementsB, "DC"];
[deltaLeak, gains] = errorCorrectionCost(params.announcementsA,announcementsBobDC,...
    expectationsJointDC,params.keyMap,params.f);

%N = params.N;
%nSift = gains*N; 


%[deltaLeak, gains] = errorCorrectionCost(params.announcementsA,params.announcementsB,...
%    expectationsJoint,params.keyMap,params.f);

debugInfo.storeInfo("deltaLeak",deltaLeak);
debugInfo.storeInfo("expectationsJoint",expectationsJoint);
debugInfo.storeInfo("expectationsJointDC",expectationsJointDC);
%% Error correction

mathSolverInput.krausOps = params.krausOps;
mathSolverInput.keyProj = params.keyProj;
% also include rhoA from the description if it was given
if ~isnan(params.rhoA)
    mathSolverInput.rhoA = params.rhoA;
end

pGen = params.pGen; 
epsilon = params.epsilon;
epsAT = epsilon.AT;

numObs = numel(params.observablesJoint);
%% ADAPTIVE LENGTH CONSTRAINTS %%
constraintOps = params.observablesJoint; 
sigmaAlphabet = numel(constraintOps) + 2; 
N = params.N; 
fObs = params.expectationsJoint; 

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

% now we call the math solver function on the formulated inputs, with
% it's options.
[relEnt,~] = mathSolverFunc(mathSolverInput,mathSolverOptions,debugMathSolver);
%store the key rate (even if negative)
debugInfo.storeInfo("relEnt", relEnt); 
%keyRate = computeFiniteKeyRate(relEnt, deltaLeak, gains, params); 
%keyRate = relEnt - deltaLeak;


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

function mapping = activeSquashingMap()
    mapping = [[0;0;0;0;1],[1;0;0;0;0],[0;1;0;0;0],[1/2;1/2;0;0;0],[0;0;0;0;1],[0;0;1;0;0],[0;0;0;1;0],[0;0;1/2;1/2;0]];
end
function mapping = activeSquashingMapDCFlag()
    mapping = [[0;0;0;0;1;0],[1;0;0;0;0;0],[0;1;0;0;0;0],[0;0;0;0;0;1],[0;0;0;0;1;0],[0;0;1;0;0;0],[0;0;0;1;0;0],[0;0;0;0;0;1]];
end
function mapping = activeSquashingMapDCVac()
    mapping = [[0;0;0;0;1],[1;0;0;0;0],[0;1;0;0;0],[0;0;0;0;1],[0;0;0;0;1],[0;0;1;0;0],[0;0;0;1;0],[0;0;0;0;1]];
end
function mapping = passiveSquashingMap()
% squashes detector patterns from H,V,R,L (as bit string patterns 0000,
% 1000, ..., 1111) to signal+vac values H,V,R,L, (double clicks), CC, vac.

function mapping = quickMap(mapping,pattern,remaping)
    mapping(:,sub2indPlus(2*ones(1,numel(pattern)),pattern+1)) = remaping;
end

mapping = zeros(5,16);

% The vast majority of the squashed bits are cross clicks that are mapped
% to vac for discarding. We will replace patterns that don't represent
% cross clicks in later steps.
%mapping(5,:) = 1;

%order S_( H V D A ) D_( HV RL ) CC_(ANY) NON

% zero clicks to vac
mapping = quickMap(mapping,[0,0,0,0],[0,0,0,0,1]);

% single clicks to single clicks
mapping = quickMap(mapping,[1,0,0,0],[1,0,0,0,0]); % H
mapping = quickMap(mapping,[0,1,0,0],[0,1,0,0,0]); % V
mapping = quickMap(mapping,[0,0,1,0],[0,0,1,0,0]); % D
mapping = quickMap(mapping,[0,0,0,1],[0,0,0,1,0]); % A

% double clicks
mapping = quickMap(mapping,[1,1,0,0],[1/2,1/2,0,0,0]); % Z (HV)
mapping = quickMap(mapping,[0,0,1,1],[0,0,1/2,1/2,0]); % X (RL)

%cross clicks map to vac
% 2 clicks
mapping = quickMap(mapping,[1,0,1,0],[0,0,0,0,1]); %HR
mapping = quickMap(mapping,[1,0,0,1],[0,0,0,0,1]); %HL
mapping = quickMap(mapping,[0,1,1,0],[0,0,0,0,1]); %VR
mapping = quickMap(mapping,[0,1,0,1],[0,0,0,0,1]); %VL
% 3 clicks
mapping = quickMap(mapping,[1,1,1,0],[0,0,0,0,1]); %HV R
mapping = quickMap(mapping,[1,1,0,1],[0,0,0,0,1]); %HV L
mapping = quickMap(mapping,[1,0,1,1],[0,0,0,0,1]); %H RL
mapping = quickMap(mapping,[0,1,1,1],[0,0,0,0,1]); %V RL
%4 clicks
mapping = quickMap(mapping,[1,1,1,1],[0,0,0,0,1]); %HV RL

end
function observablesAndDimensionsMustBeTheSame(observables,dimA,dimB)
if ~allCells(observables,@(x) size(x,1) == dimA*dimB)
    throwAsCaller(MException("BasicBB84_4DAliceKeyRateFunc:ObservablesAndDimensionsMustBeTheSame","The Observables must have the same dimensions as Alice and Bob multiplied together."));
end
end

function mustBeSizedLikeAnnouncements(jointExpectations,announcementsA,announcementsB)
if ~isequal(size(jointExpectations),[numel(announcementsA),numel(announcementsB)])
    throwAsCaller(MException("BasicBB84_4DAliceKeyRateFunc:jointKeyDoesNotHaveSameSizeAsAnnouncements",...
        "The joint key distribution must have size numel(announcementsA) by numel(announcementsB)."))
end
end