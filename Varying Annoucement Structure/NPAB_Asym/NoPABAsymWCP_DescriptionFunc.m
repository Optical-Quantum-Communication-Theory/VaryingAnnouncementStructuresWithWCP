function [newParams,modParser] = NoPABAsymWCP_DescriptionFunc(params, options, debugInfo)
% Description function for NPAB BB84 with WCPs
% Input parameters:
% * pz: The probability that Alice measures in the Z-basis (for this protocol,
%   it's also the probability that Bob measures in the Z-basis aswell). It
%   must be between 0 and 1.
% Output parameters:
% * observablesJoint: The joint observables for Alice and Bob's measurement
%   of the signals.
% * dimA: dimension of Alice's system.
% * dimB: dimension of Bob's system.
% * rhoA: Alice's reduced density matrix for prepare and measure based
%   protocols.
% * POVMA: Alice's set of POVM operators which she measures her state with
%   in the source replacement scheme.
% * POVMB: Bob's set of POVM operators which he measures his state with.
% * announcementsA: Alice's announcements for each of her POVM operators.
%   Can be integers or strings.
% * announcementsB: Bob's announcements for each of his POVM operators.
%   Can be integers or strings.
% * keyMap: An array of KeyMapElement objects that contain pairs of accepted
%   announcements and an array dictating the mapping of Alice's measurement
%   outcome to key bits (May be written with Strings).
% * krausOps: A cell array of matrices. The Kraus operators that form the G
%   map on Alice and Bob's joint system. These should form a completely
%   postive trace non-increasing linear map. Each Kraus operator must be
%   the same size.
% * keyProj:  A cell array of projection operators that extract the key from
%   G(\rho). These projection operators should sum to identity.
% Options:
% * none
% DebugInfo:
% * krausSum: sum_i K^\dagger_i*K_i which should be <= identity for
%   a CPTNI map.
%
% See also QKDDescriptionModule, BasicBB84_4DAliceKeyRateFunc, makeGlobalOptionsParser
arguments
    params (1,1) struct
    options (1,1) struct
    debugInfo (1,1) DebugInfo
end

%% options parser
optionsParser = makeGlobalOptionsParser(mfilename);
optionsParser.parse(options);
options = optionsParser.Results;


%% module parser
modParser = moduleParser(mfilename);
modParser.addRequiredParam("pz",@(x) mustBeInRange(x,0,1));
% modParser.addRequiredParam("prob0");
modParser.addRequiredParam("lamb");
modParser.addRequiredParam("ncutoff");

%modParser.addRequiredParam("epsSound", @(e) mustBeNonnegative(e));
%modParser.addRequiredParam("epsATfrac", @(f) mustBeInRange(f, 0, 1));

modParser.parse(params);
params = modParser.Results;


%% simple setup
newParams = struct();

%% after doing schmidt decomposition
pz = params.pz;

ketP = [1;1]/sqrt(2);
ketM = [1;-1]/sqrt(2);


nCut = params.ncutoff;
lambda = params.lamb; 
dimB = 3;

[~, ~, rhoAsA, POVMsA, sqrtPOVMsA] = generateTwoDimPsi(nCut, lambda, [pz/2, pz/2, (1-pz)/2, (1-pz)/2], dimFlag = 4);
dimA = size(rhoAsA, 1); 
newParams.rhoA = rhoAsA;
debugInfo.storeInfo("rhoA",rhoAsA);
newParams.dimA = dimA;
newParams.dimB = dimB;

POVMsB = {pz*diag([1,0,0]),pz*diag([0,1,0]),(1-pz)*([ketP;0]*[ketP;0]'),(1-pz)*([ketM;0]*[ketM;0]'), diag([0,0,1])};

newParams.POVMA = POVMsA;
newParams.POVMB = POVMsB;

% Set up a cell array to contain all of the joint observables from
% Alice and Bob's measurments.
observablesJoint = cell(numel(POVMsA),numel(POVMsB));

for indexA = 1:numel(POVMsA)
    for indexB = 1:numel(POVMsB)
        observablesJoint{indexA,indexB} = kron(POVMsA{indexA},POVMsB{indexB});
    end
end

newParams.observablesJoint = observablesJoint;
debugInfo.storeInfo("observablesJoint",observablesJoint);


%% Kraus Ops (for G map)
dimC = 2;

%% Kraus Ops (for G map)
%since Alice is 4D, sqrt(GammaA) = GammaA

krausOpN = kron(zket(2,1), kron(zket(dimC,1), sqrtPOVMsA{1})) + kron(zket(2,2), kron(zket(dimC,1), sqrtPOVMsA{2}));
krausOpN = krausOpN + kron(zket(2,1), kron(zket(dimC, 2), sqrtPOVMsA{3})) + kron(zket(2,2), kron(zket(dimC, 2), sqrtPOVMsA{4}));
krausOpN = kron(krausOpN, [diag([1,1]),[0;0]]); 
krausOps = {krausOpN};
%

newParams.announcementsA = ["N", "N", "N", "N"];
newParams.announcementsB = [newParams.announcementsA, "VAC"]; 
newParams.keyMap = [KeyMapElement("N","N",[1,2,1,2])];


 %Here we compute sum_i K^\dagger_i*K_i. Which should satisfy sum_i
%K^\dagger_i*K_i <= I. A.K.A. the Kraus operators represent a
%completely positive, trace non-increasing linear map.
krausSum = 0;
for index = 1:numel(krausOps)
    krausSum = krausSum+krausOps{index}'*krausOps{index};
end
debugInfo.storeInfo("krausSum",krausSum);

%% key projection
proj0 = kron(diag([1,0]),eye(dimA*(dimB-1)*dimC)); %#
 proj1 = kron(diag([0,1]),eye(dimA*(dimB-1)*dimC)); %#

keyProj = {proj0,proj1};


%% set key map, kraus ops, and announcements in new parameters
newParams.krausOps = krausOps;
newParams.keyProj = keyProj;
end
