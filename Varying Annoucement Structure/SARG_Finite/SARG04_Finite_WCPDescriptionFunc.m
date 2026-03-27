function [newParams,modParser] = SARG04_Finite_WCPDescriptionFunc(params, options, debugInfo)
% Input parameters:
% * pz: The probability that Alice measures in the Z-basis (for this protocol,
%   it's also the probability that Bob measures in the Z-basis as well). It
%   must be between 0 and 1.
% * lamb: The intensity of weak coherent pulse. It must be positive.
% * ncutoff: The photon number cutoff. It will calculate key rate
%   contributions for photon number less than ncutoff.

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
%   positive trace non-increasing linear map. Each Kraus operator must be
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
%WCP parameters
modParser = moduleParser(mfilename);
modParser.addRequiredParam("pz",@(x) mustBeInRange(x,0,1));
modParser.addRequiredParam("lamb");
modParser.addRequiredParam("ncutoff");

modParser.parse(params);
params = modParser.Results;


%% simple setup
newParams = struct();
lamb = params.lamb;
pz = params.pz;

ketP = [1;1]/sqrt(2);
ketM = [1;-1]/sqrt(2);

%POVMsA,POVMsB and rhoA
[rhoA,POVMsA,dimA]=Generate_prepared_state_test(lamb,params.ncutoff,pz);

POVMsB = {pz*diag([1,0,0]),pz*diag([0,1,0]),(1-pz)*([ketP;0]*[ketP;0]'),(1-pz)*([ketM;0]*[ketM;0]'), diag([0,0,1])};
dimB = 3;
%we want to ouput these as new parameters to merge onto the list of all other parameters
newParams.rhoA =rhoA;
debugInfo.storeInfo("rhoA",rhoA);
newParams.dimA = dimA;
newParams.dimB = dimB;

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

%% Kraus Ops
krausOps = {Improved_Kraus2(1,POVMsA),Improved_Kraus2(2,POVMsA),Improved_Kraus2(3,POVMsA),Improved_Kraus2(4,POVMsA)};
%Here we compute sum_i K^\dagger_i*K_i. Which should satisfy sum_i
%K^\dagger_i*K_i <= I. A.K.A. the Kraus operators represent a
%completely positive, trace non-increasing linear map.
krausSum = 0;
for index = 1:numel(krausOps)
    krausSum = krausSum+krausOps{index}'*krausOps{index};
end
debugInfo.storeInfo("krausSum",krausSum);

%% key projection
proj0 = kron(diag([1,0]),eye(12*dimA)); %#
proj1 = kron(diag([0,1]),eye(12*dimA)); %#
keyProj = {proj0,proj1};


%% set key map, kraus ops, and announcements in new parameters
newParams.krausOps = krausOps;
newParams.keyProj = keyProj;
newParams.blockDimsA = [1,2,3,4,4];
newParams.blockDimsB = [2,1];


end
%% A function for choosing the correct POVMs 

function[kraus] = Improved_Kraus2(a,POVMsA)
kraus = 0;
pz=1/2;
ketP = [1;1]/sqrt(2);
ketM = [1;-1]/sqrt(2);
POVMsB = {pz*diag([1,0,0]),(1-pz)*([ketP;0]*[ketP;0]'),pz*diag([0,1,0]),(1-pz)*([ketM;0]*[ketM;0]')}; %same one as defined above.
POVMsB_sum = {POVMsB{3}+POVMsB{4},POVMsB{1}+POVMsB{4},POVMsB{1}+POVMsB{2},POVMsB{2}+POVMsB{3}};
for r = 1:2
    % (eq.D1)
    kraus = kraus + kron(kron(kron(zket(2,r),sqrtm(POVMsA{2*(a-1)+r})),sqrtm(POVMsB_sum{a})),zket(4,a));
end
end
