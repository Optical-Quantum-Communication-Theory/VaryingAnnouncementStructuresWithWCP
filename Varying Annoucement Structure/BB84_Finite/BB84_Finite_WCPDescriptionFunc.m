function [newParams,modParser] = BB84_Finite_WCPDescriptionFunc(params, options, debugInfo)
% Input parameters:
% * pz: The probability that Alice measures in the Z-basis (for this protocol,
%   it's also the probability that Bob measures in the Z-basis aswell). It
%   must be between 0 and 1.
% * lamb: The mean photon number for the Weak Coherent Pulse used in the
%   protocol
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
modParser.addRequiredParam("lamb");

modParser.parse(params)
params = modParser.Results;


%% simple setup
newParams = struct();

% To make it easier to write, lets just add a local variable for pz and
% lamb
% with the correct name.
lamb = params.lamb;
pz = params.pz;

ketP = [1;1]/sqrt(2);
ketM = [1;-1]/sqrt(2);
pSig = [pz/2,pz/2,(1-pz)/2,(1-pz)/2];

%POVMsA (Eq.B1)
POVMsA = {};
dim_signal = 7;
for i = 1:4
POVM = 0;
    for n = 1:3
        vec = kron(zket(4,i),zket(3,n));

        POVM = POVM + vec*vec';
    end

POVMsA = [POVMsA,{POVM}];
end

%construct the origial source repalced state and perform Schmidt Decomp
psiAAprime = 0;

for i = 1:4
signal=0;
[H,V,plus,minus]=Generate_States(1);
states_at_1 = {H,V,plus,minus};
states_diag = {zket(4,1),zket(4,2),zket(4,3),zket(4,4)};
states = {zket(dim_signal,1),[0;states_at_1{i};zeros(7-3,1)],[zeros(3,1);states_diag{i}]};
for n=1:numel(states)
     signal = signal+ sqrt(pSig(i))*sqrt(pois(n-1,lamb,numel(states)-1))*(kron(zket(4,i),zket(3,n))*states{n}');
end
psiAAprime = psiAAprime + signal;
end

% singular value decomposition
[U, S, V] = svd(psiAAprime);
% getting new POVMsA
PiA = 0;
dimnew  = size(V,1);
for index = 1:dimnew
    PiA = PiA + zket(dimnew,index)*U(:,index)';
end
for i = 1:4
    POVMsA{i} = PiA*POVMsA{i}*PiA';
end
%% getting new share state
newrhoAAprime = 0;
for index = 1:dimnew
    newrhoAAprime = newrhoAAprime + S(index,index)* kron(zket(dimnew,index),V(:,index));
end
rhoA = PartialTrace(newrhoAAprime*newrhoAAprime',2,[dimnew,dimnew]);

dimA = dimnew;
%POVMsB
POVMsB = {pz*diag([1,0,0]),pz*diag([0,1,0]),(1-pz)*([ketP;0]*[ketP;0]'),(1-pz)*([ketM;0]*[ketM;0]'), diag([0,0,1])};
dimB = 3;

newParams.rhoA =rhoA;
debugInfo.storeInfo("rhoA",rhoA);
debugInfo.storeInfo("POVMsA",POVMsA)
newParams.dimA = dimA; %we want to ouput these as new parameters to merge onto the list of all other parameters
newParams.dimB = dimB;
%%
newParams.POVMA = POVMsA;
newParams.POVMB = POVMsB;
newParams.announcementsA = ["Z","Z","X","X"];
newParams.announcementsB = ["Z","Z","X","X","0"];
newParams.keyMap = [KeyMapElement("Z","Z",[1,2,1,2]), KeyMapElement("X","X",[1,2,1,2])];
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
dimC = 2;

%% Kraus Ops (for G map)
krausOpZ = sqrt(pz)*kron((kron(zket(2,1), sqrtm(POVMsA{1})) + kron(zket(2,2), sqrtm(POVMsA{2}))), kron(diag([1;1;0]), zket(dimC,1))); %#
krausOpX = sqrt(1-pz)*kron((kron(zket(2,1), sqrtm(POVMsA{3})) + kron(zket(2,2), sqrtm(POVMsA{4}))), kron(diag([1;1;0]), zket(dimC,2))); %#


krausOps = {krausOpZ,krausOpX};
krausSum = 0;
for index = 1:numel(krausOps)
    krausSum = krausSum+krausOps{index}'*krausOps{index};
end
debugInfo.storeInfo("krausSum",krausSum);

%% key projection
 proj0 = kron(diag([1,0]),eye(dimA*dimB*dimC)); %#
 proj1 = kron(diag([0,1]),eye(dimA*dimB*dimC)); %#

keyProj = {proj0,proj1};


%% set key map, kraus ops, and announcements in new parameters
newParams.krausOps = krausOps;
newParams.keyProj = keyProj;
newParams.blockDimsA = 7;
newParams.blockDimsB = [2,1];
end

function[prob] = pois(x,lamb,n_cutoff)
if x<n_cutoff
    prob = lamb^x*exp(-lamb)/factorial(x);
else
    sum=0;
    for number = 0:(n_cutoff-1)
        sum = sum + pois(number,lamb,n_cutoff);
    end
    prob = 1-sum;
end
end
