function [newParams, modParser]= SARG04WCP2ChannelFunc(params,options,debugInfo)
% SARGWCP2ChannelFunc_depol gives the clicking patterns given the
% parameteres related to experiments.
% Input parameters:
% * eta: the transmissivity of the quantum channel; equivalent to 1 - loss.
%   Must be between 0 and 1 inclusive.
% * detectorEfficiency: the efficiency of Bob's detectors. Must be between
%   0 and 1 inclusive
% * lamb: The intensity of weak coherent pulse. It must be positve
% * misalignmentAngle: Angle Alice and Bob's bases are misaligned by around
%   the Y-axis. For example, Bob's detectors could be slightly rotated away
%   from the incoming signals. Angles must be real numbers. The angle is
%   the physical rotation angle and not th rotation angle on the bloch
%   sphere.
% * visibilityAngle: the depolarization factor. If visibilityAngle = 0,
%   then there is no depolarization.
% Output parameters:
% * expectationsConditional: The conditional expectations (as a 3D array) from 
%   Alice and Bob's measurements. This should be organized as a 4 x 16 x n 
%   array, where n = the number of intensities used in the decoy protocol.
%   In each table, these should line up with the corresponding observables 
%   at each entry.
% Options:
% * None.
% DebugInfo:
% * None.
%
% See also QKDChannelModule,  makeGlobalOptionsParser
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
modParser.addRequiredParam("pz", @(x) mustBeInRange(x, 0, 1));
modParser.addAdditionalConstraint(@isscalar,"pz");

% modParser.addRequiredParam("prob0")
% modParser.addAdditionalConstraint(@isscalar,"prob0");
modParser.addRequiredParam("lamb");

modParser.addOptionalParam("eta", 1, @(x) mustBeInRange(x, 0, 1));
modParser.addAdditionalConstraint(@isscalar,"eta");

modParser.addOptionalParam("detectorEfficiency", 1, @(x) mustBeInRange(x, 0, 1));
modParser.addAdditionalConstraint(@isscalar,"detectorEfficiency");

modParser.addOptionalParam("misalignmentAngle",0,@mustBeReal);
modParser.addAdditionalConstraint(@isscalar,"misalignmentAngle");

modParser.addOptionalParam("darkCountRate", 0, @(x) mustBeInRange(x, 0, 1));
modParser.addAdditionalConstraint(@isscalar,"darkCountRate");

modParser.addOptionalParam("visibilityAngle", 0, @mustBeReal)

modParser.parse(params);

params = modParser.Results;


lamb = params.lamb;
pz=params.pz;
visibilityAngle = params.visibilityAngle;
%% construct the signals
% We simulate the depolarization channel by assuming that the signals
% randomly misalign to the right or to the left by some fixed angle
% (visibilityAngle). This is equivalent to calculate the misalignments
% click patterns on both angles and take the average of them.
%% left
% we will use an intensity of 1 for now as we can scale that up after the
% fact.
signals = {Coherent.pauliCoherentState(1,1,1);... %H
    Coherent.pauliCoherentState(1,1,2);... %V
    Coherent.pauliCoherentState(1,2,1);... %D
    Coherent.pauliCoherentState(1,2,2)}; %A

% build the (sub) isometry transition matrix that represents the channel
% and Bob's measurement
transMat1 = simpleBB84LinearOpticsSetup(params.eta,params.misalignmentAngle+visibilityAngle,params.detectorEfficiency,1);
debugInfo.storeInfo("transMat1",transMat1);



probDetectorClickCon1 = zeros(numel(signals),size(transMat1,1),1);
expectationsCon1 = zeros(numel(signals),2^size(transMat1,1),1);

    %scale the signal states for the intensity
    signalsDecoy = cellfun(@(x) sqrt(lamb)*x,signals,"UniformOutput",false);
    [expectationsCon1(:,:,1), probDetectorClickCon1(:,:,1)] = simulateChannel(signalsDecoy,transMat1,params.darkCountRate);

debugInfo.storeInfo("detectorClickCon1",probDetectorClickCon1);

transMat2 = simpleBB84LinearOpticsSetup(params.eta,params.misalignmentAngle+visibilityAngle,params.detectorEfficiency,2);
debugInfo.storeInfo("transMat2",transMat2);



probDetectorClickCon2 = zeros(numel(signals),size(transMat2,1),numel(1));
expectationsCon2 = zeros(numel(signals),2^size(transMat2,1),numel(1));


    %scale the signal states for the intensity
    signalsDecoy = cellfun(@(x) sqrt(lamb)*x,signals,"UniformOutput",false);
    [expectationsCon2(:,:,1), probDetectorClickCon2(:,:,1)] = simulateChannel(signalsDecoy,transMat2,params.darkCountRate);

debugInfo.storeInfo("detectorClickCon2",probDetectorClickCon2);

left_expectationsCon = cat(2, pz*expectationsCon1, (1-pz)*expectationsCon2);
%% right
signals = {Coherent.pauliCoherentState(1,1,1);... %H
    Coherent.pauliCoherentState(1,1,2);... %V
    Coherent.pauliCoherentState(1,2,1);... %D
    Coherent.pauliCoherentState(1,2,2)}; %A

% build the (sub) isometry transition matrix that represents the channel
% and Bob's measurement
transMat1 = simpleBB84LinearOpticsSetup(params.eta,params.misalignmentAngle-visibilityAngle,params.detectorEfficiency,1);
debugInfo.storeInfo("transMat1",transMat1);



probDetectorClickCon1 = zeros(numel(signals),size(transMat1,1),1);
expectationsCon1 = zeros(numel(signals),2^size(transMat1,1),1);

    %scale the signal states for the intensity
    signalsDecoy = cellfun(@(x) sqrt(lamb)*x,signals,"UniformOutput",false);
    [expectationsCon1(:,:,1), probDetectorClickCon1(:,:,1)] = simulateChannel(signalsDecoy,transMat1,params.darkCountRate);

debugInfo.storeInfo("detectorClickCon1",probDetectorClickCon1);

transMat2 = simpleBB84LinearOpticsSetup(params.eta,params.misalignmentAngle-visibilityAngle,params.detectorEfficiency,2);
debugInfo.storeInfo("transMat2",transMat2);



probDetectorClickCon2 = zeros(numel(signals),size(transMat2,1),numel(1));
expectationsCon2 = zeros(numel(signals),2^size(transMat2,1),numel(1));


    %scale the signal states for the intensity
    signalsDecoy = cellfun(@(x) sqrt(lamb)*x,signals,"UniformOutput",false);
    [expectationsCon2(:,:,1), probDetectorClickCon2(:,:,1)] = simulateChannel(signalsDecoy,transMat2,params.darkCountRate);

debugInfo.storeInfo("detectorClickCon2",probDetectorClickCon2);

right_expectationsCon = cat(2, pz*expectationsCon1, (1-pz)*expectationsCon2);

% Average of both clicks patterns to model depolariztion errors
Final_expectationsCon = 1/2*(left_expectationsCon+right_expectationsCon);
newParams.expectationsConditional = Final_expectationsCon;
end

function transMat = simpleBB84LinearOpticsSetup(eta,misalignmentAngle,detectorEfficiency,rotation)

%% construct channel transition marix
%loss/transmittance
channelMat = Coherent.copyChannel(Coherent.transmittanceChannel(eta),2);

%misalignment rotation 
channelMat = Coherent.rotateStateZXY(misalignmentAngle,[0,0,1],"angleOnBlochSphere",false)*channelMat;


%% Build up Bob's detector transition matrix
% Each detector has the same efficiency so we can pull it right to the
% start.
detectorMat = Coherent.copyChannel(Coherent.transmittanceChannel(detectorEfficiency),2);

% Bob applies a rotation to convert A and D back to H and V for easier
% measurement
detectorMat = blkdiag(pauliBasis(rotation,false).')*detectorMat;


% we have to handle dark counts after we get the click probabilities for
% each detector.

transMat = detectorMat*channelMat;
end

function probDetectorClickCon = applyDarkCounts(probDetectorClickCon,darkCountRate)
probDetectorClickCon = 1-(1-probDetectorClickCon)*(1-darkCountRate);
end

function [detectorClickPatternCon, detectorClickCon] = simulateChannel(signals,transMat,darkCountRate)
    %Construct the independent detector click probabilities for each signal
    detectorClickCon = detectorClickProbabilities(signals,transMat);
    %simulate the effects of dark counts
    detectorClickCon  = applyDarkCounts(detectorClickCon,darkCountRate);
    %Construct all combinations of detector firing patterns from the
    %independent detectors.
    detectorClickPatternCon = detectorClickPatterns(detectorClickCon);
end



function probDetectorClickCon = detectorClickProbabilities(signals,transMat)

probDetectorClickCon = zeros(numel(signals),size(transMat,1));

for index = 1:numel(signals)
    bobsSignal = transMat*signals{index};
    probDetectorClickCon(index,:) = 1-Coherent.fockCoherentProb(zeros(size(transMat,1),1),bobsSignal,"combineModes",false);
end
end


function detectorClickPatternCon = detectorClickPatterns(probClickCon)
numSignals = size(probClickCon,1);
numDetectors = size(probClickCon,2);

detectorClickPatternCon = zeros(numSignals,2^numDetectors);

sizeDetectorPatterns = 2*ones(1,numDetectors);
clickProbSwitch = @(click, clickProb) (click==0).*(1-clickProb) + (click~=0).*clickProb;
for signalIndex = 1:numSignals
    for indexPat = 1:2^numDetectors
        patternVec = ind2subPlus(sizeDetectorPatterns,indexPat)-1;
        detectorClickPatternCon(signalIndex,indexPat) = prod(clickProbSwitch(patternVec,probClickCon(signalIndex,:)));
    end
end
end