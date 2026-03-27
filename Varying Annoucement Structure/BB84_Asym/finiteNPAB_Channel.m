function [newParams, modParser]= finiteNPAB_Channel(params,options,debugInfo)
% Channel function for Regular and No Public Announcement of Basis BB84. Computes
% expected protocol statistics, conditioned on Alice's choice. 
% Conversion to joint expectation values done in Key Rate Module. 

% Input parameters:
% * pz: The probability that Alice measures in the Z-basis (for this protocol,
%   it's also the probability that Bob measures in the Z-basis as well). It
%   must be between 0 and 1.
% * lamb: The mean photon number for the Weak Coherent Pulse used in the
%   protocol
% * eta: eta is the probability that a single photon that successfully passes
% through the channel and reaches Bob's detector
% * misalignmentAngle: The angle Bob's detectors are misaligned from Alice's.
% * visibilityAngle (depolarization): Amount of depolarization Bob's state experiences during
% transmission. For qubits, depolarization corresponds to a shrinking of the
% state on the Bloch sphere. With 0 depolarization a pure state remains on
% the surface of the Bloch sphere.
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
modParser.addRequiredParam("lamb")

modParser.addOptionalParam("eta", 1, @(x) mustBeInRange(x, 0, 1));
modParser.addAdditionalConstraint(@isscalar,"eta");

modParser.addOptionalParam("detectorEfficiency", 1, @(x) mustBeInRange(x, 0, 1));
modParser.addAdditionalConstraint(@isscalar,"detectorEfficiency"); %not modeled, set to 1

modParser.addOptionalParam("misalignmentAngle",0,@mustBeReal);
modParser.addAdditionalConstraint(@isscalar,"misalignmentAngle");

modParser.addOptionalParam("darkCountRate", 0, @(x) mustBeInRange(x, 0, 1));
modParser.addAdditionalConstraint(@isscalar,"darkCountRate"); %not modeled, set to 1

modParser.addRequiredParam("pz", @(p) mustBeInRange(p, 0, 1)); 
modParser.addOptionalParam("visibilityAngle", 0, @mustBeReal)
modParser.parse(params);

params = modParser.Results;

lamb = params.lamb;
visibilityAngle = params.visibilityAngle;
debugInfo.storeInfo("lamb",lamb);
% construct the signals

% we will use an intensity of 1 for now as we can scale that up after the
% fact.
signals = {Coherent.pauliCoherentState(1,1,1);... %H
    Coherent.pauliCoherentState(1,1,2);... %V
    Coherent.pauliCoherentState(1,2,1);... %D
    Coherent.pauliCoherentState(1,2,2)}; %A

% build the (sub) isometry transition matrix that represents the channel
% and Bob's measurement except for the dark counts what must be handled
% later. Active BB84.
expConPlus = expMatVis(params, lamb, signals, visibilityAngle); 
expConMin = expMatVis(params, lamb, signals, -visibilityAngle); 
Final_expectationsCon = 1/2 * (expConPlus + expConMin); %Model depolarizing error (Eq.A22)
newParams.expectationsConditional = Final_expectationsCon;

end

function expectationsCon = expMatVis(params, lamb, signals, visAng)
    transMat1 = simpleBB84LinearOpticsSetup(params.eta,params.misalignmentAngle + visAng, params.detectorEfficiency,1);
    %debugInfo.storeInfo("transMat1",transMat1);
    
    
    
    probDetectorClickCon1 = zeros(numel(signals),size(transMat1,1),1);
    expectationsCon1 = zeros(numel(signals),2^size(transMat1,1),1);
    
        %scale the signal states for the intensity
        signalsDecoy = cellfun(@(x) sqrt(lamb)*x,signals,"UniformOutput",false);
        [expectationsCon1(:,:,1), probDetectorClickCon1(:,:,1)] = simulateChannel(signalsDecoy,transMat1,params.darkCountRate);
    
    %debugInfo.storeInfo("detectorClickCon1",probDetectorClickCon1);
    
    transMat2 = simpleBB84LinearOpticsSetup(params.eta,params.misalignmentAngle + visAng, params.detectorEfficiency,2);
    %debugInfo.storeInfo("transMat2",transMat2);
    
    
    
    probDetectorClickCon2 = zeros(numel(signals),size(transMat2,1),numel(1));
    expectationsCon2 = zeros(numel(signals),2^size(transMat2,1),numel(1));
    
    
        %scale the signal states for the intensity
        signalsDecoy = cellfun(@(x) sqrt(lamb)*x,signals,"UniformOutput",false);
        [expectationsCon2(:,:,1), probDetectorClickCon2(:,:,1)] = simulateChannel(signalsDecoy,transMat2,params.darkCountRate);
    
    %debugInfo.storeInfo("detectorClickCon2",probDetectorClickCon2);
    pz = params.pz; 
    expectationsCon = cat(2, pz*expectationsCon1, (1-pz)*expectationsCon2);
end

function transMat = simpleBB84LinearOpticsSetup(eta,misalignmentAngle,detectorEfficiency,rotation)

%% construct channel transition marix
%loss/transmittance
channelMat = Coherent.copyChannel(Coherent.transmittanceChannel(eta),2);

%misalignment rotation 
channelMat = Coherent.rotateStateZXY(misalignmentAngle,[0,0,1],"angleOnBlochSphere",false)*channelMat;

%simulate visibility, mixture of rotation statistics

%% Build up Bob's detector transition matrix
% Each detector has the same efficiency so we can pull it right to the
% start.
detectorMat = Coherent.copyChannel(Coherent.transmittanceChannel(detectorEfficiency),2);

% Bob applies a beam splitter to send signals to each detector basis setup
% detectorMat = Coherent.copyChannel(Coherent.singleInputBeamSplitter(pz),...
%     2,"weaveCopies",true)*detectorMat;

% Bob applies a rotation to convert A and D back to H and V for easier
% measurement
%detectorMat = blkdiag(pauliBasis(1,false).',pauliBasis(2,false).')*detectorMat;
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