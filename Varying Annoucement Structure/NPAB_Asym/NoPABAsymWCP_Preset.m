function qkdInput = NoPABAsymWCP_Preset()
% Preset for Adaptive Length NPAB BB84 with WCPs
%% Parameters

% visibilityAngle: Amount of depolarization Bob's state experiences during
% transmission. For qubits, depolarization corresponds to a shrinking of the
% state on the Bloch sphere. With 0 depolarization a pure state remains on
% the surface of the Bloch sphere. This parameter is set separately in the
% plotting file.

% lamb: The mean photon number for weak coherent pulse. This parameter is set separately in the
% plotting file.

% eta: eta is the probability that a single photon that successfully passes
% through the channel and reaches Bob's detector. This parameter is set separately in the
% plotting file

% pz: The probability of Alice/Bob choosing the Z-basis for transmission / 
% measurement. For this BB84 protocol, if they don't choose the Z-basis
% then they transmit / measure in the X-basis. For SARG04, this is set to
% be 0.5

% misalignmentAngle: The angle Bob's detectors are misaligned from Alice's.
% This parameter is set separately in the plotting file

% ncutoff: The photon number cutoff. It will calculate key rate
% contribution for photon numbers less than ncutoff.

% detectionEfficiency: The efficiency for Bob's detectors. In this project
% we set it to be 1

qkdInput = QKDSolverInput();

%% Parameters


qkdInput.addFixedParameter("ncutoff", 3); 
qkdInput.addFixedParameter("lamb", 0.9);
% depolarization: unused for coherent simulation.
qkdInput.addFixedParameter("depolarization",0);
%Detector detection efficiency
qkdInput.addFixedParameter("detectorEfficiency",1);
% misalignmentAngle: The angle Bob's detectors are misaligned from Alice's.
qkdInput.addFixedParameter("misalignmentAngle",0);

%Channel loss
lossdB = 5; 
qkdInput.addFixedParameter("eta",10^(-lossdB/10));

%Z basis sending probability/beamsplitter ratio
% qkdInput.addOptimizeParameter("pz",struct("lowerBound",0.01,"upperBound",0.99,"initVal",0.8));
qkdInput.addFixedParameter("pz", 0.1);

%Detector Dark Count rate
qkdInput.addFixedParameter("darkCountRate", 0);

%Visibility 
qkdInput.addFixedParameter("visibilityAngle", 0.0);

% f: Efficiency of error correction
qkdInput.addFixedParameter("f",1.2);



%% modules

descriptionModule = QKDDescriptionModule(@NoPABAsymWCP_DescriptionFunc);
qkdInput.setDescriptionModule(descriptionModule);
channelModule = QKDChannelModule(@NPABWCPChannelFunc);
qkdInput.setChannelModule(channelModule);


keyRateModule = QKDKeyRateModule(@NoPABAsymWCP_KeyRateFunc);
qkdInput.setKeyRateModule(keyRateModule);

optimizerMod = QKDOptimizerModule(@coordinateDescentFunc,struct("verboseLevel",1, "maxIterations", 1, "linearResolution", 8),struct("verboseLevel",1));
qkdInput.setOptimizerModule(optimizerMod);

mathSolverOptions = struct();
mathSolverOptions.initMethod = 1; %closesest to maximally mixed.
mathSolverOptions.frankWolfeMethod = @FrankWolfe.vanilla;
mathSolverOptions.frankWolfeOptions = struct("maxIter",50,"maxGap",1e-8,"stopNegGap",true);
mathSolverOptions.linearConstraintTolerance = 1e-12;
mathSolverOptions.blockDiagonal = false;
mathSolverMod = QKDMathSolverModule(@FW2StepSolver,mathSolverOptions);
% mathSolverMod = QKDMathSolverModule(@FRGNSolver,mathSolverOptions);

qkdInput.setMathSolverModule(mathSolverMod);
qkdInput.setGlobalOptions(struct("errorHandling",ErrorHandling.CatchWarn,"verboseLevel",1,"cvxSolver","sedumi", "cvxPrecision", "high"));
%qkdInput.setGlobalOptions(struct("errorHandling", 3,"verboseLevel",1,"cvxSolver","sedumi", "cvxPrecision", "high"));
end