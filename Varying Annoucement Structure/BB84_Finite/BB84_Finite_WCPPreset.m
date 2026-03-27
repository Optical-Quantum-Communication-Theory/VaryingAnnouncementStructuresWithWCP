function qkdInput = BB84_Finite_WCPPreset()
qkdInput = QKDSolverInput();

%% Parameters
% visibilityAngle: Amount of depolarization Bob's state experiences during
% transmission. For qubits, depolarization corresponds to a shrinking of the
% state on the Bloch sphere. With 0 depolarization a pure state remains on
% the surface of the Bloch sphere. This parameter is set separately in the
% plotting file.

% lamb: The mean photon number for weak coherent pulse. This parameter is set separately in the
% plotting file.

% eta: eta is the probability that a single photon that sucessfully passes
% through the channel and reaches Bob's detector. This parameter is set separately in the
% plotting file

% pz: The probability of Alice/Bob choosing the Z-basis for transmission / 
% measurement. For this BB84 protocol, if they don't choose the Z-basis
% then they transmit / measure in the X-basis. This parameter is set separately in the
% plotting file

% misalignmentAngle: The angle Bob's detectors are misaligned from Alice's.

% detectionEfficiency: The efficiency for Bob's detectors. In this project
% we set it to be 1
qkdInput.addFixedParameter("detectorEfficiency",1);

% we do not consider dark count in this project. This is always set to 0
qkdInput.addFixedParameter("darkCountRate",0);


qkdInput.addFixedParameter("pz", 0.9); 
%finite size parameters

epsSound = 1e-9; 
epsilon.AT = (1/4)*epsSound; %failure probability for parameter estimation 
epsilon.EC = (1/2)*epsSound; %failure probability for error-correction
epsilon.PA = (1/4)*epsSound; %failure probability for privacy amplification
qkdInput.addFixedParameter("epsilon", epsilon);

% probability of key gen round
qkdInput.addFixedParameter("pGen", 0.85);
% f: Efficiency of error correction. Real error correction protocols don't
% reach the Shannon limit. f is a scalar muliple, that scales up the amount
% of information leaked to fix a single bit. f=1, is the Shannon limit.
% f=1.16 is a more realistic value.
qkdInput.addFixedParameter("f", 1.2);

%% modules

descriptionModule = QKDDescriptionModule(@BB84_Finite_WCPDescriptionFunc);
qkdInput.setDescriptionModule(descriptionModule);

channelModule = QKDChannelModule(@finiteNPAB_Channel);

% The experiment setups are the same for the three protocols. Thus, we
% reuse the channel fron NPAB. We use this to simulate the expected behavior
% of the channel.
qkdInput.setChannelModule(channelModule);

keyRateModule = QKDKeyRateModule(@BB84_Finite_WCPKeyRateFunc);
qkdInput.setKeyRateModule(keyRateModule);

optimizerMod = QKDOptimizerModule(@coordinateDescentFunc,struct("verboseLevel",2,"maxIterations",1,"linearResolution",10),struct("verboseLevel",1));
qkdInput.setOptimizerModule(optimizerMod);

mathSolverOptions = struct();
mathSolverOptions.initMethod = 1; %closesest to maximally mixed.
mathSolverOptions.frankWolfeMethod = @FrankWolfe.vanilla;
mathSolverOptions.frankWolfeOptions = struct("maxIter",50,"maxGap",1e-6,"stopNegGap",true);
mathSolverOptions.linearConstraintTolerance = 1e-14;
mathSolverOptions.blockDiagonal = 0;
mathSolverMod = QKDMathSolverModule(@FW2StepSolver,mathSolverOptions);

qkdInput.setMathSolverModule(mathSolverMod);

%% global options
qkdInput.setGlobalOptions(struct("errorHandling",ErrorHandling.CatchWarn,"verboseLevel",1,"cvxSolver","Mosek", "cvxPrecision", "high"));
end