function qkdInput = BB84WCPPreset()
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
% This parameter is set separately in the plotting file


% detectionEfficiency: The efficiency for Bob's detectors. In this project
% we set it to be 1
qkdInput.addFixedParameter("detectorEfficiency",1);

% we do not consider dark count in this project. This is always set to 0
qkdInput.addFixedParameter("darkCountRate",0);



% f: Efficiency of error correction. Real error correction protocols don't
% reach the Shannon limit. f is a scalar muliple, that scales up the amount
% of information leaked to fix a single bit. f=1, is the Shannon limit.
% f=1.16 is a more realistic value. In the asymptotic case, we set this to
% 1
qkdInput.addFixedParameter("f",1);

%% modules
descriptionModule = QKDDescriptionModule(@BB84WCPDescriptionFunc);
qkdInput.setDescriptionModule(descriptionModule);

% The experiment setups are the same for the three protocols. Thus, we
% reuse the channel fron NPAB.
channelModule = QKDChannelModule(@finiteNPAB_Channel);
qkdInput.setChannelModule(channelModule);

% The key rate module contains the proof techniques to produce a safe lower
% bound on the key rate.
keyRateModule = QKDKeyRateModule(@BB84WCPKeyRateFunc);
qkdInput.setKeyRateModule(keyRateModule);

optimizerMod = QKDOptimizerModule(@coordinateDescentFunc,struct("verboseLevel",1,"maxIterations",1),struct("verboseLevel",1));
qkdInput.setOptimizerModule(optimizerMod);

mathSolverOptions = struct();
mathSolverOptions.initMethod = 1; %closesest to maximally mixed.
mathSolverOptions.linearConstraintTolerance = 1e-14;
mathSolverOptions.frankWolfeMethod = @FrankWolfe.vanilla;
mathSolverOptions.blockDiagonal = 1;
mathSolverOptions.frankWolfeOptions = struct("maxIter",50,"maxGap",1e-6,"stopNegGap",true);
mathSolverMod = QKDMathSolverModule(@FW2StepSolver,mathSolverOptions);
% mathSolverMod = QKDMathSolverModule(@FRGNSolver,mathSolverOptions);

qkdInput.setMathSolverModule(mathSolverMod);

%% global options
% options used everywhere. Many modules may have options given that
% overwrite these, but these are the basic options used everywhere.
% From the documentation on makeGlobalOptionsParser:
% Currently the global options are:
% * cvxSolver (default "SDPT3"): String naming the solver CVX should use.
% * cvxPrecision (default "high"): String that CVX uses to set the solver
%   precision.
% * verboseLevel (default 1): Non-negative integer telling the program how
%   much information it should display in the command window. 0, minimum; 1
%   basic information; 2, full details, including CVX output.
% * errorHandling (default 1): Integer 1, 2, or 3, detailing how the
%   program should handle run time errors. 1: catch and warn the user. The
%   key rate for the point is set to 0. 2: catch but don't warn the user.
%   The key rate for the point is set to 0. 3: don't catch the error and
%   let it up the stack.
qkdInput.setGlobalOptions(struct("errorHandling",ErrorHandling.CatchWarn,"verboseLevel",1,"cvxSolver","mosek", "cvxPrecision", "high"));
end
