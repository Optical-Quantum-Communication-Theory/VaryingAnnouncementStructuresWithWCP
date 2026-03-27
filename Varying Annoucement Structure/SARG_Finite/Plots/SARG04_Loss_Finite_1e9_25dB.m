qkdInput = SARG04_Finite_Adaptive_WCPPreset();
% Some points in the plot fail due to numerics stability. Thus this
% function generate the failed points by changing the initial points for
% optimization a bit.
%run the QKDSolver with this input
lossdB = 25;
lossEta = 10.^(-lossdB/10);
qkdInput.addFixedParameter("eta",lossEta);
qkdInput.addOptimizeParameter("lamb",struct("lowerBound",0,"initVal",0.1,"upperBound",0.2));
qkdInput.addFixedParameter("misalignmentAngle",0);
qkdInput.addFixedParameter("visibilityAngle",0);
qkdInput.addFixedParameter("N", 1e9);
results = MainIteration(qkdInput);

% save("Final_SARG04_Finite_Adaptive_WCP_scan_ncut7_1e12.mat",'qkdInput','results')
save("Archived_SARG04_Loss_Finite_1e9_missingpt_at_25dB.mat","qkdInput","results");