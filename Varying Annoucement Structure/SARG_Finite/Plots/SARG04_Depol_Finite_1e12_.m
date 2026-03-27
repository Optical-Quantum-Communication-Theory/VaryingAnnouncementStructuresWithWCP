qkdInput = SARG04_Finite_Adaptive_WCPPreset();
%run the QKDSolver with this input
lossdB = 10;
lossEta = 10.^(-lossdB/10);
qkdInput.addFixedParameter("eta",lossEta);
qkdInput.addOptimizeParameter("lamb",struct("lowerBound",0.4,"initVal",0.8,"upperBound",1.2));
qkdInput.addFixedParameter("misalignmentAngle",0);
qkdInput.addScanParameter("visibilityAngle",num2cell(linspace(0,0.25,11)));
qkdInput.addFixedParameter("N", 1e12);
results = MainIteration(qkdInput);

% save("Final_SARG04_Finite_Adaptive_WCP_scan_ncut7_1e12.mat",'qkdInput','results')
save("Archived_SARG04_Depol_Finite_1e12.mat","qkdInput","results");