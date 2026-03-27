% WARNING: It takes more than an hour to generate the plots since the dim of this
% problem is large
qkdInput = SARG04WCP2Preset();
lossdB = 10;
lossEta = 10.^(-lossdB/10);
qkdInput.addFixedParameter("eta",lossEta);
qkdInput.addOptimizeParameter("lamb",struct("lowerBound",0.4,"initVal",0.8,"upperBound",1.2));
qkdInput.addFixedParameter("misalignmentAngle",0);
qkdInput.addScanParameter("visibilityAngle",num2cell(linspace(0,0.25,11)));
results = MainIteration(qkdInput);

% save("Final_SARG04_Finite_Adaptive_WCP_scan_ncut7_1e12.mat",'qkdInput','results')
save("Archived_SARG04_Depol_Asym.mat","qkdInput","results");