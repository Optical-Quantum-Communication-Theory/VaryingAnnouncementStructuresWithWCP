qkdInput = NoPABFiniteWCP_Preset();
%run the QKDSolver with this input
lossdB = 10;
lossEta = 10.^(-lossdB/10);
qkdInput.addFixedParameter("eta",lossEta);
qkdInput.addOptimizeParameter("lamb",struct("lowerBound",0.01,"initVal",0.1,"upperBound",0.2));
qkdInput.addScanParameter("misalignmentAngle",num2cell(linspace(0,0.3,11)));
qkdInput.addFixedParameter("visibilityAngle",0);
qkdInput.addFixedParameter("N", 1e12);
results = MainIteration(qkdInput);

save("NPAB_Misalign_1e12.mat","qkdInput","results");