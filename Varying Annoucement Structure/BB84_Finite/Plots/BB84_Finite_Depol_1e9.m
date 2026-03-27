qkdInput = BB84_Finite_WCPPreset();

qkdInput.addFixedParameter("pz",0.1);
qkdInput.addOptimizeParameter("lamb",struct("lowerBound",0,"initVal",0.1,"upperBound",0.3));

lossdB = 10;
lossEta = 10.^(-lossdB/10);
qkdInput.addFixedParameter("eta",lossEta);

qkdInput.addFixedParameter("misalignmentAngle",0);

qkdInput.addScanParameter("visibilityAngle",num2cell(linspace(0,0.25,11)));
qkdInput.addFixedParameter("N", 1e9);

results = MainIteration(qkdInput);
save("BB84_Depol_Finite_1e9.mat","qkdInput","results");