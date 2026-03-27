qkdInput = BB84WCPPreset();

%run the QKDSolver with this input
qkdInput.addFixedParameter("pz",0.1);
qkdInput.addOptimizeParameter("lamb",struct("lowerBound",0,"initVal",0.1,"upperBound",0.3));

lossdB = 10;
lossEta = 10.^(-lossdB/10);
qkdInput.addFixedParameter("eta",lossEta);

qkdInput.addFixedParameter("visibilityAngle",0);

qkdInput.addScanParameter("misalignmentAngle",num2cell(linspace(0,0.6,21)));


results = MainIteration(qkdInput);
%%
%save the results and preset to a file.
save("BB84WCP_Asym_misalignment_10dB.mat","results","qkdInput")

QKDPlot.plotParameters({results},"misalignmentAngle","lamb","yScaleStyle","log")
QKDPlot.simple1DPlot(qkdInput,results,"yScaleStyle","log")