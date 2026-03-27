qkdInput = NoPABFiniteWCP_Preset();

%run the QKDSolver with this input
qkdInput.addFixedParameter("pz",0.1);
qkdInput.addOptimizeParameter("lamb",struct("lowerBound",0,"initVal",0.1,"upperBound",0.3));

lossdB = 10;
lossEta = 10.^(-lossdB/10);
qkdInput.addFixedParameter("eta",lossEta);

qkdInput.addFixedParameter("misalignmentAngle",0);

qkdInput.addScanParameter("visibilityAngle",num2cell(linspace(0,0.25,11)));


results = MainIteration(qkdInput);
%%
%save the results and preset to a file.
save("Archive_BB84WCP_Asym_Depol_10dB.mat","results","qkdInput")


QKDPlot.plotParameters({results},"visibilityAngle","lamb","yScaleStyle","log")
QKDPlot.simple1DPlot(qkdInput,results,"yScaleStyle","log")