qkdInput = NoPABAsymWCP_Preset();

%run the QKDSolver with this input
qkdInput.addFixedParameter("pz",0.1);
qkdInput.addOptimizeParameter("lamb",struct("lowerBound",0,"initVal",0.1,"upperBound",1));

lossdB = linspace(0,40,10);
lossEta = 10.^(-lossdB/10);
qkdInput.addScanParameter("eta",num2cell(lossEta));

qkdInput.addFixedParameter("misalignmentAngle",4.5 * pi/180);

qkdInput.addFixedParameter("visibilityAngle",0);


results = MainIteration(qkdInput);
%%
%save the results and preset to a file.
save("NPABWCP_Asym_loss.mat","results","qkdInput")
QKDPlot.simple1DPlot(qkdInput,results,"yScaleStyle","log")