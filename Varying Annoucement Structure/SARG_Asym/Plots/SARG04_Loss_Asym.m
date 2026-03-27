lossdB = linspace(0,40,41);
lossEta = 10.^(-lossdB/10);
delta = 0.2; %range of the intensity you want to search around last optimal intensity
%initializing

keyrate = [];
mu = [];
loss = [];
full_results = {};
%initial value and initial bounds for lamb
lamb = 1.2;
lamb_up = lamb+delta;
lamb_low = lamb-delta;

for i = 1: numel(lossdB)
    eta = lossEta(i);
    qkdInput = SARG04WCP2Preset();
    qkdInput.addFixedParameter("eta",eta);
    qkdInput.addOptimizeParameter("lamb",struct("lowerBound",lamb_low, ...
        "initVal",lamb,"upperBound",lamb_up));
    qkdInput.addFixedParameter("misalignmentAngle",0);
    qkdInput.addFixedParameter("visibilityAngle",0);
    qkdInput.addFixedParameter("pz",0.5);
    try
        results = MainIteration(qkdInput);
    catch ME
     
    end 
    %get the intensity and key
    intensity = results.currentParams.lamb;
    key = [results.keyRate];

    %store the results
    keyrate = [keyrate;key];
    mu = [mu,intensity];
    full_results = [full_results,{results}];
    loss = [loss,lossdB(i)];
    %updating the bounds and initial point for next itration
    lamb = intensity;
    lamb_up = intensity*(1+ delta);
    lamb_low = max(0,intensity*(1- delta)); 

    %check consecutive zero keyrate
    if sum(keyrate<0)>2
        break
    end
end

%%
save("SARG04_Loss_Asym.mat","keyrate","mu","loss","full_results")

%% rerun points that fails
lossEta = 10.^(-21/10);
 qkdInput = SARG04WCP2Preset();
    qkdInput.addFixedParameter("eta",lossEta);
    qkdInput.addOptimizeParameter("lamb",struct("lowerBound",0.2, ...
        "initVal",0.24,"upperBound",0.28));
    qkdInput.addFixedParameter("misalignmentAngle",0);
    qkdInput.addFixedParameter("visibilityAngle",0);
    qkdInput.addFixedParameter("pz",0.5);
point_results = MainIteration(qkdInput);
