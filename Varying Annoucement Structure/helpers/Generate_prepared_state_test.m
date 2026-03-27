function [rhoA,POVMsA,dimA] = Generate_prepared_state_test(lamb,n_cutoff,pz)
%% A function to calculate rhoA, POVMA and simplify it.
%get the correct dimA and dim_signal
dimA = (n_cutoff+1)*4;
if n_cutoff<=4
dim_signal = 4+(n_cutoff+1)*(n_cutoff)/2; 
% # of flag + dim of signals (vac is 1 dim, single photon is 2 dim ...)
% we add up to (exlucding) signals with n_cutoff number of photons)
else
dim_signal = 10+4*(n_cutoff-3);
end

%First we obtain the original source replaced POVMsA and state rhoAA'
POVMsA = Generate_POVMsA(n_cutoff); % becasue each state sending out has 2 possible announcements
pa = [pz/2,pz/2,(1-pz)/2,(1-pz)/2]; % probability of send each state

% the source replaced state
psiAAprime = 0;

for i = 1:4
signal=0;
states_all = {};
counts = 1;
    for pn = 1:(n_cutoff-1)
    [H,V,plus,minus]=Generate_States_test(pn);
    states_at_pn = {H,V,plus,minus};

    states_at_pn_i = [zeros(counts,1);states_at_pn{i}];
    counts = counts + min(4,(1+pn));
    states_at_pn_i = [states_at_pn_i;zeros(dim_signal-counts,1)];

    states_all = [states_all,{states_at_pn_i}];

    end
states = [zket(dim_signal,1),states_all,zket(dim_signal,dim_signal-(4-i))];

    for n=1:numel(states)
          signal = signal+ sqrt(pa(i))*sqrt(pois(n-1,lamb,n_cutoff))*(kron(zket(numel(states),n),zket(4,i))*states{n}');
    end
psiAAprime = psiAAprime + signal;
end


% Use Schidmt decomposition to cut dimensions of the source replaced state
[U, S, V] = svd(psiAAprime);
% getting new POVMsA
PiA = 0;
dimnew  = size(V,1);
dimA = dimnew;
for index = 1:dimnew
    PiA = PiA + zket(dimnew,index)*U(:,index)';
end
for i = 1:8
    POVMsA{i} = PiA*POVMsA{i}*PiA';
end
% getting new share state rhoAA'
newrhoAAprime = 0;
for index = 1:dimnew
    newrhoAAprime = newrhoAAprime + S(index,index)* kron(zket(dimnew,index),V(:,index));
end
rhoA = PartialTrace(newrhoAAprime*newrhoAAprime',2,[dimnew,dimnew]);
rhoA = (rhoA+rhoA')/2;


%% function to calculate pois pdf.
function[prob] = pois(x,lamb,n_cutoff)
if x<n_cutoff
    prob = lamb^x*exp(-lamb)/factorial(x);
else
    sum=0;
    for number = 0:(n_cutoff-1)
        sum = sum + pois(number,lamb,n_cutoff);
    end
    prob = 1-sum;
end
end
end