function[POVMsA] = Generate_POVMsA(n_cutoff)
% A function to generate POVMs for Alice after source replacement.
% POVMsA are arranged to correct order (based on the announcement structure) for convenience.
POVMsA = {};
for i = 1:4
POVM = 0;
    for n = 1:(n_cutoff+1)
        vec = kron(zket(n_cutoff+1,n),zket(4,i));
        POVM = POVM + vec*vec';
    end

POVMsA = [POVMsA,{1/2*POVM}];
end
POVMsA = {POVMsA{1},POVMsA{3},POVMsA{2},POVMsA{3},POVMsA{2},POVMsA{4},POVMsA{1},POVMsA{4}};
end