% Generate the matrix maps p(a,b,x,y) to p(r,a,b,x,y) based on SARG04 protocol
% R:keybit
% X:Alice's measurement
% Y:Bob's measurement
% A:Alice's announcement
% B:Bob's announcement
% The matrix has the form NB = sum p(r|a,b,x,y)|r>|a>|b>|x><a|<b|<x|otimes|y><y| (a,y are untouched)
function[G] = Produce_G()
G = 0;
for r = 1:3 % key0 (r=1), key1 (r=2), discard (r=3)
    for b = 1:2 % b=1 (inconclusive announcement) b=2 (conclusive announcement)
        for x = 1:4 % 4 different states Alice can send
% To ensure the ordering: |r>|b>|a>|x>|y>
            G  = G + prob_r_given_b_x(r,b,x)*kron(kron(kron(kron(zket(3,r),zket(2,b)*zket(2,b)'),eye(4)),zket(4,x)*zket(4,x)'),eye(5));
        end
    end
end

    function [prob] = prob_r_given_b_x(r,b,x)
        basis={[1,2],[3,4]}; % there are two bases
    if r==1 % Alice maps x to key0 (r=1) if Bob annouces an conclusive result (b=2) and she actually sends out a state in Z basis
        if b==2
            if any(x==basis{1}) % Z basis
                prob = 1;
            else
                prob = 0;
            end
        else 
            prob = 0;
        end
    end
    
    if r==2 %Alice maps x to key1 (r=2) if Bob annouces an conclusive result (b=2) and she actually sends out a state in X basis
        if b==2
            if any(x==basis{2})
                prob = 1;
            else
                prob = 0;
            end
        else 
            prob = 0;
        end
    end
    
    if r==3 % the key gets discarded if Bob announces an inconclusive result (b=1) 
        if b==1
            prob = 1;
        else
            prob = 0;
        end
    end
    end
end