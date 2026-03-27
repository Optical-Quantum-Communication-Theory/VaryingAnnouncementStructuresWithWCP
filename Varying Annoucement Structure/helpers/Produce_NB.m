% Generate the matrix maps p(a,x,y) to p(a,b,x,y) based on SARG04 protocol

% The matrix has the form NB = sum p(b|a,y)|b>|a>|y><a|<y|otimes|x><x| (x is untouched)
% where S is the squashing map, mapping the raw click patterns y_unsquashed
% to y.

% Y:Bob's measurement
% A:Alice's announcement
% B:Bob's announcement
function [NB] = Produce_NB()
NB=0;
for b = 1:2 % b=1 (inconclusive announcement) b=2 (conclusive announcement)
    for a = 1:4 % 4 different groups of state Alice can annonce
        for y = 1:5 % 5 different outcomes that Bob can obtain after squashing
             %  To ensure the ordering:  |b>|a>|x>|y>
            NB = NB+prob_b_given_a_y(b,a,y)*kron(kron(kron(zket(2,b),zket(4,a)*zket(4,a)'),eye(4)),zket(5,y)*zket(5,y)');
    
        end
    end
end
    
    function [prob_b] = prob_b_given_a_y(b,a,y)
    groups = {[1,3],[2,3],[2,4],[1,4]}; % four possible announcements from Alice (1->H 2->V 3->D 4->A) [1,3] means Alices says she sends either H or D.
    if y ~=5 % cases where y is NOT vac/noclick
        if b == 2 % in the case Bob announces that he has a conclusive measurement (b=2)
            if ~(any(y == groups{a})) % Bob only makes this annoucement when in the case where Bob's measurement outcome is NOT one of Alice's announced states
                prob_b = 1; 
            else
                prob_b = 0; 
            end
        end    
    
        if b == 1 %in the case Bob announces that he has an inconclusive measurement(b=1)
            if any(y == groups{a}) % Bob only makes this annoucement when in the case where Bob's measurement outcome isone of Alice's announced states
                prob_b = 1;
            else
                prob_b = 0;
            end
        end 
    end

    if y ==5  % cases where y is vac/noclick, Bob only annouces vac/inconclusive (b=1)
        if b == 2
             prob_b = 0;
        end
       
    if b == 1
        prob_b = 1;
    end 
    end

    end      
end