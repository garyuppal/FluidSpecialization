

function bact = deterministicMutation(bact,numMutate,mdiff,binaryMutation)

    Nb =  length(bact(bact(:,3)>=0)); % number alive 
    
    Nm = min(numMutate,Nb); % mutate number of mutants upto number alive
    
    if binaryMutation == 1
        bact(1:Nm,3) = 0;
    else
        bact(1:Nm,3) = abs(bact(1:Nm,3) + mdiff*2*(rand(Nm,1)-1));
    end

end