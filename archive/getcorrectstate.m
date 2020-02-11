function statevector = getcorrectstate(seqMatrix, tempstate, nexttraj, trajfollowing)
%statevector = getcorrectstate(seqMatrix, state, nexttraj, trajfollowing)

[r c] = size(seqMatrix);
for k=1:(r-1) %length of sequence matrix /num rows
    if trajfollowing == 0 & nexttraj == 0
        if correctf == 1 || correctf == 2
            if seqMatrix(k,3) == 0 %if incorrect traj
                tempstate( seqMatrix(k,1):(seqMatrix(k+1,1)-1) ) = -2;
            end
        elseif correctf==0
            if seqMatrix(k,3) ~= 0 %if correct traj
                tempstate( seqMatrix(k,1):(seqMatrix(k+1,1)-1) ) = -2;
            end
        end
    elseif trajfollowing == 1 & k>1
        if seqMatrix(k,3) == 1 & seqMatrix(k-1,3) ==1 %current traj correct and past corret: C-C
            %leave
        elseif seqMatrix(k,3) == 1 & seqMatrix(k-1,3) ==0 %current traj correct and past incorret: IC-C
            tempstate( seqMatrix(k,1):(seqMatrix(k+1,1)-1) ) = -2;
        elseif seqMatrix(k,3) == 0 & seqMatrix(k-1,3) ==1 % C-IC
            tempstate( seqMatrix(k,1):(seqMatrix(k+1,1)-1) ) = -3;
        elseif seqMatrix(k,3) == 0 & seqMatrix(k-1,3) ==0 % IC-IC
            tempstate( seqMatrix(k,1):(seqMatrix(k+1,1)-1) ) = -4;
        end
    elseif nexttraj == 1
        if seqMatrix(k,3) == 1 & seqMatrix(k+1,3) ==1 %current traj correct and nextcorret: C-C
            %leave
        elseif seqMatrix(k,3) == 1 & seqMatrix(k+1,3) ==0 %current traj correct and next incorret: C-IC
            tempstate( seqMatrix(k,1):(seqMatrix(k+1,1)-1) ) = -2;
        elseif seqMatrix(k,3) == 0 & seqMatrix(k+1,3) ==1 % IC-C
            tempstate( seqMatrix(k,1):(seqMatrix(k+1,1)-1) ) = -3;
        elseif seqMatrix(k,3) == 0 & seqMatrix(k+1,3) ==0 % IC-IC
            tempstate( seqMatrix(k,1):(seqMatrix(k+1,1)-1) ) = -4;
        end
    end    
end
%last row
if seqMatrix(r,3) == 0 %if incorrect traj
    tempstate( seqMatrix(r,1):end ) = -2;
end

statevector = tempstate;
% statevector(undefind)=-1;

end