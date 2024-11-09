function [GoldSeq] = generateGoldLfsrSequenceCA(n, ciVecA, ciVecB, a0VecA, a0VecB, G2Delay)
    % Generate a Gold sequence for GPS PRN codes using two LFSRs.
    % Adjusts the Gold sequence output to match the GPS PRN -1/1 format.

    % Initialize the LFSR states
    lfsrStateA = a0VecA;
    lfsrStateB = a0VecB;
    
    % Sequence Length
    m = 2^n - 1;
    GoldSeq = zeros(m, 1);  % Preallocate Gold sequence
    
    % Generate LFSR sequences and combine them into a Gold sequence
    for i = 1:m
        % Apply delay to the B sequence (G2) for the specified PRN
        delayedBitB = lfsrStateB(G2Delay);  % G2 delayed output
        
        % XOR the outputs of the two LFSRs
        GoldSeq(i) = mod(lfsrStateA(end) + delayedBitB, 2);
        
        % Update LFSR A
        feedbackA = mod(sum(lfsrStateA(ciVecA)), 2);  % XOR feedback for A
        lfsrStateA = [feedbackA; lfsrStateA(1:end-1)];  % Shift LFSR A state
        
        % Update LFSR B
        feedbackB = mod(sum(lfsrStateB(ciVecB)), 2);  % XOR feedback for B
        lfsrStateB = [feedbackB; lfsrStateB(1:end-1)];  % Shift LFSR B state
    end
    
    % Convert from 0/1 to -1/+1 format for GPS PRN compatibility
    GoldSeq = 2 * GoldSeq - 1;
end

