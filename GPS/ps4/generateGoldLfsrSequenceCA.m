function [GoldSeq] = generateGoldLfsrSequenceCA(n,ciVecA,ciVecB,a0VecA,...
    a0VecB,G2tab)
%
% Generate a 1/0-valued linear feedback shift register (LFSR) sequence.
%
% INPUTS
%
% n ------ Number of stages in the linear feedback shift register.
%
% ciVec -- Nc-by-1 vector whose elements give the indices of the Nc nonzero
% connection elements. For example, if the characteristic polynomial
% of an LFSR is f(D) = 1 + D^2 + D^3, then ciVec = [2,3]’ or [3,2]’.
%
% a0Vec -- n-by-1 1/0-valued initial state of the LFSR, where a0Vec = [a(-1),
% a(-2), ..., a(-n)]’. In defining the initial LFSR state, a
% Fibonacci LFSR implementation is assumed.
%
% OUTPUTS
%
% GoldSeq -- m-by-1 vector whose elements are the 1/0-valued LFSR sequence
% corresponding to n, ciVec, and a0Vec, where m = 2^n - 1. If the
% sequence is a maximal-length sequence, then there is no
% repetition in the m sequence elements.
%


lfrsStateA=a0VecA;
lfrsStateB=a0VecB;
% Sequence Length
m=2^n-1;
% Preallocate GoldSeq
GoldSeq = zeros(m, 1);
for i=1:m
    G2         = mod(lfrsStateB(G2tab(1))+lfrsStateB(G2tab(2)),2);
    % shiftB     = circshift(lfrsStateB,G2Delay);
    GoldSeq(i) = mod(lfrsStateA(end)+G2,2);
    a=0;
    b=0;
    for jj = 1:length(ciVecA)
        a=mod(a+lfrsStateA(ciVecA(jj)),2);
    end

    for kk = 1:length(ciVecB)
        b=mod(b+lfrsStateB(ciVecB(kk)),2);
    end
    lfrsStateA=[a;lfrsStateA(1:end-1)];
    lfrsStateB=[b;lfrsStateB(1:end-1)];

end

