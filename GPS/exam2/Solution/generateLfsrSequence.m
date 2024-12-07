function [lfsrSeq] = generateLfsrSequence(n,ciVec,a0Vec)
% [lfsrSeq] = generateLfsrSequence(n,ciVec,a0Vec)
%
% Generate a 1/0−valued linear feedback shift register (LFSR) sequence.
%
% INPUTS
%
% n −−−−−− Number of stages in the linear feedback shift register.
%
% ciVec −− Nc−by−1 vector whose elements give the indices of the Nc nonzero
%          connection elements. For example, if the characteristic polynomial
%          of an LFSR is f(D) = 1 + Dˆ2 + Dˆ3, then ciVec = [2,3]' or [3,2]'.
%
% a0Vec −− n−by−1 1/0−valued initial state of the LFSR, where a0Vec = [a(−1),
%          a(−2), ..., a(−n)]'. In defining the initial LFSR state, a
%          Fibonacci LFSR implementation is assumed.
%
% OUTPUTS
%
% lfsrSeq −− m−by−1 vector whose elements are the 1/0−valued LFSR sequence
%            corresponding to n, ciVec, and a0Vec, where m = 2ˆn − 1. If the
%            sequence is a maximal−length sequence, then there is no
%            repetition in the m sequence elements.
%

m = 2^n - 1;
lfsrSeq = zeros(m, 1);

cVec = zeros(n, 1);
cVec(ciVec) = 1;

aVec = a0Vec;
for ii = 1 : m
    lfsrSeq(ii) = aVec(end);
    axci = aVec .* cVec;
    if mod(sum(axci), 2)
        aVec = [1; aVec(1 : end-1)];
    else
        aVec = [0; aVec(1 : end-1)];
    end
end