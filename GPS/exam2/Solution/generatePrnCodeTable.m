function [prnCodeTable] = generatePrnCodeTable(fs, Ta)
N = 1023;
delChip = 1.023e6 / fs;
Ns = Ta * fs;

n = 10;
a0Vec = ones(10, 1);

ciVec_G1 = [3, 10];
G1 = generateLfsrSequence(n, ciVec_G1, a0Vec);

ciVec_G2 = [2, 3, 6, 8, 9, 10];
G2 = generateLfsrSequence(n, ciVec_G2, a0Vec);

delayTable = [5
    876
    17
    18
    139
    140
    141
    251
    252
    254
    255
    256
    257
    258
    469
    470
    471
    472
    473
    474
    509
    512
    513
    514
    515
    516
    859
    860
    861
    862
    863
    950
    947
    948
    950];

prnCodeTable = zeros(length(delayTable), Ns);

for ii = 1 : length(delayTable)
    G2i = [G2(end-delayTable(ii)+1 : end); G2(1 : end-delayTable(ii))];
    prnCode = bitxor(G1, G2i);
    prnCodeOS = oversampleSpreadingCode(sign(prnCode-0.5), delChip, 0, Ns, N);
    prnCodeTable(ii, :) = prnCodeOS;
end