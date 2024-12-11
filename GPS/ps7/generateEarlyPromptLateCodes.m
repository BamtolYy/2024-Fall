function [codeEarly, codePrompt, codeLate] = generateEarlyPromptLateCodes(prn, ts, teml, fs, Nk)
    % INPUTS
    % prn -------- PRN number.
    % ts --------- Code start time (seconds).
    % teml ------- Early-minus-late spacing (chips).
    % fs --------- Sampling frequency.
    % Nk --------- Number of samples in one accumulation interval.
    %
    % OUTPUTS
    % codeEarly -- Oversampled early C/A code.
    % codePrompt - Oversampled prompt C/A code.
    % codeLate --- Oversampled late C/A code.

    % Chip interval for GPS L1 C/A
    Tc = 1e-3 / 1023; % Seconds per chip

    % Generate the C/A Code (Gold Code) for the PRN
    nStages = 10; % Number of stages for LFSR
    ciVec1 = [10, 3]'; % Feedback taps for G1
    ciVec2 = [10, 9, 8, 6, 3, 2]'; % Feedback taps for G2
    a0Vec1 = ones(nStages, 1); % Initial state for G1
    a0Vec2 = ones(nStages, 1); % Initial state for G2
    G2tab = [2, 6; 3, 7; 4, 8; 5, 9; 1, 9; 2, 10; 1, 8; 2, 9; 3, 10; ... 
             2, 3; 3, 4; 5, 6; 6, 7; 7, 8; 8, 9; 9, 10; 1, 4; 2, 5; ...
             3, 6; 4, 7; 5, 8; 6, 9; 1, 3; 4, 6; 5, 7; 6, 8; 7, 9; ...
             8, 10; 1, 6; 2, 7; 3, 8; 4, 9; 5, 10; 4, 10; 1, 7; 2, 8; ...
             4, 10];
    goldCode = generateGoldLfsrSequenceCA(nStages, ciVec1, ciVec2, a0Vec1, a0Vec2, G2tab(prn, :));
    goldCode = 2 * goldCode - 1; % Convert to +1/-1

    % Sampling interval in chips
    delChip = 1 / (fs * Tc);

    % Compute offsets for Early, Prompt, and Late codes
    delOffsetPrompt = 0; % Start time in chips
    delOffsetEarly = delOffsetPrompt - teml / 2; % Early code offset
    delOffsetLate = delOffsetPrompt + teml / 2; % Late code offset

    % Generate Oversampled Codes
    codePrompt = oversampleSpreadingCode(goldCode, delChip, -delOffsetPrompt, Nk, length(goldCode));
    codeEarly = oversampleSpreadingCode(goldCode, delChip, -delOffsetEarly, Nk, length(goldCode));
    codeLate = oversampleSpreadingCode(goldCode, delChip,-delOffsetLate, Nk, length(goldCode));
% figure,
% plot(codePrompt(1:100), 'b'); hold on;
% plot(codeEarly(1:100), 'r--'); hold on;
% plot(codeLate(1:100), 'o');
% legend('Prompt', 'Early', 'Late');
% title('Generated Early, Prompt, Late Codes');

end
