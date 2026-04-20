function [phasechange, max_change] = arg_change_calc_stat(f_samp,inv_delta,boxSize,bw)
%% Function calculating the phase change around each point
%    of a complex-valued function
% Input:
% f_samp -- function values of function F on a grid as a matrix,
% inv_delta -- number of grid points per unit length
% boxSize -- parameter;  distance at which charge is calculated
% bw -- additional samples as boundary, required to correct phase
%
% Output:
% phasechange -- matrix of calculated phase change, of the same size as f_samp
% max_change -- maximal phase difference in the calculation of phase change
%                note that the values of phasechange and max_change on a boundary
%                of boxSize elements are meaningless

  argc = arg(f_samp);
  [N,M] = size(argc);

%

  %% Approximate the path integral around each point
  phasechange=zeros(N,M);
  max_change=zeros(N,M);
  for ii=-boxSize:(boxSize-1)
    phase_bottom = mod(argc((1+2*boxSize):end,(1+boxSize+ii):(end-boxSize+ii)) ...
        - argc((1+2*boxSize):end,(1+boxSize+ii+1):(end-boxSize+ii+1))+pi, 2*pi)-pi;
    phase_top = mod(-argc((1):(end-2*boxSize),(1+boxSize+ii):(end-boxSize+ii)) ...
        + argc((1):(end-2*boxSize),(1+boxSize+ii+1):(end-boxSize+ii+1))+pi, 2*pi)-pi;
    phase_right = mod(-argc((1+boxSize+ii):(end-boxSize+ii),(1+2*boxSize):end) ...
        + argc((1+boxSize+ii+1):(end-boxSize+ii+1),(1+2*boxSize):end)+pi, 2*pi)-pi;
    phase_left = mod(argc((1+boxSize+ii):(end-boxSize+ii),(1):(end-2*boxSize)) ...
        - argc((1+boxSize+ii+1):(end-boxSize+ii+1),(1):(end-2*boxSize))+pi, 2*pi)-pi;
    phasechange((1+boxSize):(end-boxSize),(1+boxSize):(end-boxSize)) ...
      = phasechange((1+boxSize):(end-boxSize),(1+boxSize):(end-boxSize)) ...
        + phase_bottom + phase_top + phase_right + phase_left;
    max_change((1+boxSize):(end-boxSize),(1+boxSize):(end-boxSize)) = ...
                max(max_change((1+boxSize):(end-boxSize),(1+boxSize):(end-boxSize)), abs(phase_bottom));
    max_change((1+boxSize):(end-boxSize),(1+boxSize):(end-boxSize))  ...
              = max(max_change((1+boxSize):(end-boxSize),(1+boxSize):(end-boxSize)), abs(phase_top));
    max_change((1+boxSize):(end-boxSize),(1+boxSize):(end-boxSize))  ...
              = max(max_change((1+boxSize):(end-boxSize),(1+boxSize):(end-boxSize)), abs(phase_right));
    max_change((1+boxSize):(end-boxSize),(1+boxSize):(end-boxSize))  ...
              = max(max_change((1+boxSize):(end-boxSize),(1+boxSize):(end-boxSize)), abs(phase_left));
  end
phasechange = -phasechange/(2*pi);

