function ind_step1 = arg_jumps_coarse_step1(f_samp,inv_delta)
  %% Function calculating Step 1 of Argjumps coarse:
  % Find points w such that all points mu, mu' at distance deltastar and with
  % |mu-mu'|=1/inv_delta satisfy
  % |F_w(mu)| > 2 |F_w(mu)- F_w(mu')|
  % Input:
  % f_samp -- function values of function F on a grid as a matrix
  % inv_delta -- inverse of delta (distance between samples of F)
  % Output:
  % ind_step1 -- binary matrix of the same size as f_samp indicating samples that
  %              passed Step 1 of Argjumps coarse

  abss = abs(f_samp);
  [N,M] = size(f_samp);
  ind_step1 = ones(N,M);

  boxSize = ceil(sqrt(inv_delta));
  for ii=-boxSize:boxSize-1
    test = ( circshift(abss,[ii,boxSize]) > ...
                   2*abs(exp(-pi*1I/inv_delta^2*(0:M-1)).*circshift(f_samp,[ii+1,boxSize]) - ...
                         circshift(f_samp,[ii,boxSize])) );
    ind_step1=ind_step1.*test;
    test = ( circshift(abss,[ii,-boxSize]) > ...
                   2*abs(exp(-pi*1I/inv_delta^2*(0:M-1)).*circshift(f_samp,[ii+1,-boxSize]) - ...
                         circshift(f_samp,[ii,-boxSize])) );
    ind_step1=ind_step1.*test;
  end

  for jj=-boxSize:boxSize-1
    test = ( circshift(abss,[boxSize,jj]) > ...
                   2*abs(circshift(f_samp,[boxSize,jj+1]).*exp(pi*1I/inv_delta^2*(0:N-1)') - ...
                         circshift(f_samp,[boxSize,jj])) );
    ind_step1=ind_step1.*test;
    test = ( circshift(abss,[-boxSize,jj]) > ...
                   2*abs(circshift(f_samp,[-boxSize,jj+1]).*exp(pi*1I/inv_delta^2*(0:N-1)') - ...
                         circshift(f_samp,[-boxSize,jj])) );
    ind_step1=ind_step1.*test;
  end
