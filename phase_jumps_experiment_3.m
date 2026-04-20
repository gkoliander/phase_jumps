% Script to calculate empirical charge distribution
% of the STFT of white noise
% with a first Hermite window function
% with and without phase corrections

% Load required packages
pkg load ltfat; % version 2.6.0

% Load seed used in experiments
load('rand_seed.mat');
randn ('state', v);


% Indicate whether output is exported or plotted automatically
export_indicator = 0;
plot_indicator = 1;


% Parameter choices
subsamp = 128;
L = 32;            % can at most be subsamp/4
M = 4*subsamp*L;   % number of freq channels per dgt
a = 1;			       % Decimation factor

inv_delta = 128;   % inverse of final time and frequency resolution
tot_length = 4*L*inv_delta; % total signal length

iters = 5;         % number of realizations to average over


%Design window
  mytimeseq = circshift(((-tot_length/2):(tot_length/2-1)),tot_length/2,2);
  mywind = comp_hermite(1,mytimeseq/inv_delta*(2*pi)^(1/2));
  mywind = mywind/norm(mywind)*sqrt(inv_delta); % normalize to L2 norm 1
%comp_hermite is the same as
  %mywindalt2 = 2^(1/2)*pi^(-1/4)*(mytimeseq/time_subsamp).*exp(-1/2*(mytimeseq/time_subsamp).^2);
  %mywindalt2 = mywindalt2/norm(mywindalt2);

% Prepare arrays for estimated zero distribution
estposzerofct = zeros(4*L^2*inv_delta^2,1);
estposzerofct_stat = zeros(4*L^2*inv_delta^2,1);
estnegzerofct = zeros(4*L^2*inv_delta^2,1);
estnegzerofct_stat = zeros(4*L^2*inv_delta^2,1);

% Specify downsampling factor (note that sampling at a lower subsamp
% is not equivalent to downsampling as L is restricted by subsamp)
ds_exp = 2;
ds_factor = 2^ds_exp;
inv_delta_ds = inv_delta/ds_factor;
N = 2*L*inv_delta/ds_factor;
S = make_spiral(N);


for kk=1:(iters)
  printf ("Started realizaton %d of %d.\n",
        kk, iters);
  noise=1/sqrt(2)*(randn(M,1)+1i*randn(M,1)); % standard normal
  noise=noise/sqrt(subsamp);
  f = noise;



  % Compute coefficients
  c_orig = zeros(subsamp*inv_delta,tot_length);
  for f_ind = 1:(subsamp/M*inv_delta)
    for t_ind = 1:(inv_delta/subsamp)
      c_orig(f_ind:(subsamp/M*inv_delta):end,t_ind:(inv_delta/subsamp):end) = ...
           dgt(f'.*exp(-2*pi*1I*(f_ind-1)*(0:M-1)/(inv_delta*subsamp)),...
                  mywind((inv_delta/subsamp)-t_ind+1:(inv_delta/subsamp):end),a,M);
    end
  end


  [N2,M2] = size(c_orig);
  % Correct phase
  c_orig = c_orig.*exp(pi*1I/inv_delta^2*(0:N2-1)'*(0:M2-1));
  % for image and ignored phase correction shift center to M2/2 M2/2
  c_img = exp((-M2/2+1)*pi*1I/inv_delta^2*(0:M2-1)).*...
           c_orig.*...
           exp((M2/2+1)*pi*1I/inv_delta^2*(0:N2-1)');
  % shift center to M2/4 M2/4
  c_orig = exp((-M2/4+1)*pi*1I/inv_delta^2*(0:M2-1)).*...
           c_orig.*...
           exp((M2/4+1)*pi*1I/inv_delta^2*(0:N2-1)');


  boundary_width = 2*max(2, ceil(sqrt(inv_delta_ds)));
  f_samp_ds = c_orig(tot_length/4+1-ds_factor*boundary_width:ds_factor:3*tot_length/4+ds_factor*boundary_width,...
                     tot_length/4+1-ds_factor*boundary_width:ds_factor:3*tot_length/4+ds_factor*boundary_width);
  [N3,M3] = size(f_samp_ds);


  % Calculate phasechange at distance 2 delta
  [phasechange, maxchange] = ...
                  arg_change_calc(f_samp_ds, inv_delta_ds,2,boundary_width);

  % Calculate charges Argjumps
  charge_matrix = phasechange.* (maxchange<0.9*pi);
  clstsize = 6;
  [lmxposfilt,lmyposfilt,lmxnegfilt,lmynegfilt]= ...
                         resolve_clusters(f_samp_ds,charge_matrix,clstsize);
  % Remove zeros in the boundary region
  valid_charge = find((lmxposfilt > boundary_width).*(lmyposfilt > boundary_width)...
                .*(lmxposfilt < N3-boundary_width+1).*(lmyposfilt < M3-boundary_width+1));
  lmxposfilt = lmxposfilt(valid_charge)-boundary_width;
  lmyposfilt = lmyposfilt(valid_charge)-boundary_width;

  valid_charge = find((lmxnegfilt > boundary_width).*(lmynegfilt > boundary_width)...
                .*(lmxnegfilt < N3-boundary_width+1).*(lmynegfilt < M3-boundary_width+1));
  lmxnegfilt = lmxnegfilt(valid_charge)-boundary_width;
  lmynegfilt = lmynegfilt(valid_charge)-boundary_width;

  lmxposfilt_ts = lmxposfilt;
  lmyposfilt_ts = lmyposfilt;
  lmxnegfilt_ts = lmxnegfilt;
  lmynegfilt_ts = lmynegfilt;


  % Collect all zeros of Argjumps in a function
  % Go through the zeros in a spiral
  lmnegvec = zeros(4*L^2*inv_delta^2,1);
  lmnegvec(ds_factor^2*S(sub2ind([N,N], lmxnegfilt,lmynegfilt)))=1;
  lmposvec = zeros(4*L^2*inv_delta^2,1);
  lmposvec(ds_factor^2*S(sub2ind([N,N], lmxposfilt,lmyposfilt)))=1;

  estnegzerofct = estnegzerofct + cumsum(lmnegvec);
  estposzerofct = estposzerofct + cumsum(lmposvec);

  % Ignore phase correction
  f_samp_ds = c_img(tot_length/4+1-ds_factor*boundary_width:ds_factor:3*tot_length/4+ds_factor*boundary_width,...
                     tot_length/4+1-ds_factor*boundary_width:ds_factor:3*tot_length/4+ds_factor*boundary_width);
  [N3,M3] = size(f_samp_ds);
  % Calculate phasechange at distance 2 delta
  [phasechange_stat, maxchange] = ...
                  arg_change_calc_stat(f_samp_ds, inv_delta_ds,2,boundary_width);

  % Calculate charges Argjumps
  charge_matrix = phasechange_stat.* (maxchange<0.9*pi);
  clstsize = 6;
  [lmxposfilt,lmyposfilt,lmxnegfilt,lmynegfilt]= ...
                         resolve_clusters(f_samp_ds,charge_matrix,clstsize);
  % Remove zeros in the boundary region
  valid_charge = find((lmxposfilt > boundary_width).*(lmyposfilt > boundary_width)...
                .*(lmxposfilt < N3-boundary_width+1).*(lmyposfilt < M3-boundary_width+1));
  lmxposfilt = lmxposfilt(valid_charge)-boundary_width;
  lmyposfilt = lmyposfilt(valid_charge)-boundary_width;

  valid_charge = find((lmxnegfilt > boundary_width).*(lmynegfilt > boundary_width)...
                .*(lmxnegfilt < N3-boundary_width+1).*(lmynegfilt < M3-boundary_width+1));
  lmxnegfilt = lmxnegfilt(valid_charge)-boundary_width;
  lmynegfilt = lmynegfilt(valid_charge)-boundary_width;


  % Collect all zeros of Argjumps in a function
  % Go through the zeros in a spiral
  lmnegvec = zeros(4*L^2*inv_delta^2,1);
  lmnegvec(ds_factor^2*S(sub2ind([N,N], lmxnegfilt,lmynegfilt)))=1;
  lmposvec = zeros(4*L^2*inv_delta^2,1);
  lmposvec(ds_factor^2*S(sub2ind([N,N], lmxposfilt,lmyposfilt)))=1;

  estnegzerofct_stat = estnegzerofct_stat + cumsum(lmnegvec);
  estposzerofct_stat = estposzerofct_stat + cumsum(lmposvec);


end;

% Average over realizatons
estnegzerofct = estnegzerofct/iters;
estnegzerofct_stat = estnegzerofct_stat/iters;
estposzerofct = estposzerofct/iters;
estposzerofct_stat = estposzerofct_stat/iters;

  proposedzerofct = ((1:(4*L^2*inv_delta^2))/inv_delta^2)*5/3;
  estzerofct = estposzerofct + estnegzerofct;
  estzerofct_stat = estposzerofct_stat + estnegzerofct_stat;
  proposedchargefct = ((1:(4*L^2*inv_delta^2))/inv_delta^2);
  estchargefct = estposzerofct - estnegzerofct;
  estchargefct_stat = estposzerofct_stat - estnegzerofct_stat;

f_samp_ds = c_img(tot_length/4+1:ds_factor:3*tot_length/4,...
                       tot_length/4+1:ds_factor:3*tot_length/4);


if plot_indicator
  figure(31);
    plot(proposedzerofct);
    hold on;
      plot(estzerofct);
      plot(estzerofct_stat);
    hold off;

  figure(32);
    plot(proposedchargefct);
    hold on;
      plot(estchargefct);
      plot(estchargefct_stat);
    hold off;

  figure(33); imagesc(phasechange_stat);
  figure(34); imagesc(maxchange>0.9*pi);

  pos_diff = zeros(N);
  pos_diff(sub2ind([N,N], lmxposfilt_ts,lmyposfilt_ts)) += 1;
  pos_diff(sub2ind([N,N], lmxposfilt,lmyposfilt)) += -1;
  [pos_diff_x,pos_diff_y] = find(pos_diff==1);
  neg_diff = zeros(N);
  neg_diff(sub2ind([N,N], lmxnegfilt_ts,lmynegfilt_ts)) += 1;
  neg_diff(sub2ind([N,N], lmxnegfilt,lmynegfilt)) += -1;
  [neg_diff_x,neg_diff_y] = find(neg_diff==1);

  figure(51);   I = abs(f_samp_ds)/max(abs(f_samp_ds(:)));
  imagesc(I);
  hold on;
    scatter(pos_diff_y ,pos_diff_x,100, 'w',style='x');
    scatter(neg_diff_y ,neg_diff_x,100, 'w',style='o');
  hold off;

  [pos_diff_x,pos_diff_y] = find(pos_diff==-1);
  [neg_diff_x,neg_diff_y] = find(neg_diff==-1);
  figure(52);   I = abs(f_samp_ds)/max(abs(f_samp_ds(:)));
  imagesc(I);
  hold on;
    scatter(pos_diff_y ,pos_diff_x,100, 'w',style='x');
    scatter(neg_diff_y ,neg_diff_x,100, 'w',style='o');
  hold off;

end


% build and save export
export_x = (1:(4*L^2*inv_delta^2))'/inv_delta^2;
export_x = reshape(export_x,4*L^2*inv_delta^2/1024,1024);
export_x = mean(export_x,1)';
if export_indicator
  export_dir = 'D:\Experiments\';

  export_mat = [export_x,proposedchargefct(1:4*L^2*inv_delta^2/1024:end)'];
  save([export_dir, 'proposedchargefct_stat.dat'],'export_mat', '-ascii')
  export_mat = [export_x,proposedzerofct(1:4*L^2*inv_delta^2/1024:end)'];
  save([export_dir, 'proposedzerofct_stat.dat'],'export_mat', '-ascii')
  export_mat = [export_x,estzerofct(1:4*L^2*inv_delta^2/1024:end)];
  save([export_dir, 'estzerofct_ts.dat'],'export_mat', '-ascii')
  export_mat = [export_x,estzerofct_stat(1:4*L^2*inv_delta^2/1024:end)];
  save([export_dir, 'estzerofct_stat.dat'],'export_mat', '-ascii')
  export_mat = [export_x,estchargefct(1:4*L^2*inv_delta^2/1024:end)];
  save([export_dir, 'estchargefct_ts.dat'],'export_mat', '-ascii')
  export_mat = [export_x,estchargefct_stat(1:4*L^2*inv_delta^2/1024:end)];
  save([export_dir, 'estchargefct_stat.dat'],'export_mat', '-ascii')
end


