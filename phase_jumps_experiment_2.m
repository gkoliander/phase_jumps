% Script to calculate empirical charge distribution
% of the STFT of white noise plus a Gaussian mean
% with a first Hermite window function

% Load required packages
pkg load ltfat;  % version 2.6.0

% Load seed used in experiments
load('rand_seed.mat');
randn ('state', v);


% Indicate whether output is exported or plotted automatically
export_indicator = 0;
plot_indicator = 1;

% Parameter choices
subsamp = 32;
L = 2;             % effective domain size [-L,L]^2, L can at most be subsamp/4
M = 4*subsamp*L;   % number of freq channels used per dgt
a = 1;			       % Decimation factor, should be 1 for consistency

inv_delta = 32;  % inverse of final time and frequency resolution
tot_length = 4*L*inv_delta; % total signal length

iters = 500;         % number of realizations to average over

%Design window
  mytimeseq = circshift(((-tot_length/2):(tot_length/2-1)),tot_length/2,2);
  mywind = comp_hermite(1,mytimeseq/inv_delta*(2*pi)^(1/2));
%comp_hermite is the same as
  %mywindalt2 = 2^(1/2)*pi^(-1/4)*(mytimeseq/inv_delta*(2*pi)^(1/2)).*exp(-1/2*(mytimeseq/inv_delta*(2*pi)^(1/2)).^2);

%Design mean signal
  mysignal = 0.2*comp_hermite(0,((-M/2):(M/2-1))/subsamp*(2*pi)^(1/2))'...
             .*exp((-2*L*(2*pi)^(-1/2))*2*pi*1I*((-M/2):(M/2-1))'/subsamp*(2*pi)^(1/2));
  c_sig = zeros(subsamp*inv_delta,tot_length);
  for f_ind = 1:(subsamp/M*inv_delta)
    for t_ind = 1:(inv_delta/subsamp)
      c_sig(f_ind:(subsamp/M*inv_delta):end,t_ind:(inv_delta/subsamp):end) = ...
           dgt(mysignal'.*exp(-2*pi*1I*(f_ind-1)*(0:M-1)/(inv_delta*subsamp)),...
                  mywind((inv_delta/subsamp)-t_ind+1:(inv_delta/subsamp):end),a,M);
    end
  end

  [N2,M2] = size(c_sig);
  % Correct phase
  c_sig = c_sig.*exp(pi*1I/inv_delta^2*(0:N2-1)'*(0:M2-1));
  % apply twisted shift so (M2/2 M2/2) corresponds to the origin
  c_sig = exp((-M2/2+1)*pi*1I/inv_delta^2*(0:M2-1)).*...
           c_sig.*...
           exp((M2/2+1)*pi*1I/inv_delta^2*(0:N2-1)');


% Prepare arrays for estimated zero distribution
% Results for argjumps algorithm
estposzerofct = zeros(4*L^2*inv_delta^2, 1);    %positively charged zeros
estnegzerofct = zeros(4*L^2*inv_delta^2, 1);    %negatively charged zeros
% Results for argjumps coarse algorithm
estposzerofct_c = zeros(4*L^2*inv_delta^2, 1);  %positively charged zeros
estnegzerofct_c = zeros(4*L^2*inv_delta^2, 1);  %negatively charged zeros

% Specify downsampling factor
ds_exp = 0;
ds_factor = 2^ds_exp;
inv_delta_ds = inv_delta/ds_factor;

% Construct vector S that indexes an N x N matrix from the center in a spiral
N = 2*L*inv_delta/ds_factor;
S = make_spiral(N);


%Start iterations
for kk=1:(iters)
  printf ("Started realizaton %d of %d.\n",
        kk, iters);


  %Generate white noise
  noise=1/sqrt(2)*(randn(M,1)+1i*randn(M,1));
  noise=noise/sqrt(inv_delta);
  %Add signal
  f = noise;

  % Compute random field values by applying dgts
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
  % apply twisted shift so (M2/2 M2/2) corresponds to the origin
  c_orig = exp((-M2/2+1)*pi*1I/inv_delta^2*(0:M2-1)).*...
           c_orig.*...
           exp((M2/2+1)*pi*1I/inv_delta^2*(0:N2-1)');

  %Add signal
  c_orig = c_orig + c_sig;
  c_img = c_orig(tot_length/4+1:3*tot_length/4,...
                       tot_length/4+1:3*tot_length/4);

  % apply twisted shift so (M2/4 M2/4) corresponds to the origin
  c_orig = exp((M2/4+1)*pi*1I/inv_delta^2*(0:M2-1)).*...
           c_orig.*...
           exp((-M2/4+1)*pi*1I/inv_delta^2*(0:N2-1)');



  boundary_width = 2*max(2, ceil(sqrt(inv_delta_ds)));
  f_samp_ds = c_orig(tot_length/4+1-ds_factor*boundary_width:ds_factor:3*tot_length/4+ds_factor*boundary_width,...
                     tot_length/4+1-ds_factor*boundary_width:ds_factor:3*tot_length/4+ds_factor*boundary_width);
  [N3,M3] = size(f_samp_ds);

  % Calculate Step 1 of Argjumps coarse
  ind_step1 = arg_jumps_coarse_step1(f_samp_ds,inv_delta_ds);

  % Calculate phasechange at distance delta**
  [phasechange_coarse, maxchange_coarse] = ...
                  arg_change_calc(f_samp_ds, inv_delta_ds,ceil(sqrt(inv_delta_ds)),boundary_width);

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

  % Calculate charges Argjumps coarse
    charge_matrix = phasechange.* ind_step1;
    clstsize = 5*ceil(sqrt(inv_delta_ds));
    [lmxposfilt_c,lmyposfilt_c,lmxnegfilt_c,lmynegfilt_c]= ...
                           resolve_clusters(f_samp_ds,charge_matrix,clstsize);
    % Remove zeros in the boundary region
    valid_charge = find((lmxposfilt_c > boundary_width).*(lmyposfilt_c > boundary_width)...
                  .*(lmxposfilt_c < N3-boundary_width+1).*(lmyposfilt_c < M3-boundary_width+1));
    lmxposfilt_c = lmxposfilt_c(valid_charge)-boundary_width;
    lmyposfilt_c = lmyposfilt_c(valid_charge)-boundary_width;

    valid_charge = find((lmxnegfilt_c > boundary_width).*(lmynegfilt_c > boundary_width)...
                  .*(lmxnegfilt_c < N3-boundary_width+1).*(lmynegfilt_c < M3-boundary_width+1));
    lmxnegfilt_c = lmxnegfilt_c(valid_charge)-boundary_width;
    lmynegfilt_c = lmynegfilt_c(valid_charge)-boundary_width;

  % Collect all zeros of Argjumps in a function
  % Go through the zeros in a spiral
    lmnegvec = zeros(4*L^2*inv_delta^2,1);
    lmnegvec(ds_factor^2*S(sub2ind([N,N], lmxnegfilt,lmynegfilt)))=1;
    lmposvec = zeros(4*L^2*inv_delta^2,1);
    lmposvec(ds_factor^2*S(sub2ind([N,N], lmxposfilt,lmyposfilt)))=1;

  estposzerofct(:,ds_exp+1) = estposzerofct(:,ds_exp+1) + cumsum(lmposvec);
  estnegzerofct(:,ds_exp+1) = estnegzerofct(:,ds_exp+1) + cumsum(lmnegvec);

  % Collect all zeros of Argjumps coarse in a function
  % Go through the zeros in a spiral
    lmnegvec = zeros(4*L^2*inv_delta^2,1);
    lmnegvec(ds_factor^2*S(sub2ind([N,N], lmxnegfilt_c,lmynegfilt_c)))=1;
    lmposvec = zeros(4*L^2*inv_delta^2,1);
    lmposvec(ds_factor^2*S(sub2ind([N,N], lmxposfilt_c,lmyposfilt_c)))=1;

  estposzerofct_c(:,ds_exp+1) = estposzerofct_c(:,ds_exp+1) + cumsum(lmposvec);
  estnegzerofct_c(:,ds_exp+1) = estnegzerofct_c(:,ds_exp+1) + cumsum(lmnegvec);

endfor

% Average over realizatons
estposzerofct = estposzerofct/iters;
estnegzerofct = estnegzerofct/iters;
estposzerofct_c = estposzerofct_c/iters;
estnegzerofct_c = estnegzerofct_c/iters;
estzerofct = estposzerofct + estnegzerofct;
estchargefct = estposzerofct - estnegzerofct;
estzerofct_c = estposzerofct_c + estnegzerofct_c;
estchargefct_c = estposzerofct_c - estnegzerofct_c;

%Calculate first intensity of the charge
  mywind = comp_hermite(1,mytimeseq/inv_delta*(2*pi)^(1/2));
  mywind = mywind/norm(mywind)*sqrt(inv_delta); % normalize to L2 norm 1

  mysignal = 0.2*comp_hermite(0,((-M/2):(M/2-1))/subsamp*(2*pi)^(1/2))'...
             .*exp((-2*L*(2*pi)^(-1/2))*2*pi*1I*((-M/2):(M/2-1))'/subsamp*(2*pi)^(1/2));
  c_sig = zeros(subsamp*inv_delta,tot_length);
  for f_ind = 1:(subsamp/M*inv_delta)
    for t_ind = 1:(inv_delta/subsamp)
      c_sig(f_ind:(subsamp/M*inv_delta):end,t_ind:(inv_delta/subsamp):end) = ...
           dgt(mysignal'.*exp(-2*pi*1I*(f_ind-1)*(0:M-1)/(inv_delta*subsamp)),...
                  mywind((inv_delta/subsamp)-t_ind+1:(inv_delta/subsamp):end),a,M);
    end
  end
  %Derivative D_1
  mywind_D1 = comp_hermite(2,mytimeseq/inv_delta*(2*pi)^(1/2));
  mywind_D1 = mywind_D1/norm(mywind_D1)*sqrt(inv_delta); % normalize to L2 norm 1

  D1_sig = zeros(subsamp*inv_delta,tot_length);
  for f_ind = 1:(subsamp/M*inv_delta)
    for t_ind = 1:(inv_delta/subsamp)
      D1_sig(f_ind:(subsamp/M*inv_delta):end,t_ind:(inv_delta/subsamp):end) = ...
           dgt(mysignal'.*exp(-2*pi*1I*(f_ind-1)*(0:M-1)/(inv_delta*subsamp)),...
                  mywind_D1((inv_delta/subsamp)-t_ind+1:(inv_delta/subsamp):end),a,M);
    end
  end
  D1_sig = sqrt(2)*D1_sig;

%Derivative D_2
  mywind_D2 = comp_hermite(0,mytimeseq/inv_delta*(2*pi)^(1/2));
  mywind_D2 = mywind_D2/norm(mywind_D2)*sqrt(inv_delta); % normalize to L2 norm 1

  D2_sig = zeros(subsamp*inv_delta,tot_length);
  for f_ind = 1:(subsamp/M*inv_delta)
    for t_ind = 1:(inv_delta/subsamp)
      D2_sig(f_ind:(subsamp/M*inv_delta):end,t_ind:(inv_delta/subsamp):end) = ...
           dgt(mysignal'.*exp(-2*pi*1I*(f_ind-1)*(0:M-1)/(inv_delta*subsamp)),...
                  mywind_D2((inv_delta/subsamp)-t_ind+1:(inv_delta/subsamp):end),a,M);
    end
  end

  sig2 = 1;

  rho_pm = 1/sig2*exp(-abs(c_sig.^2)/sig2).*(abs(D1_sig).^2-abs(D2_sig).^2+sig2);
  rho_pm = rho_pm(tot_length/4+1:3*tot_length/4,tot_length/4+1:3*tot_length/4);
  rho_pm(S(:)) = rho_pm(:);
  proposedchargefct_sig = cumsum(rho_pm(:))/inv_delta^2;

if plot_indicator
  figure(31);
  plot(estzerofct(:));

  figure(32);
  plot(estzerofct_c(:));

  figure(33);
  plot(proposedchargefct_sig);
  hold on;
    plot(estchargefct(:));
  hold off;

  figure(34);
  plot(proposedchargefct_sig);
  hold on;
    plot(estchargefct_c(:));
  hold off;


  % Specify figure size
  startx = 20;
  starty = 20;
  pwidth = 800;
  pheight = 800;


  figure(43, 'position',[startx,starty,pwidth,pheight]);
  I = (abs(c_img))/max(abs(c_img(:)));
  imagesc(I);
  set(gca, 'Position', [0 0 1 1]); axis off;
  hold on;
  scatter(lmyposfilt ,lmxposfilt,100, 'w',style='x', 'linewidth',1);
  scatter(lmynegfilt ,lmxnegfilt,100, 'w',style='o', 'linewidth',1);
  hold off;

  figure(44, 'position',[startx,starty,pwidth,pheight]);
  I = (abs(c_img))/max(abs(c_img(:)));
  imagesc(I);
  set(gca, 'Position', [0 0 1 1]); axis off;
  hold on;
  scatter(lmyposfilt_c ,lmxposfilt_c,100, 'w',style='x', 'linewidth',1);
  scatter(lmynegfilt_c ,lmxnegfilt_c,100, 'w',style='o', 'linewidth',1);
  hold off;

  figure(45, 'position',[startx,starty,pwidth,pheight]);
  my_map = circular_rgb_map(512, [0.5, 0.5, 0.5], 0.6, [1, -1, 1]);
  imagesc(arg(c_img)); colormap(my_map);
  set(gca, 'Position', [0 0 1 1]); axis off;
  hold on;
  scatter(lmyposfilt ,lmxposfilt,100, 'k',style='x', 'linewidth',1);
  scatter(lmynegfilt ,lmxnegfilt,100, 'k',style='o', 'linewidth',1);
  hold off;
end


% build and save export
if export_indicator
  export_dir = 'D:\Experiments\';

  export_x = (1:(4*L^2*inv_delta^2))'/inv_delta^2;
  export_x = reshape(export_x,4*L^2*inv_delta^2/1024,1024);
  export_x = mean(export_x,1)';
  export_mat = [export_x,proposedchargefct_sig(1:4*L^2*inv_delta^2/1024:end)];
  save([export_dir, 'proposedchargefct_sig.dat'],'export_mat', '-ascii')

    export_y = (estchargefct(1:end));
    export_y = reshape(export_y,4*L^2*inv_delta^2/1024,1024);
    export_y = mean(export_y,1)';
    export_mat = [export_x,export_y];
    save([export_dir, 'estchargefct_mean', num2str(1),'.dat'],'export_mat', '-ascii')


    export_y = (estzerofct(1:end));
    export_y = reshape(export_y,4*L^2*inv_delta^2/1024,1024);
    export_y = mean(export_y,1)';
    export_mat = [export_x,export_y];
    save([export_dir, 'estzerofct_mean', num2str(1),'.dat'],'export_mat', '-ascii')

    export_y = (estchargefct_c(1:end));
    export_y = reshape(export_y,4*L^2*inv_delta^2/1024,1024);
    export_y = mean(export_y,1)';
    export_mat = [export_x,export_y];
    save([export_dir, 'estchargefct_c_mean', num2str(1),'.dat'],'export_mat', '-ascii')

    export_y = (estzerofct_c(1:end));
    export_y = reshape(export_y,4*L^2*inv_delta^2/1024,1024);
    export_y = mean(export_y,1)';
    export_mat = [export_x,export_y];
    save([export_dir, 'estzerofct_c_mean', num2str(1),'.dat'],'export_mat', '-ascii')
end
