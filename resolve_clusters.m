function [lmxposfilt,lmyposfilt,lmxnegfilt,lmynegfilt]= ...
                              resolve_clusters(f_samp,charge_matrix,clstsize)
  %% Function resolving clusters of zeros to a single zero
  % Input:
  % f_samp -- function values of function F on a grid as a matrix
  % charge_matrix -- -1/0/1 matrix indicating potential zero positions and
  %                  associated charge; same size as f_samp
  % clstsize -- desired separation of zeros with the same charge
  % Output:
  % lmxposfilt,lmyposfilt,lmxnegfilt,lmynegfilt
  %                             -- x/y coordinates of positively / negatively
  %                                charged zeros after resolving clusters;

  abss = abs(f_samp);
  % Determine whether positive, negative, or zero charge
  jacthresh = 0.5;
  locminmatpos = (charge_matrix>jacthresh);
  locminmatneg = (charge_matrix<-jacthresh);
  locminmatzero = (charge_matrix>-jacthresh).* (charge_matrix<jacthresh);
  [lmxpos,lmypos]=find(locminmatpos);
  [lmxneg,lmyneg]=find(locminmatneg);

  lmxposfilt = zeros(0,0);
  lmyposfilt = zeros(0,0);
  while (length(lmxpos) > 0)
    current_min = min(abss(sub2ind(size(abss),lmxpos,lmypos)));
    [x_min, y_min] = find(abss == current_min);
    intodel = (lmxpos<x_min+clstsize).*(lmxpos>x_min-clstsize);
    intodel = intodel.*(lmypos<y_min+clstsize).*(lmypos>y_min-clstsize);

    lmxposfilt = [lmxposfilt, x_min];
    lmyposfilt = [lmyposfilt, y_min];

    lmxpos(find(intodel)) = 0;
    lmypos = lmypos(find(lmxpos>0));
    lmxpos = lmxpos(find(lmxpos>0));

  endwhile

  lmxnegfilt = zeros(0,0);
  lmynegfilt = zeros(0,0);
  while (length(lmxneg) > 0)
    current_min = min(abss(sub2ind(size(abss),lmxneg,lmyneg)));
    [x_min, y_min] = find(abss == current_min);
    intodel = (lmxneg<x_min+clstsize).*(lmxneg>x_min-clstsize);
    intodel = intodel.*(lmyneg<y_min+clstsize).*(lmyneg>y_min-clstsize);
    lmxnegfilt = [lmxnegfilt, x_min];
    lmynegfilt = [lmynegfilt, y_min];

    lmxneg(find(intodel)) = 0;
    lmyneg = lmyneg(find(lmxneg>0));
    lmxneg = lmxneg(find(lmxneg>0));

  endwhile
