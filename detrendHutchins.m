function [u_de,u_ls] = detrendHutchins(vel,fs,nmin)

for i = 1:size(vel,2)
    % pad velocity data to avoid edge effects 
    v_raw = vel(:,i);
    v = [flipud(v_raw); v_raw; flipud(v_raw)];
    % first find the large-scale velocity data by filtering
    v_prime_temp = v - mean(v);
    Wn = 1/(60*nmin); 
    N = 2; % nmin minutes
    [b, a] = butter(N, 2*Wn/fs) ;
    v_prime_temp;
    v_ls = filtfilt(b, a, v_prime_temp);
    
    % find detrended velocity
    v_de (:,i) = v - v_ls;

    % find fluctuations 
    v_prime(:,i) = v_de (:,i) - mean(v_de (:,i)); 

    % extract only the middle piece back (get rid of padding)
    u_de(:,i) = v_de(length(v_raw)+1:end-length(v_raw),i); 
    u_ls(:,i) = v_ls(length(v_raw)+1:end-length(v_raw));
  %  u_prime(:,i) = v_prime((length(v_raw)+1:end-length(v_raw)));
    
end
end 