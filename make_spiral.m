function S = make_spiral(n)
%% Function generating a matrix of size n x n with elements from 1 to n^2
% starting at n/2,n/2 and counting upwards in a spiral
% Input:
% n -- dimension of matrix (must be even)
%
% Output:
% S -- matrix

   S = zeros(n, n);

    % Start coordinates
    r = n / 2;
    c = n / 2;

    % Start value
    val = 1;
    S(r, c) = val;

    % Step counters
    step_size = 1;
    step_count = 0;
    direction = 0; % 0=Right, 1=Down, 2=Left, 3=Up

    while val < n^2
        % Perform moves for the current leg of the spiral
        for k = 1:step_size
            if val >= n^2, break; end

            % Update coordinates based on direction
            switch direction
                case 0 % Right
                    c = c + 1;
                case 1 % Down
                    r = r + 1;
                case 2 % Left
                    c = c - 1;
                case 3 % Up
                    r = r - 1;
            end

            val = val + 1;
            S(r, c) = val;
        end

        % Change direction (0->1->2->3->0...)
        direction = mod(direction + 1, 4);

        % Increase step size every 2 legs (Right/Down=1, Left/Up=2, R/D=3...)
        step_count = step_count + 1;
        if mod(step_count, 2) == 0
            step_size = step_size + 1;
        end
    end
end
