function cmap = circular_rgb_map(N, center, radius, normal_vec)
    % CIRCULAR_RGB_MAP Creates a circular colormap in 3D RGB space.
    %
    % Usage:
    %   cmap = circular_rgb_map(N, center, radius, normal_vec)
    %
    % Inputs:
    %   N          - Number of colors (default: 64)
    %   center     - [R, G, B] vector for the center of the circle (default: [0.5, 0.5, 0.5])
    %   radius     - Scalar radius of the circle (default: 0.4)
    %   normal_vec - [x, y, z] vector perpendicular to the circle's plane (default: [1, 1, 1])
    %
    % Returns:
    %   cmap       - N x 3 matrix of RGB values

    % --- Default Arguments ---
    if nargin < 1 || isempty(N), N = 64; end
    if nargin < 2 || isempty(center), center = [0.5, 0.5, 0.5]; end
    if nargin < 3 || isempty(radius), radius = 0.4; end
    if nargin < 4 || isempty(normal_vec), normal_vec = [1, 1, 1]; end

    % Ensure inputs are column/row vectors as needed
    center = reshape(center, 1, 3);
    normal_vec = reshape(normal_vec, 1, 3);

    % --- 1. Create Coordinate System (Basis Vectors) ---
    % Normalize the normal vector
    n = normal_vec / norm(normal_vec);

    % Find an arbitrary vector 'a' that is NOT parallel to n
    % We try [1,0,0]. If n is close to X-axis, we use [0,1,0].
    if abs(n(1)) > 0.9
        a = [0, 1, 0];
    else
        a = [1, 0, 0];
    end

    % Create vector u: perpendicular to n
    u = cross(n, a);
    u = u / norm(u);

    % Create vector v: perpendicular to n and u
    v = cross(n, u);
    v = v / norm(v);

    % --- 2. Generate Angles ---
    % Generate theta from 0 to 2*pi (exclusive of last point to close circle smoothly)
    theta = linspace(0, 2*pi, N+1)';
    theta(end) = []; % Remove the duplicate end point for a perfect cycle

    % --- 3. Calculate RGB Coordinates ---
    % Circle equation: P = C + r*cos(theta)*u + r*sin(theta)*v
    circle_points = center + radius .* (cos(theta)*u + sin(theta)*v);

    % --- 4. Clamp and Return ---
    % Ensure all values stay within legitimate RGB bounds [0, 1]
    cmap = max(0, min(1, circle_points));

    % Warning if clamping occurred (optional, good for debugging geometry)
    if any(circle_points(:) > 1) || any(circle_points(:) < 0)
        warning('circular_rgb_map:clipping', ...
        'The specified circle extends outside the RGB cube. Colors have been clamped.');
    end
end
