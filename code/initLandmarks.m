function [landmarks] = initLandmarks(poses, observations)

    % struct built to save relevant data for triangulating landmarks
    landmarks_info = struct();
    % land_id_mapping(i) will contain the index of info concerning
    % landmark with land_id=i in the struct landmarks_info
    % (thankfully land_ids are not sparse)
    land_id_mapping = [];
    % counting the found landmarks, just for convenience
    uniq_landmarks_found = 0;

    % fill landmarks_info with each landmark's observations
    for i=1:length(observations)
        curr_pose = searchById(poses, observations(i).pose_id);
        if isnumeric(curr_pose)  % pose of current observation set doesn't exist
            continue
        endif
        curr_x = curr_pose.x;
        curr_y = curr_pose.y;
        for j=1:length(observations(i).observation)
            land_id = observations(i).observation(j).id;
            range_obs = observations(i).observation(j).range;
            if land_id <= length(land_id_mapping) && land_id_mapping(land_id)
                % landmark is already observed
                landmarks_info(land_id_mapping(land_id)).observations(end+1, :) = [curr_x, curr_y, range_obs];
            else
                % new landmark
                uniq_landmarks_found += 1;
                land_id_mapping(land_id) = uniq_landmarks_found;
                landmarks_info(uniq_landmarks_found) = landinfo(land_id, range_obs, curr_x, curr_y);
            endif
        endfor
    endfor

    % initialize landmarks
    landmarks = struct();

    % for each landmark use the obtained observations to compute a likely position
    for i=1:length(landmarks_info)
        [x_pose, y_pose] = trilaterateLandmark(landmarks_info(i).observations);
        if x_pose != inf
            landmarks(end+1) = landmark(landmarks_info(i).landmark_id, [x_pose, y_pose]);
        endif
    endfor

    % remove first empty and meaningless field of the struct
    landmarks(1) = [];
endfunction

function out = landinfo(landmark_id, range_obs, x, y)
    out.landmark_id = landmark_id;
    out.observations = [x, y, range_obs];
endfunction

function [x, y] = trilaterateLandmark(observations)

    % Discard under-determined cases
    n = length(observations);
    if n < 4
        x = y = inf;
        return
    endif

    % (x-xi)^2 + (y-yi)^2 = ri^2 is the circle equation
    % => x^2 + y^2 - 2xix - 2yiy = ri^2 - xi^2 - yi^2
    % We have n_obs of such equations. In order to obtain a set of linear
    % equations (in the unknownks) we proceed subtracting the last equation
    % with i == n to all the others obtaining:
    % 2x(xn-xi) + 2y(yn-yi) = ri^2 - rn^2 - xi^2 + xn^2 - yi^2 + yn^2

    % Therefore, we need xn, yn and rn
    xn = observations(end, 1);
    yn = observations(end, 2);
    rn = observations(end, 3);

    % Initializing the matrices which will encode the n-1 equations
    A = zeros(n-1, 2);
    b = zeros(n-1, 1);

    for i=1:n-1
        xi = observations(i, 1);
        yi = observations(i, 2);
        ri = observations(i, 3);

        % Filling A and b with relevant coefficients
        A(i, :) = 2*[xn-xi yn-yi];
        b(i) = ri^2 - rn^2 - xi^2 + xn^2 - yi^2 + yn^2;
    endfor

    % Computing the pseudoinverse to deal with redundancy
    land_pose = pinv(A)*b;
    x = land_pose(1);
    y = land_pose(2);

    % Residue computation for sanity checking
    residue = 0;

    for i=1:n
        xi = observations(i, 1);
        yi = observations(i, 2);
        ri = observations(i, 3);

        % The larger the more inconsistent the set of equations
        residue += (norm([x y] - [xi yi]) - ri)/n;
    endfor

    % if residue in absolute value is higher than average of radii
    % discard the result
    if abs(residue) > mean(observations(:, 3))
        x = y = inf;
    endif
endfunction
