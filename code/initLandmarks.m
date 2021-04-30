function [landmarks] = initLandmarks(poses, observations)

    % struct built to save relevant data for triangulating landmarks
    landmarks_info = struct();
    % land_id_mapping(i) will contain the index of info concerning
    % landmark with land_id=i in the struct landmarks_info
    % (thankfully land_ids are not sparse)
    land_id_mapping = [];
    % counting the found landmarks, just for convenience
    uniq_landmarks_found = 0;
    % poses(i) is mapped to id 1099 + i
    % in particular, we can obtain pose = poses(pose_id - poses_id_offset);
    poses_id_offset = 1099;

    for i=1:length(observations)
        curr_pose_id = observations(i).pose_id;
        curr_x = poses(curr_pose_id - poses_id_offset).x;
        curr_y = poses(curr_pose_id - poses_id_offset).y;
        for j=1:length(observations(i).observation)
            land_id = observations(i).observation(j).id;
            range_obs = observations(i).observation(j).range;
            if land_id <= length(land_id_mapping) && land_id_mapping(land_id)
                % landmark is already observed
                landmarks_info(land_id_mapping(land_id)).observations(end+1, :) = [range_obs, curr_x, curr_y];
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

    % now, for each landmark use the obtained observations to compute a likely position
    for i=1:length(landmarks_info)
        [x_pose, y_pose] = trilaterateLandmark(landmarks_info(i).observations);
        if x_pose != null
            landmarks(end+1) = landmark(landmarks_info(i).landmark_id, [x_pose, y_pose]);
        endif
    endfor

    % remove first empty and meaningless field of the struct
    landmarks = landmarks(2:end);

end

function out = landinfo(landmark_id, range_obs, x, y)
    out.landmark_id = landmark_id;
    out.observations = [range_obs, x, y];
end

function [x, y] = trilaterateLandmark(observations, res_threshold)

    % Discard under-determined cases
    n = length(observations);
    if n < 4
        x = null;
        y = null;
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
        % Therefore, we need xn, yn and rn
        xi = observations(i, 1);
        yi = observations(i, 2);
        ri = observations(i, 3);

        % Filling A and b with relevant coefficients
        A(i, :) = 2[xn-xi yn-yi];
        b(i) = ri^2 - rn^2 - xi^2 + xn^2 - yi^2 + yn^2;
    endfor

    % Computing the pseudoinverse to deal with redundancy
    land_pose = pinv(A)*b;
    x = land_pose[1];
    y = land_pose[2];

    % Residue for sanity checking
    residue = 0;

    for i=1:n
        % Therefore, we need xn, yn and rn
        xi = observations(i, 1);
        yi = observations(i, 2);
        ri = observations(i, 3);

        % The larger the more inconsistent the set of equations
        residue += (norm([x y] - [xi yi]) - ri)/n;
    endfor

    residue

    if residue > res_threshold
        x = null;
        y = null;
    endif

end
