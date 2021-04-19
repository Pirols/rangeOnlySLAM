function [landmarks] = initLandmarks(poses, observations)

    landmarks_info = struct();
    % land_id_mapping(i) will contain the index of info concerning
    % landmark with land_id=i in the struct landmarks_info
    land_id_mapping = [];
    % just for convenience
    uniq_landmarks_found = 0;

    for i=1:length(observations)
        curr_pose_id = observations(i).pose_id;
        for j=1:length(observations(i).observation)
            land_id = observations(i).observation(j).id;
            range_obs = observations(i).observation(j).range;
            if land_id <= length(land_id_mapping) && land_id_mapping(land_id)
                % landmark is already observed
                landmarks_info(land_id_mapping(land_id)).observations(end+1, :) = [range_obs, curr_pose_id];
            else
                % new landmark
                uniq_landmarks_found += 1;
                land_id_mapping(land_id) = uniq_landmarks_found;
                landmarks_info(end+1) = landinfo(land_id, range_obs, curr_pose_id);
            endif
        endfor
    endfor

    % removing empty first element of struct
    landmarks_info = landmarks_info(2:end);

    % poses(i) is mapped to id 1099 + i
    % in particular, we can obtain pose = poses(pose_id - poses_id_offset);
    poses_id_offset = 1099;

    % initialize landmarks
    landmarks = landmarks_info;

    % now, for each landmark use the obtained observations to compute a likely position
end
