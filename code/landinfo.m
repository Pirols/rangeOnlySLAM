function out = landinfo(landmark_id, range_info, pose_id)
    out.landmark_id = landmark_id;
    out.observations = [range_info, pose_id];
end
