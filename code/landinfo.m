function out = landinfo(landmark_id, range_obs, x, y)
    out.landmark_id = landmark_id;
    out.observations = [range_obs, x, y];
end
