addpath 'g2o_wrapper'

function mse = mse_land(XL_guess, landmarks_ig, landmarks_gt)
	global LAND_NUM;
	mse = 0;

	for i=1:LAND_NUM
		id_guess = landmarks_ig(i).id;
		land_gt = searchById(landmarks_gt, id_guess);
		if isnumeric(land_gt)
			disp("Error: landmark is present in landmarks_ig but not in landmarks_gt");
			return
		endif
		mse += norm([
			XL_guess(:, i)(1) - land_gt.x_pose;
			XL_guess(:, i)(2) - land_gt.y_pose
		]);
	end
	mse /= LAND_NUM;
end

function mse = mse_pose(XR_guess, poses_ig, poses_gt)
	global POSE_NUM;
	mse = 0;

	for i=1:POSE_NUM
		id_guess = poses_ig(i).id;
		pose_gt = searchById(poses_gt, id_guess);
		if isnumeric(pose_gt)
			disp("Error: pose is present in poses_ig but not in poses_gt");
		endif
		mse += norm([
			XR_guess(:, :, i)(1, 3) - pose_gt.x;
			XR_guess(:, :, i)(2, 3) - pose_gt.y
		]);
	end
	mse /= POSE_NUM;
end
