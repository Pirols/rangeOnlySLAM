addpath 'g2o_wrapper'

function mse = mse_land(lands_ground, lands_guess)
	mse = 0;
	n_lands = length(lands_guess);

	for i=1:n_lands
		id_guess = lands_guess(i).id;
		land_gt = searchById(lands_ground, id_guess);
		mse += norm([
			lands_guess(i).x_pose - land_gt.x_pose;
			lands_guess(i).y_pose - land_gt.y_pose
		]);
	end
	mse /= n_lands;
end

function mse = mse_pose(poses_ground, poses_guess)
	mse = 0;
	n_poses = length(poses_guess);

	for i=1:n_poses
		id_guess = poses_guess(i).id;
		pose_gt = searchById(poses_ground, id_guess);
		mse += norm([
			poses_guess(i).x - pose_gt.x;
			poses_guess(i).y - pose_gt.y
		]);
	end
	mse /= n_poses;
end
