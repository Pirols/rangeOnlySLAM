close all
clear
clc

addpath 'g2o_wrapper'
addpath 'visualization'
source 'mse.m'

% To disable singularity warnings (couldn't find a way to disable just that kind of warnings)
% warning("off", "all")

%%% load datasets
% initial guess
[~, poses, transitions, observations] = loadG2o('../datasets/slam2d_range_only_initial_guess.g2o');
% ground truth
[landmarks_gt, poses_gt, transitions_gt, observations_gt] = loadG2o('../datasets/slam2d_range_only_ground_truth.g2o');

%%% remove first elements of structs which are empty and meaningless
%%% plus compute initial guess of landmarks (trilatering relevant observations)
% initial guess
poses(1) = [];
transitions(1) = [];
observations(1) = [];
landmarks = initLandmarks(poses, observations);
% ground truth
landmarks_gt(1) = [];
poses_gt(1) = [];
transitions_gt(1) = [];
observations_gt(1) = [];

%%% least squares optimization
[poses_ls, landmarks_ls] = slamLS(
    landmarks,
    poses(1),
    transitions,
    observations,
    10     % algorithm number of iterations
);

%%% printing some evaluation metrics
mse_lands = mse_land(landmarks_gt, landmarks);
mse_lands_ls = mse_land(landmarks_gt, landmarks_ls);
printf("MSE over landmarks went from %d to %d (after LS optimization)\n",
        mse_lands, mse_lands_ls);
mse_poses_xy = mse_pose(poses_gt, poses);
mse_poses_xy_ls = mse_pose(poses_gt, poses_ls);
printf("MSE over poses(only xy) went from %d to %d (after LS optimization)\n",
        mse_poses_xy, mse_poses_xy_ls);

% % render graphics
% figure(1); title("rangeonly-slam");
% plotState(landmarks, mu);
