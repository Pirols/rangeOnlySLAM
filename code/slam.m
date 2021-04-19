close all
clear
clc

addpath 'tools/g2o_wrapper'
addpath 'tools/visualization'
source "tools/utilities/geometry_helpers_2d.m"

%load your own dataset dataset
[_, poses, transitions, observations] = loadG2o('../03-RangeOnlySLAM/slam2d_range_only_initial_guess.g2o');

% remove first elements of structs which are empty and meaningless
poses = poses(2:end);
transitions = transitions(2:end);
observations = observations(2:end);

%% init stuff

% initial pose
mu = [0; 0; 0];

% initialize landmarks using odometry edges
landmarks = initLandmarks(poses, observations);

% init trajectory
trajectory = [mu(1), mu(2)];

% least squares optimization

% render graphics
% figure(1); title("rangeonly-slam");
% plotState(landmarks, mu);
