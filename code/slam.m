close all
clear
clc

addpath 'tools/g2o_wrapper'
addpath 'tools/visualization'
source "tools/utilities/geometry_helpers_2d.m"

%load your own dataset dataset
[landmarks, poses, transitions, observations] = loadG2o('../03-RangeOnlySLAM/slam2d_range_only_initial_guess.g2o');

% Check struct fields
% landmarks
% poses
% transitions
% observations
% observations(2).observation(1)

% just adding a landmark for demonstration purposes
landmarks(end) = landmark(1, [0, 2])

%% init stuff
%initial pose
mu = rand(3,1)*20-10;
mu(3) = normalizeAngle(mu(3));
printf('Random initial pose: [%f, %f, %f]\n', mu(1),mu(2), mu(3));
fflush(stdout);

%init covariance
sigma = eye(3)*0.001;

%init graphics
figure(1); title("rangeonly-slam");
plotState(landmarks, mu);
