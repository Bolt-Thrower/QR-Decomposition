% -------------------------------------------------------------------------
%  benchmark.m - Benchmarking different QR decomposition algorithms

%  Copyright 2018 Seedo Eldho Paul <seedoeldhopaul@gmail.com>
 
%  Licensed under the Apache License, Version 2.0 (the "License");
%  you may not use this file except in compliance with the License.
%  You may obtain a copy of the License at
%  http://www.apache.org/licenses/LICENSE-2.0
 
%  Unless required by applicable law or agreed to in writing, software
%  distributed under the License is distributed on an "AS IS" BASIS,
%  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
%  See the License for the specific language governing permissions and
%  limitations under the License.
%  ------------------------------------------------------------------------

clc;
close all;

% Compile all the required source files
compile;

N = 100;                    % Maximum size of the matrix
M = 200;                    % Number of repetitions
time = zeros(N, 4);         % Array to store the times for each size
inst_time = zeros(M, 1);    % For storing instantaneous times

% First, built in qr()
for i = 1:N
    % Test M times
    for j = 1:M
        % Create a matrix
        A = randn(i, i);
        % Initialize the clock
        tic;
        % Do the decomposition
        %---------------------
        [~, R] = qr(A);
        %---------------------
        % Calculate time elapsed
        inst_time(j) = toc;
    end
    % Average the times and store in the array
    time(i, 1) = mean(inst_time);
end

% Now, householder QR (Armadillo)
for i = 1:N
    % Test M times
    for j = 1:M
        % Create a matrix
        A = randn(i, i);
        % Initialize the clock
        tic;
        % Do the decomposition
        %---------------------
        R = qr_house(A);
        %---------------------
        % Calculate time elapsed
        inst_time(j) = toc;
    end
    % Average the times and store in the array
    time(i, 2) = mean(inst_time);
end

% Givens QR, matlab
for i = 1:N
    % Test M times
    for j = 1:M
        % Create a matrix
        A = randn(i, i);
        % Initialize the clock
        tic;
        % Do the decomposition
        %---------------------
        R = qr_givens_m(A);
        %---------------------
        % Calculate time elapsed
        inst_time(j) = toc;
    end
    % Average the times and store in the array
    time(i, 3) = mean(inst_time);
end

% Givens QR, arma
for i = 1:N
    % Test M times
    for j = 1:M
        % Create a matrix
        A = randn(i, i);
        % Initialize the clock
        tic;
        % Do the decomposition
        %---------------------
        R = qr_givens(A);
        %---------------------
        % Calculate time elapsed
        inst_time(j) = toc;
    end
    % Average the times and store in the array
    time(i, 4) = mean(inst_time);
end

figure('Name', 'Comparison of QR Decomposition');
plot(1:N, time(:, 1));
xlabel('Size of the matrix'), ylabel('Elapsed time');
title('Size of the matrix vs time consumed');
hold on;

plot(1:N, time(:, 2));
plot(1:N, time(:, 3));
plot(1:N, time(:, 4));
legend('MATLAB QR', 'Householder QR (Arma)', 'Givens MATLAB', 'Givens ARMA');
