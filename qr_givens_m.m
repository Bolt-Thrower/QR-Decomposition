% -------------------------------------------------------------------------
%  qr_givens_m.m - MATLAB implementation of QR Decomposition using Givens
%  rotations

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

function R = qr_givens_m(A, J)

% Initialize matrices
R = A;
tol = 0.00001;

% Start rotations
r = 0;  % Column index
c = 0;  % Row index

% Size of A
[m, n] = size(A);

% We only care about real matrices
while c < min(n, m)
    % Increment row and column index
    c = c + 1;
    r = r + 1;
    
    % Next row pointer. Initialize with the current pointer now
    r1 = r;
    r2 = r1;
    
    while r2 < m
        % Increment r2
        r2 = r2 + 1;
        
        % Get the elements
        a1 = R(r1, c);
        a2 = R(r2, c);
        
        if abs(a2) <= tol * abs(a1)
            % Already in the required structure
            R(r2, c) = 0;
        elseif abs(a1) <= tol * abs(a2)
            % Reverse structure, swap rows
            s = R(r1, :);
            R(r1, :) = R(r2, :);
            R(r2, :) = s;
            R(r2, c) = 0;
        elseif abs(a1) > abs(a2)
            t = a2 / a1;
            z = sqrt(1 + t^2);
            R(r1, :) = (R(r1, :) + t * R(r2, :)) / z;
            R(r2, :) = -t * R(r1, :) + z * R(r2, :);
            R(r2, c) = 0;
        else
            t = a1 / a2;
            z = sqrt(1 + t^2);
            R(r1, :) = (t * R(r1, :) + R(r2, :)) / z;
            R(r2, :) = (-R(r1, :) + z * R(r2, :)) / t;
            R(r2, c) = 0;
        end           
    end
end