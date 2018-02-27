//-------------------------------------------------------------------------
// qr_givens.cc - Armadillo implementation of QR decomposition using Givens
// rotations
//
// Copyright 2018 Seedo Eldho Paul <seedoeldhopaul@gmail.com>
// 
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
// http://www.apache.org/licenses/LICENSE-2.0
// 
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.
// ------------------------------------------------------------------------

#include "armaMex.hpp"


static const double TOL = 0.000001;


void qr_givens(mat& R)
{
	// Get size of the matrix
	int m = R.n_rows;
	int n = R.n_cols;

	int S = std::min(m, n);

	// Start the rotations
	for (int r = 0, c = 0; c < S - 1; r++, c++) {
		// Iterate over the rows
		for (int r2 = r + 1; r2 < m; r2++) {
			// Get the elements
			// Get first element in the row
			double a1 = R(r, c);
			double a2 = R(r2, c);
			double a1_abs = std::fabs(a1);
			double a2_abs = std::fabs(a2);

			if (a2_abs <= TOL * a1_abs) {
				// Good, already in the required structure
				R(r2, c) = 0;
			} else if (a1_abs <= TOL * a2_abs) {
				// Swap, then we have the required structure
				R.swap_rows(r, r2);
				R(r2, c) = 0;
			} else if (a1_abs > a2_abs) {
				// Mixed downdating
				double t = a2 / a1;
				double z = std::sqrt(1 + t * t);
				R.row(r) = (R.row(r) + t * R.row(r2)) / z;
				R.row(r2) = -t * R.row(r) + z * R.row(r2);
				R(r2, c) = 0;
			} else {
				// Mixed downdating
				double t = a1 / a2;
				double z = std::sqrt(1 + t * t);
				R.row(r) = (t * R.row(r) + R.row(r2)) / z;
				R.row(r2) = (-R.row(r) + z * R.row(r2)) / t;
				R(r2, c) = 0;
			}
		}
	}
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	mat A = armaGetPr(prhs[0]);
	mat B = A;

	// Perform QR
	qr_givens(B);

	// Create output arguments and return
	plhs[0] = armaCreateMxMatrix(B.n_rows, B.n_cols);
	armaSetPr(plhs[0], B);
}