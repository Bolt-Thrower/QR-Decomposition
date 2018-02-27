//-------------------------------------------------------------------------
// qr_house.cc - Armadillo implementation of QR decomposition using
// Householder reflections
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

using namespace arma;

// Housholder QR, own implementation
void qr_house(mat& A)
{
	int m = A.n_rows;
	int n = A.n_cols;

	for (int i = 0; i < m - 1; i++) {
		// Get the relevant column
		vec x = A(span(i, n-1), i);

		// Find the vector `v`.
		// For convenience, we use the vector x
		x(0) -= norm(x);

		// Update A
		A(span(i, n-1), span(i, n-1)) -= 2 / as_scalar(x.t() * x) * x * (x.t() * A(span(i, n-1), span(i, n-1)));
	}
}


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	mat A = armaGetPr(prhs[0]);
	mat B = A;

	// Perform QR
	qr_house(B);

	// Create output arguments and return
	plhs[0] = armaCreateMxMatrix(B.n_rows, B.n_cols);
	armaSetPr(plhs[0], B);
}