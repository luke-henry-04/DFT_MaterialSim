#pragma once
#include <vector>
#include <complex>
#include <algorithm>



class LinearOperator {
public:
	template <typename T>
	static std::vector<std::vector<std::vector<T>>> Laplacian(std::vector<std::vector<std::vector<T>>>& input) {
		//kernel convoltion, dont make entire matrix, implement sum directly at each point
		//TODO
		std::vector<std::vector<std::vector<T>>> result = input;

		for (int x = 0;x < input.size();x++) {
			for (int y = 0;y < input[0].size();y++) {
				for (int z = 0;z < input[0][0].size();z++) {

					int xl = std::min(std::max(x - 1, 0), (int)input.size() - 1);
					int xr = std::max(0, std::min(x + 1, (int)input.size() - 1));
					int yl = std::min(std::max(y - 1, 0), (int)input[0].size() - 1);
					int yr = std::max(0, std::min(y + 1, (int)input[0].size() - 1));
					int zl = std::min(std::max(z - 1, 0), (int)input[0][0].size() - 1);
					int zr = std::max(0, std::min(z + 1, (int)input[0][0].size() - 1));

					result[x][y][z] *= -6.0;
					result[x][y][z] += input[xl][y][z] + input[xr][y][z]
						+ input[x][yl][z] + input[x][yr][z]
						+ input[x][y][zl] + input[z][y][zr]						
					;


				}
			}
		}

		return result;
	}

	static std::vector<std::vector<std::vector<std::complex<double>>>> Hamiltonian(std::vector<std::vector<std::vector<std::complex<double>>>>& input, std::vector<std::vector<std::vector<double>>>& V_s) {
		//will need to call laplacian, and also add in component-wise multiplication by V_s. 
		//TODO

		return std::vector<std::vector<std::vector<std::complex<double>>>>();
	}

};