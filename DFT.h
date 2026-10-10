#pragma once
#include <vector>
#include <complex>
#include <numbers>
#include "LinearOperator.h"
#include "Fields.h"
#include "fftw3.h"
#include <Eigen/Core>
#include <Spectra/HermEigsSolver.h>
#include <omp.h>



/// Takes references from Fields class, and directly updates each respective field.
class DFT {
public:
	inline static int NumGridPoints = 0;

	static void calculateV_s(Fields& fields);

	//TOD)
	//needs references to respective objects
	static void calculateOrbitals(Fields& fields, fftw_plan& FFTplan, fftw_plan& IFFTplan);
	static void calculateN(Fields& fields);

	class Hamiltonian
	{
	public:
		using Scalar = std::complex<double>;   // A typedef named "Scalar" is required
		int rows() const { return DFT::NumGridPoints; }
		int cols() const { return DFT::NumGridPoints; }
		// y_out = M * x_in
		void perform_op(const Scalar* x_in, Scalar* y_out) const;
		static void setFields(Fields* F) { 
			fields = F;
			NumGridPoints = fields->w * fields->h * (fields->l / 2 + 1);
		}
	private:
		inline static Fields* fields;
		static Fields* getFields() { return fields; }
	};
	

private:
	
	static void calculateV_xc(std::vector<std::vector<std::vector<double>>>& n, std::vector<std::vector<std::vector<double>>>& V_xc);
	static void calculateV_ht(fftw_plan& FFT_plan, fftw_plan& FFT_inv_plan, std::vector<std::complex<double>>& FFT_out, 
		std::vector<double>& n_flat, std::vector<std::vector<std::vector<double>>>& V_ht, double cellSize);

	//helper vars and methods:
	static constexpr double VWN_A	 = 0.0310907;
	static constexpr double VWN_b	 = 3.72744;
	static constexpr double VWN_c	 = 12.9352;
	static constexpr double VWN_Q	 = 6.1520;		// = sqrt(4c-b^2)
	static constexpr double VWN_x0	 = -0.10498;
	static constexpr double VWN_X_x0 = 12.5549; // = VWN_X(VWN_x0)

	static double VWN_X(double x) {
		return x*x + VWN_b * x + VWN_c;
	}
	static double VWN_eps_c(double r_s);

	/*
	static void PrintHamiltonianStats(Fields* fields, std::vector<std::complex<double>> x_in)  {
		//TODO
		//needs to preform laplaican and multiply Vs into input vector (arbitrary entire flattened grid)
			//this being in frequency space, means i take the inverse fft to real space to multiply vs, then fft back and add to laplacian

		//Grab fields from static var
		
		//std::cout << fields->V_s[32][32][32] << "\n";

		//okay you have everything you need to go forwards..

		//x_in is the entire orbital field, flattened as with FFTW format 
		// (but is this time in frequency space, symetric along z, so half as many z idxs)
		//y_out is flattened result vector
		// (this is same size as x_in [q.v 'operator'], so also is in frequency space.
		std::vector<std::complex<double>> y_out = std::vector<std::complex<double>>(x_in.size());


		//laplacian 
		//Calculate Lx Ly Lz = bohr volume of whole grid
		int Nx = fields->V_s.size();
		int Ny = fields->V_s[0].size();
		int Nz = fields->V_s[0][0].size();

		double Lx = ((double)Nx) * fields->cellSize;
		double Ly = ((double)Ny) * fields->cellSize;
		double Lz = ((double)Nz) * fields->cellSize;

		const double tau = 2.0 * std::numbers::pi;

		std::cout << "  len =" << fields->fft_orbital.size();
		std::complex<double> Kinetic_DotProduct = 0.0;
		#pragma omp parallel for collapse(3)
		for (int i_x = 0; i_x < Nx; i_x++)
		{
			for (int i_y = 0; i_y < Ny; i_y++)
			{

				for (int i_z = 0; i_z < Nz/2+1; i_z++)
				{
					double Gx = tau / Lx * (double)((i_x <= (int)Nx / 2) ? i_x : (i_x - (int)Nx));
					double Gy = tau / Ly * (double)((i_y <= (int)Ny / 2) ? i_y : (i_y - (int)Ny));
					double Gz = tau / Lz * (double)(i_z);

					double sqMagG = (Gx * Gx) + (Gy * Gy) + (Gz * Gz);
					int idx =
						i_x * (Nz/2+1)*Ny
						+ i_y * (Nz/2+1)
						+i_z;


					double ecut = 30.0; // Hartree
					if (0.5 * sqMagG > ecut) {
						y_out[idx] = 0;
						fields->fft_orbital[idx] = 0;
						continue;
					}
					//this sets y = lacplacian(x)
						//next step will add the result of (x*Vs) to y
					y_out[idx] = (x_in[idx] * 0.5 * sqMagG) * x_in[idx];
					fields->fft_orbital[idx] = x_in[idx];
					Kinetic_DotProduct += y_out[idx];
					
				}
			}
		}
		std::cout << "Kinetic: " << Kinetic_DotProduct << "\n";


		//N * Vs, component-wise in real space

		//step 1, ifft x to real space.
		fftw_execute(fields->FFT_inv_plan_KohnSham); //results in fields->ifft_orbital[] (real valued full)
		
		double orbitalnorm = 0.0;
		for (int i = 0;i < fields->ifft_orbital.size();i++) {
			orbitalnorm += fields->ifft_orbital[i].real() * fields->ifft_orbital[i].real() + fields->ifft_orbital[i].imag() * fields->ifft_orbital[i].imag();

		}
		for (int i = 0;i < x_in.size();i++) {
			fields->ifft_orbital[i] /= sqrt(orbitalnorm);
			fields->ifft_orbital[i] /= sqrt(fields->cellSize * fields->cellSize * fields->cellSize);
		}

		//multiply Vs
		std::complex<double> V_ext_contrib = 0.0;
		std::complex<double> V_ht_contrib = 0.0;
		std::complex<double> V_xc_contrib = 0.0;
		#pragma omp parallel for collapse(3)
		for (int x = 0;x < Nx;x++) {
			for (int y = 0;y < Ny;y++) {
				for (int z = 0;z < Nz;z++) {
					int idx =
						x * (Nz)*Ny
						+ y * (Nz)
						+z;
					y_out[idx] = fields->ifft_orbital[idx] * fields->V_ext[x][y][z];
					V_ext_contrib += y_out[idx] * fields->ifft_orbital[idx]* fields->cellSize * fields->cellSize * fields->cellSize;

					y_out[idx] = fields->ifft_orbital[idx] * fields->V_hartree[x][y][z];
					V_ht_contrib += y_out[idx] * fields->ifft_orbital[idx]* fields->cellSize * fields->cellSize * fields->cellSize;

					y_out[idx] = fields->ifft_orbital [idx] * fields->V_xc[x][y][z];
					V_xc_contrib += y_out[idx] * fields->ifft_orbital[idx]* fields->cellSize * fields->cellSize * fields->cellSize;
				}
			}

		}
		std::cout << "Vext: " << V_ext_contrib << "\n";
		std::cout << "Vht: " << V_ht_contrib << "\n";
		std::cout << "Vxc: " << V_xc_contrib << "\n";
		//fft back to FreqSpace
		fftw_execute(fields->FFT_plan_KohnSham);


		double norm = Nx * Ny * Nz;


		//Add To laplacian term to get final y_out
		#pragma omp parallel for
		for (int i = 0;i < Nx * Ny * (Nz);i++) {
			y_out[i] += fields->fft_orbital[i] / norm;
		}

		#pragma omp parallel for collapse(3)
		for (int i_x = 0; i_x < Nx; i_x++)
		{
			for (int i_y = 0; i_y < Ny; i_y++)
			{
				for (int i_z = 0; i_z < Nz; i_z++)
				{
					double Gx = tau / Lx * (double)((i_x <= (int)Nx / 2) ? i_x : (i_x - (int)Nx));
					double Gy = tau / Ly * (double)((i_y <= (int)Ny / 2) ? i_y : (i_y - (int)Ny));
					double Gz = tau / Lz * (double)((i_z <= (int)Nz / 2) ? i_z : (i_z - (int)Nz));

					double sqMagG = (Gx * Gx) + (Gy * Gy) + (Gz * Gz);
					int idx =
						i_x * (Nz)*Ny
						+ i_y * (Nz)
						+i_z;

					double ecut = 30.0; // Hartree
					if (0.5 * sqMagG > ecut) {
						y_out[idx] = 0;
					}
				}
			}
		}

	}
	*/
	
	
};