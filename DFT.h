#pragma once
#include <vector>
#include <complex>
#include <numbers>
#include "LinearOperator.h"
#include "Fields.h"
#include "fftw3.h"
#include <Eigen/Core>
#include <Spectra/HermEigsSolver.h>



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
		static void setFields(Fields* F) { fields = F; }
	private:
		inline static Fields* fields;
		static Fields* getFields() { return fields; }
	};
	

private:
	
	static void calculateV_xc(std::vector<std::vector<std::vector<double>>>& n, std::vector<std::vector<std::vector<double>>>& V_xc);
	static void calculateV_ht(fftw_plan& FFT_plan, fftw_plan& FFT_inv_plan, std::vector<std::complex<double>>& FFT_out, 
		std::vector<double>& n_flat, std::vector<std::vector<std::vector<double>>>& V_ht, double cellSize);

	//helper vars and methods:
	static constexpr double VWN_A	 = 0.5 / (std::numbers::pi * std::numbers::pi);
	static constexpr double VWN_b	 = 3.72744;
	static constexpr double VWN_c	 = 12.9352;
	static constexpr double VWN_Q	 = 6.1520;		// = sqrt(4c-b^2)
	static constexpr double VWN_x0	 = -0.10498;
	static constexpr double VWN_X_x0 = 12.5549; // = VWN_X(VWN_x0)

	static double VWN_X(double r_s) {
		return r_s * r_s + VWN_b * r_s + VWN_c;
	}
	static double VWN_eps_c(double r_s);


	
	
};