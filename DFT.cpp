#include "DFT.h"
#include <iostream>

double DFT::VWN_eps_c(double r_s) {
	
	//Vosko Wilk Nusair paper, with VWN5 constants
	double atanQ_rb = std::atan(VWN_Q / (2.0 * r_s + VWN_b));

	//term 3
	double eps_c = std::log((r_s - VWN_x0) * (r_s - VWN_x0) / VWN_X(r_s));
	eps_c += atanQ_rb * (2.0 * (VWN_b + 2.0 * VWN_x0) / VWN_Q);
	//term 3 coeff
	eps_c *= -VWN_b * VWN_x0 / VWN_X_x0;

	//add in term 1 and term 2
	eps_c += std::log(r_s * r_s / VWN_X(r_s));
	eps_c += (2.0 * VWN_b / VWN_Q) * atanQ_rb;

	//whole expression coeff = (1/(2pi^2))
	eps_c *= VWN_A;

	return eps_c;
}

void DFT::calculateV_xc(std::vector<std::vector<std::vector<double>>>& n, std::vector<std::vector<std::vector<double>>>& V_xc) {
	//Local Density Approximation:
	
	//V_x(r) = ex_coeff * n(r)^(1/3) for exchange functional
	double ex_coeff = -pow(3.0 / std::numbers::pi, 1.0 / 3.0);

	//convert n(r) to Wigner-Seitz radius for correlation functional
		// via r_s = r_s_coeff * n^(-1/3)
	double r_s_coeff = pow(3.0 / (4.0 * std::numbers::pi), 1.0 / 3.0);

	//Iterate over n(r)
	for (int x = 0;x < n.size();x++) {
		for (int y = 0;y < n[x].size();y++) {
			for (int z = 0;z < n[x][y].size();z++) {
				//beware divide by 0
				if (std::abs(n[x][y][z]) < 1e-10) {
					V_xc[x][y][z] = 0.0;
					continue;
				}

				//Start by setting V_xc(r) = V_x(r)
				V_xc[x][y][z] = ex_coeff * pow(n[x][y][z], 1.0 / 3.0);

				//Wigner-Seitz radius for VWN correlation functional
				double r_s = r_s_coeff * pow(1.0 / n[x][y][z], 1.0 / 3.0);

				//Calculate correlation potental and derivative according to VWN LDA
				double eps_c = VWN_eps_c(r_s);
				double delta_r = r_s * 1e-7;
				double d_eps_c = (VWN_eps_c(r_s + delta_r) - eps_c) / delta_r;
				
				//Increment V_xc by V_c
				V_xc[x][y][z] += eps_c - (r_s * d_eps_c / 3.0);
				 
			}
		}
	}	
}

void DFT::calculateV_ht(fftw_plan& FFT_plan, fftw_plan& FFT_inv_plan, std::vector<std::complex<double>>& FFT_out, 
	std::vector<double>& n_flat, std::vector<std::vector<std::vector<double>>>& V_ht, double cellSize) {
	//TODO
	//Poisson equation, use LinearOperator::Laplacian<double>
	//May well need to reference a different solver? In case CG is not ideal here??
	
	//FFT on n(r) --> frequency space
	fftw_execute(FFT_plan);
	std::cout << "\nfirst fft done\n\n";

	std::cout << "FFT[1][0][0]: " << FFT_out[1 * V_ht[0].size() * (V_ht[0][0].size() / 2 + 1)] << "\n";
	std::cout << "FFT[0][1][0]: " << FFT_out[1 * (V_ht[0][0].size() / 2 + 1)] << "\n";
	std::cout << "FFT[0][0][1]: " << FFT_out[1] << "\n\n";
	std::cout << "phase: " << std::arg(FFT_out[1 * V_ht[0].size() * (V_ht[0][0].size() / 2 + 1)]) << "\n";

	//Calculate Lx Ly Lz = bohr volume of whole grid
	double Lx = ((double)V_ht.size()) * cellSize;
	double Ly = ((double)V_ht[0].size()) * cellSize;
	double Lz = ((double)V_ht[0][0].size()) * cellSize;

	const double tau = 2.0 * std::numbers::pi;

	// solve for the complex plane wave inputs to my reverse FFT, output is V_ht.
	for (int i_x = 0; i_x < V_ht.size(); i_x++) {
		//
		for (int i_y = 0; i_y < V_ht[0].size(); i_y++) {

			for (int i_z = 0; i_z < (V_ht[0][0].size() / 2 + 1); i_z++) {
				
				//grab n_tilde(G) in flat array
				std::complex<double>& n_tilde = FFT_out[
					  i_x * (V_ht[0][0].size()/2+1) * V_ht[0].size()
					+ i_y * (V_ht[0][0].size()/2+1)
					+ i_z
				];

		
				//get |G|^2 , Laplacian in frequency space, ternary op corrects for FFT negative wrapping
				//not putting the int casts before V_ht.size() gave me weekslong headache. Unsigned overflow ruined all negative frequencies.
				// a lifelong caution around unsigned ints was built here. 
				double Gx = tau / Lx * (double)( (i_x <= (int)V_ht.size()    / 2) ? i_x : (i_x - (int)V_ht.size()    ) );
				double Gy = tau / Ly * (double)( (i_y <= (int)V_ht[0].size() / 2) ? i_y : (i_y - (int)V_ht[0].size() ) );
				double Gz = tau / Lz * (double)   i_z;
				double sqMagG = (Gx * Gx) + (Gy * Gy) + (Gz * Gz);
				if (i_x == 1 && i_y == 0 && i_z == 0) std::cout << "G2(1,0,0): " << sqMagG << "\n";
				if (i_x == 0 && i_y == 1 && i_z == 0) std::cout << "G2(0,1,0): " << sqMagG << "\n";
				if (i_x == 0 && i_y == 0 && i_z == 1) std::cout << "G2(0,0,1): " << sqMagG << "\n";
				if (i_x == 1 && i_y == 1 && i_z == 1) std::cout << "G2(1,1,1): " << sqMagG << "\n";
				//no div by 0
				int zeroed = 0;
				if (sqMagG <= 1e-10) {
					n_tilde = std::complex<double>(0.0,0.0);
					zeroed++;
					std::cout << "zerod ";
					continue;
				}

				//solve poisson equation in frequency space -- ooh so clean :)
				n_tilde = std::complex<double>(n_tilde.real() * -2.0 * tau / sqMagG, n_tilde.imag() * -2.0 * tau / sqMagG);



			}
		}
		
	}
	
	//reverse transform (output in n_flat)

	std::cout << "FFT[1][0][0]: " << FFT_out[1 * V_ht[0].size() * (V_ht[0][0].size() / 2 + 1)] << "\n";
	std::cout << "FFT[0][1][0]: " << FFT_out[1 * (V_ht[0][0].size() / 2 + 1)] << "\n";
	std::cout << "FFT[0][0][1]: " << FFT_out[1] << "\n\n";


	fftw_execute(FFT_inv_plan);
	std::cout << "second fft done\n";

	//normalize
	double norm = V_ht.size()* V_ht[0].size()* V_ht[0][0].size();

	//set V_ht as result of reverse ftt ( from flattened, reusing n_flat as output )
	for (int x = 0;x < V_ht.size();x++) {
		for (int y = 0;y < V_ht[0].size();y++) {
			for (int z = 0;z < V_ht[0][0].size();z++) {
				
			
				V_ht[x][y][z] = (n_flat[
					x * V_ht[0][0].size() * V_ht[0].size()
						+ y * V_ht[0][0].size()
						+ z
				]
				* (1.0/norm) );
				
			}
		}

	}

}

void DFT::calculateV_s(Fields& fields) {

	//V_s = V_ext + V_ht + V_xc
		//V_ext is set beforehand
		//V_xc is simple and from n(r)
		//V_ht is hard and from poisson/FFT problem

	//Start with V_xc
	
	DFT::calculateV_xc(fields.n, fields.V_xc);
	
	//Flatten n for FFT
	Fields::FlattenField(fields.n, fields.n_flat);

	//Then solve V_ht
	std::cout << "starting HT";
	
	
	DFT::calculateV_ht(fields.FFT_plan_n, fields.FFT_inv_plan_n, fields.FFT_out, fields.n_flat, fields.V_hartree, fields.cellSize);
	std::cout << "finsihed HT";
	//Fields::FlattenN(fields.V_hartree, fields.n_flat);


	//std::cout << "n center: " << fields.V_ext[128][128][128] << "\n";
	//std::cout << "V center: " << fields.V_hartree[128][128][128] << "\n";
	//std::cout << "ratio: " << fields.V_hartree[128][128][128] / fields.V_ext[128][128][128] << "\n";
	//std::cout << "Nx*Ny*Nz: " << 257 * 257 * 257 << "\n";
	//DFT::calculateV_ht(fields.FFT_plan, fields.FFT_inv_plan, fields.FFT_out, fields.n_flat, fields.V_xc, fields.cellSize);

	//Set V_s by element-wise addition
	for (int x = 0;x < fields.n.size();x++) {
		for (int y = 0;y < fields.n[x].size();y++) {
			for (int z = 0;z < fields.n[x][y].size();z++) {
				fields.V_s[x][y][z] = fields.V_ext[x][y][z] + fields.V_hartree[x][y][z] + fields.V_xc[x][y][z];
			}
		}
	}

}

//TODO:


void DFT::Hamiltonian::perform_op(const Scalar* x_in, Scalar* y_out) const{
	//TODO
	//needs to preform laplaican and multiply Vs into input vector (arbitrary entire flattened grid)
		//this being in frequency space, means i take the inverse fft to real space to multiply vs, then fft back and add to laplacian
	
	//Grab fields from static var
	Fields* fields = getFields();


	//okay you have everything you need to go forwards..

	//x_in is the entire orbital field, flattened as with FFTW format 
	// (but is this time in frequency space, symetric along z, so half as many z idxs)
	//y_out is flattened result vector
	// (this is same size as x_in [q.v 'operator'], so also is in frequency space.
	

	
	//laplacian 
	//Calculate Lx Ly Lz = bohr volume of whole grid
	int Nx = fields->V_s.size();
	int Ny = fields->V_s[0].size();
	int Nz = fields->V_s[0][0].size();

	double Lx = ((double)Nx) * fields->cellSize;
	double Ly = ((double)Ny) * fields->cellSize;
	double Lz = ((double)Nz) * fields->cellSize;

	const double tau = 2.0 * std::numbers::pi;

	std::cout << "len =" << fields->fft_orbital.size();

	for (int i_x = 0; i_x < Nx; i_x++)
	{
		for (int i_y = 0; i_y < Ny; i_y++)
		{
			for (int i_z = 0; i_z < Nz/2+1; i_z++)
			{
				double Gx = tau / Lx * (double)((i_x <= (int)Nx / 2) ? i_x : (i_x - (int)Nx));
				double Gy = tau / Ly * (double)((i_y <= (int)Ny/ 2) ? i_y : (i_y - (int)Ny));
				double Gz = tau / Lz * (double)i_z;

				double sqMagG = (Gx * Gx) + (Gy * Gy) + (Gz * Gz);
				int idx = 
					i_x * (Nz / 2 + 1) * Ny
					+ i_y * (Nz / 2 + 1)
					+ i_z;

				//this sets y = lacplacian(x)
					//next step will add the result of (x*Vs) to y
				y_out[idx] = x_in[idx] * 0.5 * sqMagG;
				fields->fft_orbital[idx] = x_in[idx];
			}
		}
	}

	//N * Vs, component-wise in real space

	//step 1, ifft x to real space.
	fftw_execute(fields->FFT_inv_plan_KohnSham); //results in fields->ifft_orbital[] (real valued full)


	//multiply Vs
	int idx_ = 0;
	for (int x = 0;x < Nx;x++) {
		for (int y = 0;y < Ny;y++) {
			for (int z = 0;z < Nz;z++) {

				fields->ifft_orbital[idx_] *= fields->V_s[x][y][z];
					idx_++;
			}
		}
		
	}

	//fft back to FreqSpace
	fftw_execute(fields->FFT_plan_KohnSham);


	double norm = Nx*Ny*Nz;

	//Add To laplacian term to get final y_out
	for (int i = 0;i < Nx * Ny * (Nz/2+1);i++) {
		y_out[i] += fields->fft_orbital[i]/norm;
	}

}


void DFT::calculateOrbitals(Fields& fields, fftw_plan& FFTplan, fftw_plan& IFFTplan) {
	
	//orbital plane waves in fourier space:
	
	


	//sample code from docs:
	DFT::Hamiltonian op;
	op.setFields(&fields);
	Spectra::HermEigsSolver<DFT::Hamiltonian> eigs(op, fields.nFreeElectrons, fields.nFreeElectrons*2);//need to figure out how many eigenvalues i need? each one is an electron, kinda?? mybe??
	eigs.init();


	//run eigensolver:
	int nconv = eigs.compute(Spectra::SortRule::SmallestAlge, 100, 1e-4);
	std::cout << nconv << " eigs \n";

	Eigen::MatrixXcd evecs = eigs.eigenvectors();
	//std::cout << "Eigenvectors:\n" << evecs << std::endl;
	
	for (int orb = 0; orb < nconv;orb++) {
		for (int idx = 0;idx < fields.orbitals[0].size(); idx++) {
			fields.orbitals[orb][idx] = evecs(idx, orb);
		}
	}

}



void DFT::calculateN(Fields& fields) {
	//TODO,
	for (int orb = 0; orb < fields.nFreeElectrons; orb++) {
		fields.fft_orbital = fields.orbitals[orb];
		fftw_execute(fields.FFT_inv_plan_KohnSham);

		//normalize and rescale orbital so wave function integrates to 1 over physical volume
		double norm = 0;
		for (int i = 0;i < fields.ifft_orbital.size();i++) {
			double mag = abs(fields.ifft_orbital[i]);
			norm += mag * mag;
		}
		for (int i = 0;i < fields.ifft_orbital.size();i++) {
			fields.ifft_orbital[i] /= sqrt(norm);
			fields.ifft_orbital[i] /= sqrt(fields.cellSize * fields.cellSize * fields.cellSize);
		}


		for (int x = 0;x < fields.n.size();x++) {
			for (int y = 0;y < fields.n[0].size();y++) {
				for (int z = 0;z < fields.n[0][0].size();z++) {
					if (orb==0) fields.n[x][y][z] *=0.8;
					double mag = abs(fields.ifft_orbital[x * fields.n[0][0].size() * fields.n[0].size() + y * fields.n[0][0].size() + z]);
					fields.n[x][y][z] += mag * mag * 0.05;
				
				}
			}
		}
	};
	
}
