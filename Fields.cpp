#include "Fields.h"


void Fields::InitMaterial(std::string filename) {
	//Initialize material properties from custom (standardized??) file
	//TODO


}

void Fields::InitMaterial(int w_, int h_, int l_, double cellSize_bohr, int atoms, int protons, double charge) {
	w = w_;
	h = h_;
	l = l_;
	cellSize = cellSize_bohr;

	//Initialize grids to  0.0f  and size  w x h x l 
	n = std::vector<std::vector<std::vector<double>>>(w,
		std::vector<std::vector<double>>(h,
			std::vector<double>(l,
				((double)(atoms * protons)) / (((double)(w * l * h)) * cellSize_bohr * cellSize_bohr * cellSize_bohr)
		)));

	//n[0][0][0] = 1.0;

	std::cout << ((double)(atoms * protons)) / (((double)(w * l * h)) * cellSize_bohr * cellSize_bohr * cellSize_bohr);
	V_ext = n;
	V_hartree = n;
	V_xc = n;
	V_s = n;

	//Generate flat N and FFTPlan, set up output array
	n_flat = std::vector<double>(w*h*l, 0.0);
	Fields::FlattenField(n, n_flat);
	FFT_out.resize(w * h * (l / 2 + 1));

	//3D FFT on n(r) --> convert to freequency space
	FFT_plan_n = fftw_plan_dft_r2c_3d(w, h,l, n_flat.data(), reinterpret_cast<fftw_complex*>(FFT_out.data()), FFTW_ESTIMATE);
	FFT_inv_plan_n = fftw_plan_dft_c2r_3d(w, h,l, reinterpret_cast<fftw_complex*>(FFT_out.data()), n_flat.data(), FFTW_ESTIMATE);

	//orbital vectors and FFT plan for solving Kohn Sham
	ifft_orbital	= std::vector<double>(n_flat.size(),0);
	fft_orbital		= std::vector<std::complex<double>>(w*h*((l/2)+1), 0);
	FFT_plan_KohnSham		= fftw_plan_dft_r2c_3d(w,h,l, ifft_orbital.data(), reinterpret_cast<fftw_complex*>(fft_orbital.data()), FFTW_ESTIMATE);
	FFT_inv_plan_KohnSham	= fftw_plan_dft_c2r_3d(w,h,l,reinterpret_cast<fftw_complex*>(fft_orbital.data()), ifft_orbital.data(), FFTW_ESTIMATE);

	
	//Create nuclei
	nFreeElectrons = 0;
	for (int i = 0; i < atoms;i++) {
		
		Nucleus nucl;

		nucl.nProtons = protons;
		nucl.nFreeElec = nucl.nProtons;
		nFreeElectrons += nucl.nFreeElec; //store num electrons in field, for now summing over each nuclei.
		
		nucl.pos_lattice[0] = ((double)(i % 5))* (((double)w) / 8.0) + 0.5 * w+0.2;
		nucl.pos_lattice[1] = ((double)(i / 5)) * (((double)h) / 8.0) + 0.5 * h+0.2;
		nucl.pos_lattice[2] = ((double)l) / 2.0+0.2;

		//TODO, add angstroms when needed?
		//TODO, add pseudopotential support (charge?)
		
		Nuclei.push_back(nucl);
	}

	//init orbitals - one per electron, same size as fft_orbital
	orbitals = std::vector<std::vector<std::complex<double>>>(nFreeElectrons, std::vector<std::complex<double>>(fft_orbital.size(), 0.0));

	//Set up V_ext
	InitExternalV();
	
}

void Fields::InitExternalV() {
	//Simple sum over coulomb potentials for each nucleus. 
		//TODO : implement pseudopotentials here...
	for (int x = 0;x < w;x++) {
		for (int y = 0;y < h;y++) {
			for (int z = 0;z < l;z++) {
				for (int atom = 0; atom < Nuclei.size(); atom++) {

					double dx = cellSize* (((double)x) - Nuclei[atom].pos_lattice[0]);
					double dy = cellSize* (((double)y) - Nuclei[atom].pos_lattice[1]);
					double dz = cellSize* (((double)z) - Nuclei[atom].pos_lattice[2]);
					
					double r2 = dx * dx + dy * dy + dz * dz;
					
					//test gaussian
					//V_ext[x][y][z] = std::exp(-r2 / (2.0 * 4 * 4));

					//actual equation
					V_ext[x][y][z] += ((double)Nuclei[atom].nProtons) / std::max(std::sqrt(dx * dx + dy * dy + dz * dz),0.00001);

				}
			}
		}
	}
}