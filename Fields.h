#pragma once
#include <vector>
#include <string>
#include <complex>
#include <iostream>
#include "fftw3.h"



class Fields {
public:
	//grid dimensions
	int w=0;
	int h=0;
	int l=0;
	
	double cellSize = 0.2; // In Bohr units
	fftw_plan FFT_plan_n; //FFTW plan (for FFT(n) in V_ht solver)
	fftw_plan FFT_inv_plan_n; //Plan for inverse
	fftw_plan FFT_plan_KohnSham; //Plans for Kohn Sham equation solver ifft(orbital) --> fft(Vs*orbital)
	fftw_plan FFT_inv_plan_KohnSham;

	std::vector<std::complex<double>> FFT_out;
	 

	//electron density field -- and flattened for FFT (eventually want to flatten everything??)
	std::vector<std::vector<std::vector<double>>> n;
	std::vector<double> n_flat;

	//Electrical potential field / components
	std::vector<std::vector<std::vector<double>>> V_ext; //external (calculated from fixed point nuclei)
	std::vector<std::vector<std::vector<double>>> V_hartree; //Coulomb repulsion 
	std::vector<std::vector<std::vector<double>>> V_xc; //exchange correlation 
	std::vector<std::vector<std::vector<double>>> V_s; //total - V_ext + V_ht + V_xc

	//Orbitals - complex valued fields over space (~1 per [valence] electron)
		//A vector of complex valued grids, of length N (one per relevant electron, plus some for potential orbital jumping)
	int nFreeElectrons;
	std::vector<std::vector< std::complex<double>>> orbitals;
	std::vector<double> ifft_orbital;
	std::vector< std::complex<double>> fft_orbital;



	struct Nucleus {
		double pos_angstroms[3]; //position in angstroms (raw position)
		double pos_lattice[3]; //position in lattice units (relative to grid)
		int nProtons; //atomic number
		int nFreeElec;
		float charge;
	};
	std::vector<Nucleus> Nuclei;
	//TODO: methods, .cpp file, potentially move nuclei to MatSim.cpp, init fields from them, init nuclei in MatSim??

	//TODO:Constructors? just blank and init separate? Feed grid reference?
	Fields() {}


	void InitMaterial(std::string filename);
	void InitMaterial(int w_, int h_, int l_, double cellSize_bohr, int atoms, int protons, double charge);
	void InitExternalV();

	static void FlattenField(std::vector<std::vector<std::vector<double>>>& N, std::vector<double>& N_flat) {
		//flatten n into a row-major order 1d array
		for (int x = 0;x < N.size();x++) {
			for (int y = 0;y < N[x].size();y++) {
				for (int z = 0;z < N[x][y].size();z++) {
					N_flat[x * N[0].size() * N[0][0].size() + y * N[0][0].size() + z] = N[x][y][z];
				}
			}
		}
	}
	

};