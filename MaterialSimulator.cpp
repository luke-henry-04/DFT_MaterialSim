// MaterialSimulator.cpp : Defines the entry point for the application.
//

/*
NOTES 9/30/26

convergence is not working
seems to oscillate each iteration
looks like film grain
only shows spike in one point when nucleus is offset from the exact grid point (no dist=0 to divide by)
NEED LERP FOR N!!!!

now that its "done" lmao, i need to just go through the whole program top to bottom. look for issues

NOTES 10/7/26
It "works" if you can call it that. Results look promising, but eigenvalues are off by large factors
Check V_HT solver for missing factors, 
Check/confirm units are in Hartree/Bohr and that those are compatible??
Confirm orbital counting is correct?
Energy cutoffs are by vibe right now, will want to rework entire system into plane wave basis, (tedious but not hard)
*/

#include "MaterialSimulator.h"


using namespace std;

//set up grid? need matsim class? or do i just do the file?

int main(int argc, char* argv[])
{

	cout << "Starting.." << endl;

	Fields fields = Fields();

	//temporary debug initializations
	//eventually switch to input file parsing -> commandline filename -> GUI import file
	
	
	/*fields.w = 50;
	fields.h = 50;
	fields.l = 50;
	fields.cellSize = 0.05;
	fields.InitAtom(1, 1, -(28 * 0.5), 0, 0);
	fields.InitAtom(1, 1, (28 * 0.5), 0, 0);*/

	std::vector<std::vector<std::vector<double>>> ReferenceH2 = std::vector< std::vector< std::vector<double>>>(90, std::vector< std::vector<double>>(90,  std::vector<double>(90,0)));
	std::vector<std::vector<std::vector<double>>> ResidualH2 = ReferenceH2;
	std::ifstream file("../../../h2_lda_reference.txt");

	// Always check if the file opened successfully
	if (!file.is_open()) {
		std::cerr << "Error: Could not open the file!" << std::endl;
		return 1;
	}
	else std::cout << "Reference file opened!";

	std::string line;
	// Read line-by-line until the end of the file
	
	for (int x = 0; x < 90;x++){
		for (int y = 0; y < 90;y++) {
			for (int z = 0; z < 90;z++) {
				std::getline(file, line);
				double val = stod(line);
				ReferenceH2[x][y][z] = val;
			}
		}
	}


	if (argc == 7) {
		fields.InitMaterial(
			stoi(argv[1]),
			stoi(argv[2]),
			stoi(argv[3]),
			stof(argv[4]),
			stoi(argv[5]),
			stoi(argv[6]),
			1
		);
	}
	else fields.InitMaterial(90, 90, 90, 0.175, 2, 1, 1);
	
	DFT::NumGridPoints = fields.w * fields.h * (fields.l/2+1);
	for (int i = 0;i < 1;i++) {
		DFT::calculateV_s(fields);
		DFT::calculateOrbitals(fields, fields.FFT_plan_n, fields.FFT_inv_plan_n);
		DFT::calculateN(fields);
		std::cout << "\niteration# " << i << "\n";
	}


	for (int x = 0; x < 90;x++) {
		for (int y = 0; y < 90;y++) {
			for (int z = 0; z < 90;z++) {
				ResidualH2[x][y][z] = ReferenceH2[x][y][z] - fields.n[x][y][z];
			}
		}
	}

	//Start with initial guess for n(r) - DONE
	
	//then solve for V_s = V_ext + V_ht + V_xc - DONE
		//V_ext only solved once at start, a simple function of nuclei charges/pos ☺ 
		//V_ht is a poisson equation solution, probably want to implement a reusable conjugate gradient like solver,
			// or use a standard library??
		//V_xc is approximated as simple function of n(r) [LDA or other means]
		
	//then solve the Kohn–Sham equations for the φi orbitals - TODO
		// It is a literal matrix eigenvalue equation H*phi_i = epsilon_i * phi_i
		// H is sparce, H*phi_i(r) = (V_s + laplacian)phi_i
		// where V_s is multiplied component-wise and laplacian is a kernel convolved with phi_i
		// never actually write matrix (w*h*l x w*h*l sized), takes in whole fields as flattened vectors
		// instead use a CG or jacobi or etc solver that only needs to know certain {H*v_1, H*v_2, H*v_3, ...}
		// orthonormalization - need a gram schimdt process (??) to stop phi_i's from occupying same orbitals?
	
	//then calculate n(r) from orbitals (straightforward integral i believe)
	
	//iterate until convergence (pray for convergence)
	//dont just replace n(r) with new field, gradual lerp. (t ~ 0.1 or 0.2 ish)

	//This yields, a material - as electron orbitals/potentials --> can be saved, used in other software.
	
	//next step is to figure out inital configurations, how to set up molecules and how to group them into layers/sheets
	//make a GUI and allow editing and saving.

	//This is the DFT step.
		//this alone should be software. This alone should be polished and a project.

	//Next part will be rendering useages. 
	// Either a tweak to my material sim,
	// or a BRDF lookup table generator with a library for doing the material lookup (?)
		//Both will still need to handle photons with wavelengths and interaction. (probably best to do brdf table?)


	//DEBUGGING WINDOW -- NOT ACTUAL GUI - SWITCH TO QT EVENTUALLY
		//Inner loop shows field values with sigmoid rescaling.
	const int W = fields.w*10;
	const int H = fields.h*10;
	StartPixelWindow(W, H);
	std::vector<uint32_t> pixels(W * H);

//	int cx =128, cy = 128, cz = 128;

	/*
	// test x symmetry
	std::cout << "n x+50: " << fields.n[cx + 50][cy][cz] << " x-50: " << fields.n[cx - 50][cy][cz] << "\n";
	std::cout << "n y+50: " << fields.n[cx][cy + 50][cz] << " y-50: " << fields.n[cx][cy - 50][cz] << "\n";
	std::cout << "n z+50: " << fields.n[cx][cy][cz + 50] << " z-50: " << fields.n[cx][cy][cz - 50] << "\n";

	std::cout << "Vht x+50: " << fields.V_hartree[cx + 50][cy][cz] << " x-50: " << fields.V_hartree[cx - 50][cy][cz] << "\n";
	std::cout << "Vht y+50: " << fields.V_hartree[cx][cy + 50][cz] << " y-50: " << fields.V_hartree[cx][cy - 50][cz] << "\n";
	std::cout << "Vht z+50: " << fields.V_hartree[cx][cy][cz + 50] << " z-50: " << fields.V_hartree[cx][cy][cz - 50] << "\n";

	std::cout << fields.V_hartree[64][64][20] << " " << fields.V_hartree[64][64][236] << "\n";
	*/

	//fields.V_hartree = LinearOperator::Laplacian<double>(fields.V_hartree);
	while (true)
	{
		//z-layer loop
		for (int z = 0;z < fields.l;z++) {
			
			//Pixel loop
			for (int y = 0; y < H; y++)
			{
				for (int x = 0; x < W; x++)
				{
					//if (5000.0 * abs(fields.n[x / 4][y / 4][z]) >= 0xD3 && 5000.0 * abs(fields.n[x / 4][y / 4][z]) <= 0xD5) {
					pixels[y * W + x] =
							//(  (int) (0x0F * (abs(fields.V_hartree[y / 2][x / 2][z]/ fields.V_ext[y / 2][x / 2][z]) <15 ? abs(fields.V_hartree[y / 2][x / 2][z]) : 0)) % 0x100 );
							(((int)(000.0 * abs(ReferenceH2[x / 10][y / 10][z])) % 0xFF));

					pixels[y * W + x] +=
						//(  (int) (0x0F * (abs(fields.V_hartree[y / 2][x / 2][z]/ fields.V_ext[y / 2][x / 2][z]) <15 ? abs(fields.V_hartree[y / 2][x / 2][z]) : 0)) % 0x100 );
						(((int)(000.0 * abs(fields.n[x / 10][y / 10][z])) % 0xFF))<<8;

					pixels[y * W + x] +=
						//(  (int) (0x0F * (abs(fields.V_hartree[y / 2][x / 2][z]/ fields.V_ext[y / 2][x / 2][z]) <15 ? abs(fields.V_hartree[y / 2][x / 2][z]) : 0)) % 0x100 );
						(((int)(255.0 * abs(1.0/(1.0+exp( ResidualH2[x / 10][y / 10][z]/ ReferenceH2[x / 10][y / 10][z] )))) % 0xFF))<<16;
					//}
					//else pixels[y * W + x] = 0xFF;
				}
			}

			//draw frame and sleep 60ms
			UpdatePixelWindow(pixels.data());
			ProcessPixelWindowEvents();
			
			_sleep(20);

		}
	}


	return 0;
}
