#include <iostream>
#include <iomanip>
#include <fstream>
#include <sys/stat.h>
#include <sstream>
#include <cstdlib>
#include <algorithm>

#include "amber.h"
	
std::string const Amber::name = "Amber";

Amber::Amber() {
}

std::string const Amber::getName(){
	return "Amber";
}

void Amber::generateMdin(std::string directory){
	std::ofstream file((directory + "/mdin").c_str());
	file << " newffm Generated Geometry" << std::endl;
	file << "  &cntrl" << std::endl;
	file << "   imin=1, irest=0," << std::endl;
	file << "   ntx=1," << std::endl;
	file << "   ntb=1," << std::endl;
	file << "   cut=12.0," << std::endl;
	file << "   iamoeba=1, maxcyc=1," << std::endl;
	file << "   ntwx=1," << std::endl;
	file << "   ntwv=1" << std::endl;
	file << " /" << std::endl;
	file << "  &ewald" << std::endl;
	file << "  nfft1=80,nfft2=80,nfft3=80," << std::endl;
	file << "  skinnb=2.,nbtell=0,order=5,ew_coeff=0.45," << std::endl;
	file << " /" << std::endl;
	file << " &amoeba" << std::endl;
	file << "   do_bond=1,do_ureyb=1,do_reg_angle=1,do_trig_angle=1," << std::endl;
	file << "   do_opbend=1,do_torsion=1,do_pi_torsion=1,do_strbend=1," << std::endl;
	file << "   do_torsion_torsion=1,do_amoeba_nonbond=1," << std::endl;
	file << "   dipole_scf_tol = 0.01,dipole_scf_iter_max=100," << std::endl;
	file << "   sor_coefficient=0.6,ee_damped_cut=4.5,ee_dsum_cut=6.7," << std::endl;
	file << "   amoeba_verbose=1" << std::endl;
	file << " /" << std::endl;
	file.close();
}

void Amber::generatePrmtop(Geometry geom, std::string directory){
        int i,ii;
	mkdir(directory.c_str(),0700);
	std::ofstream outfile((directory + "/prmtop").c_str());
	outfile.setf(std::ios_base::scientific,std::ios::floatfield);

	outfile << "%VERSION VERSION_STAMP = V0001.000 DATE = 05/22/06  12:10:21" << std::endl;
	outfile << "%FLAG TITLE" << std::endl;
	outfile << "%FORMAT(20a4)" << std::endl;
	outfile << "newffm Amber Parameter File" << std::endl;

	outfile << "%FLAG POINTERS" << std::endl;
	outfile << "%FORMAT(10I8)" << std::endl;

	outfile.width(8);
	//NATOM: Total number of atoms
	outfile << geom.atoms.size();

	std::vector<std::string> types;
	bool found;
	int hydrogens = 0;
	for(std::vector<Atom>::iterator it = geom.atoms.begin(); it != geom.atoms.end(); it++){
		found = false;
		
		if(it->symbol == "H")
			hydrogens++;
		
		for(std::vector<std::string>::iterator findit = types.begin(); findit != types.end(); findit++){
			if(*findit == it->symbol){
				found = true;
				break;
			}
		}
		
		if(!found){
			types.push_back(it->symbol);
		}
	}
	
	outfile.width(8);
	//NTYPES: Total number of atom types
	outfile << types.size();

	outfile.width(8);
	//NBONH: Number of bonds containing hydrogen
	outfile << hydrogens;

	outfile.width(8);
	//MBONA: Number of bonds not containing hydrogen
	outfile << 0;

	outfile.width(8);
	//NTHETH: Number of angles containing hydrogen
	outfile << hydrogens/2;

	outfile.width(8);
	//MTHETA: Number of angles not containing hydrogen
	outfile << 0;

	outfile.width(8);
	//NPHIH: Number of dihedrals containing hydrogen
	outfile << 0;

	outfile.width(8);
	//MPHIA: Number of dihedrals not containing hydrogen
	outfile << 0;

	outfile.width(8);
	//NHPARM: Currently not used
	outfile << 0;

	outfile.width(8);
	//NPARM: Currently not used
	outfile << 0;
	outfile << std::endl;

	int nnb = 0, nres = 0;
	nnb = geom.atoms.size()-1;
	nres = geom.atoms.size()/3+1;
	outfile.width(8);
	//NNB: Number of excluded atoms
	outfile << nnb;

	outfile.width(8);
	//NRES: Number of residues;
	outfile << nres;

	outfile.width(8);
	//NBONA: MBONA + number of constraint bonds
	outfile << 0;

	outfile.width(8);
	//NTHETA: MTHETA + number of constraint angles
	outfile << 0;

	outfile.width(8);
	//NPHIA: MPHIA + number of constraint dihedrals
	outfile << 0;

	outfile.width(8);
	//NUMBND: Number of unique bond types
	outfile << 1;

	outfile.width(8);
	//NUMANG: Number of unique angle types
	outfile << 1;

	outfile.width(8);
	//NPTRA: Number of unique dihedral types
	outfile << 0;

	outfile.width(8);
	//NATYP: Number of atom types in parameter file
	outfile << types.size();

	outfile.width(8);
	//NPHB: Number of distinct 10-12 hydrogen bond pair types
	outfile << 0;
	outfile << std::endl;

	outfile.width(8);
	//IFPERT: Set to 1 if perturbation info is to be read in
	outfile << 0;

	outfile.width(8);
	//NBPER: Number of bonds to be perturbed
	outfile << 0;

	outfile.width(8);
	//NGPER: Number of angles to be perturbed
	outfile << 0;

	outfile.width(8);
	//NDPER: Number of dihedrals to be perturbed
	outfile << 0;

	outfile.width(8);
	//MBPER: Number of bonds with atoms completely in perturbed group
	outfile << 0;

	outfile.width(8);
	//MGPER: Number of angles with atoms completely in perturbed group
	outfile << 0;

	outfile.width(8);
	//MDPER: Number of dihedrals with atoms completely in perturbed groups
	outfile << 0;

	outfile.width(8);
	//IFBOX: Set to 1 if standard periodic box, 2 when truncated octahedral
	outfile << 1;

	outfile.width(8);
	//NMXRS: Number of atoms in the largest residue
	outfile << 3;

	outfile.width(8);
	//IFCAP: Set to 1 if the CAP option from edit was specified
	outfile << 0;
	outfile << std::endl;
	
	
	outfile.width(8);
	// No idea
	outfile << 0;
	outfile << std::endl;

	outfile << "%FLAG ATOM_NAME" << std::endl;
	outfile << "%FORMAT(20a4)";

	int linefactor = 0;
	for(std::vector<Atom>::iterator it = geom.atoms.begin(); it != geom.atoms.end(); it++){
		if(linefactor % 20 == 0)
			outfile << std::endl;
		outfile.width(4);
		outfile << it->symbol;
		linefactor++;
	}
	outfile << std::endl;

	outfile << "%FLAG CHARGE" << std::endl;
	outfile << "%FORMAT(5E16.8)";

	linefactor = 0;
	for(std::vector<Atom>::iterator it = geom.atoms.begin(); it != geom.atoms.end(); it++){
		if(linefactor % 5 == 0) 
			outfile << std::endl;
		outfile.width(16);
		outfile.precision(8);
		outfile << it->charge;
		linefactor++;
	}
	outfile << std::endl;

	outfile << "%FLAG MASS" << std::endl;
	outfile << "%FORMAT(5E16.8)";

	linefactor = 0;
	for(std::vector<Atom>::iterator it = geom.atoms.begin(); it != geom.atoms.end(); it++){
		if(linefactor % 5 == 0) 
			outfile << std::endl;
		outfile.width(16);
		outfile.precision(8);
		outfile << it->mass;
		linefactor++;
	}
	outfile << std::endl;

	outfile << "%FLAG RADII" << std::endl;
	outfile << "%FORMAT(5E16.8)";

	linefactor = 0;
	for(std::vector<Atom>::iterator it = geom.atoms.begin(); it != geom.atoms.end(); it++){
		if(linefactor % 5 == 0) 
			outfile << std::endl;
		outfile.width(16);
		outfile.precision(8);
		outfile << it->mass/2.0;
		linefactor++;
	}
	outfile << std::endl;

	outfile << "%FLAG ATOM_TYPE_INDEX" << std::endl;
	outfile << "%FORMAT(10I8)";

	linefactor = 0;
	for(std::vector<Atom>::iterator atomit = geom.atoms.begin(); atomit != geom.atoms.end(); atomit++){
		int typeNum = 1;
		for(std::vector<std::string>::iterator it = types.begin(); it != types.end(); it++){
			if(atomit->symbol == *it){
				if(linefactor % 10 == 0)
					outfile << std::endl;
				outfile.width(8);
				outfile << typeNum;
				linefactor++;
			}
			typeNum++;
		}
	}
	outfile << std::endl;

	outfile << "%FLAG NUMBER_EXCLUDED_ATOMS" << std::endl;
	outfile << "%FORMAT(10I8)";

	linefactor = 0;
	bool nextHyd = true;
	for(std::vector<Atom>::iterator it = geom.atoms.begin(); it != geom.atoms.end(); it++){
		if(linefactor % 10 == 0) 
			outfile << std::endl;
		outfile.width(8);
		if(it->symbol == "O")
			outfile << 2;
		else if(it->symbol == "H" && nextHyd)
			outfile << 1;
		else
			outfile << 0;

		if(it->symbol == "H")
			nextHyd = !nextHyd;

		linefactor++;
	}
	outfile << std::endl;

	outfile << "%FLAG NONBONDED_PARM_INDEX" << std::endl;
	outfile << "%FORMAT(10I8)" << std::endl;
	outfile << std::setw(8) << 1 << std::setw(8) << 2 << std::setw(8) << 4 << std::setw(8) << 2 << std::setw(8) << 3;
	outfile << std::setw(8) << 5 << std::setw(8) << 4 << std::setw(8) << 5 << std::setw(8) << 6 << std::endl;

	outfile << "%FLAG RESIDUE_LABEL" << std::endl;
	outfile << "%FORMAT(20a4)";
	std::vector<int> residuePointers;
	int waterCnt = 0, atomNum = 1;
	linefactor = 0;
	for(std::vector<Atom>::iterator it = geom.atoms.begin(); it != geom.atoms.end(); it++){
		outfile.width(4);

		if(it->symbol != "O" && it->symbol != "H"){
			if(linefactor % 20 == 0) 
				outfile << std::endl;

			outfile << it->symbol;
			linefactor++;

			residuePointers.push_back(atomNum);

			if(waterCnt != 0)
				std::cout << "Warning: File is not in amoeba format. 1" << std::endl;
		} else {
			if(waterCnt == 0)
				residuePointers.push_back(atomNum);
			waterCnt++;
		}

		if(waterCnt == 3){
			if(linefactor % 20 == 0) 
				outfile << std::endl;
			outfile << "HOH";
			linefactor++;
			waterCnt = 0;
		}

		atomNum++;
	}
	outfile << std::endl;

	outfile << "%FLAG RESIDUE_POINTER" << std::endl;
	outfile << "%FORMAT(10I8)";

	linefactor = 0;
	for(std::vector<int>::iterator it = residuePointers.begin(); it != residuePointers.end(); it++){
		if(linefactor % 10 == 0) 
			outfile << std::endl;
		outfile.width(8);
		outfile << *it;
		linefactor++;
	}
	outfile << std::endl;

	/* TODO
	outfile << "%FLAG BOND_FORCE_CONSTANT" << std::endl;
	outfile << "%FORMAT(5E16.8)";
	linefactor = 0;
	//Iterate through all Coeff 1's for each atom-atom potential pair
	for(i = 0; i < potentials.bondCoeffs.size(); i++){
		if(linefactor % 5 == 0)
			outfile << std::endl;

		outfile.width(16);
		outfile.precision(8);
		float value;
		std::istringstream iss(potentials.bondCoeffs[i][0]);
		if(iss >> value)
			outfile << value;
		else
			outfile << variables[potentials.bondCoeffs[i][0][0] - 65];

		linefactor++;
	}
	outfile << std::endl;

	outfile << "%FLAG BOND_EQUIL_VALUE" << std::endl;
	outfile << "%FORMAT(5E16.8)";
	linefactor = 0;
	//Iterate through all Coeff 2's for each atom-atom potential pair
	for(i = 0; i < potentials.bondCoeffs.size(); i++){
		if(linefactor % 5 == 0)
			outfile << std::endl;

		outfile.width(16);
		outfile.precision(8);
		float value;
		std::istringstream iss(potentials.bondCoeffs[i][1]);
		if(iss >> value)
			outfile << value;
		else
			outfile << variables[potentials.bondCoeffs[i][1][0] - 65];

		linefactor++;
	}
	outfile << std::endl;

	outfile << "%FLAG ANGLE_FORCE_CONSTANT" << std::endl;
	outfile << "%FORMAT(5E16.8)";
	linefactor = 0;
	//Iterate through all Coeff 1's for each atom-atom-atom potential set
	for(i = 0; i < potentials.angleCoeffs.size(); i++){
		if(linefactor % 5 == 0)
			outfile << std::endl;

		outfile.width(16);
		outfile.precision(8);
		float value;
		std::istringstream iss(potentials.angleCoeffs[i][0]);
		if(iss >> value)
			outfile << value;
		else
			outfile << variables[potentials.angleCoeffs[i][0][0] - 65];

		linefactor++;
	}
	outfile << std::endl;

	outfile << "%FLAG ANGLE_EQUIL_VALUE" << std::endl;
	outfile << "%FORMAT(5E16.8)";
	linefactor = 0;
	//Iterate through all Coeff 2's for each atom-atom-atom potential set
	for(i = 0; i < potentials.angleCoeffs.size(); i++){
		if(linefactor % 5 == 0)
			outfile << std::endl;

		outfile.width(16);
		outfile.precision(8);
		float value;
		std::istringstream iss(potentials.angleCoeffs[i][1]);
		if(iss >> value)
			outfile << value;
		else
			outfile << variables[potentials.angleCoeffs[i][1][0] - 65];

		linefactor++;
	}
	outfile << std::endl;
	*/

	outfile << "%FLAG DIHEDRAL_FORCE_CONSTANT" << std::endl;
	outfile << "%FORMAT(5E16.8)" << std::endl;
	outfile << std::endl;

	outfile << "%FLAG DIHEDRAL_PERIODICITY" << std::endl;
	outfile << "%FORMAT(5E16.8)" << std::endl;
	outfile << std::endl;

	outfile << "%FLAG DIHEDRAL_PHASE" << std::endl;
	outfile << "%FORMAT(5E16.8)" << std::endl;
	outfile << std::endl;

	outfile << "%FLAG SOLTY" << std::endl;
	outfile << "%FORMAT(5E16.8)" << std::endl;
	outfile << "  0.00000000e+00  0.00000000e+00  0.00000000e+00" << std::endl;

	outfile << "%FLAG LENNARD_JONES_ACOEF" << std::endl;
	outfile << "%FORMAT(5E16.8)" << std::endl;
	outfile << "  7.18764375e+11  2.92107835e+10  1.09435450e+09  2.84889794e+09  9.45223949e+07" << std::endl;
	outfile << "  6.78372168e+06" << std::endl;

	outfile << "%FLAG LENNARD_JONES_BCOEF" << std::endl;
	outfile << "%FORMAT(5E16.8)" << std::endl;
	outfile << "  5.36195627e+06  3.50065869e+05  2.19434724e+04  6.47070921e+04  3.81705983e+03" << std::endl;
	outfile << "  6.05244554e+02" << std::endl;

	outfile << "%FLAG BONDS_INC_HYDROGEN" << std::endl;
	outfile << "%FORMAT(10I8)";

	linefactor = 0;
	std::vector<int> oxyStack, hydStack;
	atomNum = 1;
	for(std::vector<Atom>::iterator it = geom.atoms.begin(); it != geom.atoms.end(); it++){
		if(it->symbol == "H"){
			hydStack.push_back(atomNum);
		} else if(it->symbol == "O"){
			oxyStack.push_back(atomNum);
			oxyStack.push_back(atomNum);
		}

		atomNum++;
	}

	while(!hydStack.empty() && !oxyStack.empty()){
		if(linefactor % 10 == 0) 
			outfile << std::endl;
		outfile.width(8);
		outfile << (oxyStack.front()-1)*3;
		linefactor++;

		if(linefactor % 10 == 0) 
			outfile << std::endl;
		outfile.width(8);
		outfile << (hydStack.front()-1)*3;
		linefactor++;

		if(linefactor % 10 == 0) 
			outfile << std::endl;
		outfile.width(8);
		outfile << 1;

		linefactor++;
		oxyStack.erase(oxyStack.begin());
		hydStack.erase(hydStack.begin());
	}
	outfile << std::endl;

	if(!hydStack.empty() || !oxyStack.empty())
		std::cout << "Warning: File is not in amoeba format. 2" << std::endl;
	
	outfile << "%FLAG BONDS_WITHOUT_HYDROGEN" << std::endl;
	outfile << "%FORMAT(10I8)" << std::endl;
	outfile << std::endl;

	outfile << "%FLAG ANGLES_INC_HYDROGEN" << std::endl;
	outfile << "%FORMAT(10I8)";

	//	linefactor = 0;
	//i = 1;
	//for(std::vector<Atom>::iterator it = geom.atoms.begin(); it != geom.atoms.end(); it++){
	//	if(it->symbol == "H" || it->symbol == "O"){
	//		if(linefactor % 10 == 0) 
	//			outfile << std::endl;
	//		outfile.width(8);
	//		outfile << (i-1)*3;
	//		linefactor++;

	//		if(linefactor % 3 == 0){
	//			if(linefactor % 10 == 0) 
	//				outfile << std::endl;
	//			outfile.width(8);
	//			outfile << 1;
	//			linefactor++;
	//		}
	//	}
	//	i++;
	//}
	//outfile << std::endl; 

	linefactor = 0;
	/* TODO
	for(i = 0; i < potentials.angleAtoms.size(); i++){
		if(potentials.angleCoeffType[potentials.angleCoeffNums[i]] != "REGULAR")
			continue;

		if(linefactor % 10 == 0)
			outfile << std::endl;

		outfile.width(8);
		outfile << (atoi(potentials.angleAtoms[i][0].c_str())-1)*3;
		linefactor++;
		
		if(linefactor % 10 == 0)
			outfile << std::endl;

		outfile.width(8);
		outfile << (atoi(potentials.angleAtoms[i][1].c_str())-1)*3;
		linefactor++;
		
		if(linefactor % 10 == 0)
			outfile << std::endl;

		outfile.width(8);
		outfile << (atoi(potentials.angleAtoms[i][2].c_str())-1)*3;
		linefactor++;

		if(linefactor % 10 == 0)
			outfile << std::endl;

		outfile.width(8);
		outfile << 1;
		linefactor++;
	}
	outfile << std::endl;
	*/


	outfile << "%FLAG ANGLES_WITHOUT_HYDROGEN" << std::endl;
	outfile << "%FORMAT(10I8)" << std::endl;
	outfile << std::endl;

	outfile << "%FLAG DIHEDRALS_INC_HYDROGEN" << std::endl;
	outfile << "%FORMAT(10I8)" << std::endl;
	outfile << std::endl;
	
	outfile << "%FLAG DIHEDRALS_WITHOUT_HYDROGEN" << std::endl;
	outfile << "%FORMAT(10I8)" << std::endl;
	outfile << std::endl;

	outfile << "%FLAG EXCLUDED_ATOMS_LIST" << std::endl;
	outfile << "%FORMAT(10I8)";
	linefactor = 0;
	i = 1;
	int writeTwo = true;
	for(std::vector<Atom>::iterator it = geom.atoms.begin(); it != geom.atoms.end(); it++){
		if(it->symbol == "H"){
			if(writeTwo){
				if(linefactor % 10 == 0) 
					outfile << std::endl;
				outfile.width(8);
				outfile << i;
				linefactor++;
			}	
			
			if(linefactor % 10 == 0) 
				outfile << std::endl;
			outfile.width(8);
			outfile << i;
			linefactor++;
			
			writeTwo = !writeTwo;
		} else if(it->symbol == "O"){
			writeTwo = false;
		}
		i++;
	}
	outfile << std::endl;
	
	outfile << "%FLAG HBOND_ACOEF" << std::endl;
	outfile << "%FORMAT(5E16.8)" << std::endl;
	outfile << std::endl;

	outfile << "%FLAG HBOND_BCOEF" << std::endl;
	outfile << "%FORMAT(5E16.8)" << std::endl;
	outfile << std::endl;

	outfile << "%FLAG HBCUT" << std::endl;
	outfile << "%FORMAT(5E16.8)" << std::endl;
	outfile << std::endl;

	outfile << "%FLAG AMBER_ATOM_TYPE" << std::endl;
	outfile << "%FORMAT(20a4)";

	linefactor = 0;
	for(std::vector<Atom>::iterator it = geom.atoms.begin(); it != geom.atoms.end(); it++){
		if(linefactor % 20 == 0)
			outfile << std::endl;
		outfile.width(4);
		if(it->symbol == "H"){
			outfile << 67;
		} else if(it->symbol == "O"){
			outfile << 66;
		} else {
			outfile << 10;
		}
		linefactor++;
	}
	outfile << std::endl;

	outfile << "%FLAG TREE_CHAIN_CLASSIFICATION" << std::endl;
	outfile << "%FORMAT(20a4)";

	linefactor = 0;
	for(std::vector<Atom>::iterator it = geom.atoms.begin(); it != geom.atoms.end(); it++){
		if(linefactor % 20 == 0)
			outfile << std::endl;
		outfile.width(4);
		outfile << "BLA";
		linefactor++;
	}
	outfile << std::endl;

	outfile << "%FLAG JOIN_ARRAY" << std::endl;
	outfile << "%FORMAT(10I8)";

	linefactor = 0;
	for(std::vector<Atom>::iterator it = geom.atoms.begin(); it != geom.atoms.end(); it++){
		if(linefactor % 10 == 0) 
			outfile << std::endl;
		outfile.width(8);
		outfile << 0;

		linefactor++;
	}
	outfile << std::endl;

	outfile << "%FLAG IROTAT" << std::endl;
	outfile << "%FORMAT(10I8)";

	linefactor = 0;
	for(std::vector<Atom>::iterator it = geom.atoms.begin(); it != geom.atoms.end(); it++){
		if(linefactor % 10 == 0) 
			outfile << std::endl;
		outfile.width(8);
		outfile << 0;

		linefactor++;
	}
	outfile << std::endl;

	outfile << "%FLAG SOLVENT_POINTERS" << std::endl;
	outfile << "%FORMAT(10I8)" << std::endl;
	outfile.width(8);
	outfile << geom.atoms.size()/3+1;
	outfile.width(8);
	outfile << geom.atoms.size()/3;
	outfile << "       2       1       2" << std::endl;
	outfile << "%FLAG ATOMS_PER_MOLECULE" << std::endl;
	outfile << "%FORMAT(10I8)";
	linefactor = 0;
	for(i = 0; i < geom.atoms.size()/3; i++){
		if(linefactor % 10 == 0) 
			outfile << std::endl;
		outfile.width(8);
		outfile << 3;

		linefactor++;
	}
	if(linefactor % 10 != 0)
		outfile << std::endl;
	outfile << "%FLAG BOX_DIMENSIONS" << std::endl;
	outfile << "%FORMAT(5E16.8)" << std::endl;
	outfile << "  9.00000000e+01  4.00000000e+01  4.00000000e+01  4.00000000e+01" << std::endl;

	outfile << "%FLAG AMOEBA_FORCEFIELD" << std::endl;
	outfile << "%COMMENT This indicates that this parm file is specific to amoeba" << std::endl;
	outfile << "%COMMENT This must be present if do_amoeba(in mdin) is 1" << std::endl;
	outfile << "%COMMENT This must NOT be present if do_amoeba is 0" << std::endl;
	outfile << "%FORMAT(i5)" << std::endl;
	outfile << "    1" << std::endl;
	
	outfile << "%FLAG AMOEBA_ATOM_TYPE_INDEX" << std::endl;
	outfile << "%COMMENT   dimension = (4)" << std::endl;
	outfile << "%FORMAT(10I8)";
	linefactor = 0;
	for(std::vector<Atom>::iterator it = geom.atoms.begin(); it != geom.atoms.end(); it++){
		if(linefactor % 10 == 0)
			outfile << std::endl;
		outfile.width(8);
		if(it->symbol == "H"){
			outfile << 285;
		} else if(it->symbol == "O"){
			outfile << 284;
		} else {
			outfile << 10;
		}
		linefactor++;
	}
	outfile << std::endl;
	
	outfile << "%FLAG AMOEBA_ATOMIC_NUMBER" << std::endl;
	outfile << "%COMMENT   dimension = (4)" << std::endl;
	outfile << "%FORMAT(10I8)";
	linefactor = 0;
	for(std::vector<Atom>::iterator it = geom.atoms.begin(); it != geom.atoms.end(); it++){
		if(linefactor % 10 == 0)
			outfile << std::endl;
		outfile.width(8);
		outfile << it->number;
		linefactor++;
	}
	outfile << std::endl;
	
	outfile << "%FLAG AMOEBA_ATOM_CLASS_INDEX" << std::endl;
	outfile << "%COMMENT   dimension = (4)" << std::endl;
	outfile << "%FORMAT(10I8)";
	linefactor = 0;
	for(std::vector<Atom>::iterator it = geom.atoms.begin(); it != geom.atoms.end(); it++){
		if(linefactor % 10 == 0)
			outfile << std::endl;
		outfile.width(8);
		if(it->symbol == "H"){
			outfile << 67;
		} else if(it->symbol == "O"){
			outfile << 66;
		} else {
			outfile << 10;
		}
		linefactor++;
	}
	outfile << std::endl;
	
	outfile << "%FLAG AMOEBA_REGULAR_BOND_NUM_LIST" << std::endl;
	outfile << "%FORMAT(I8)" << std::endl;
	outfile.width(8);
	outfile << (geom.atoms.size()-1)/3*2 << std::endl;
	
	/* TODO
	outfile << "%FLAG AMOEBA_REGULAR_BOND_LIST" << std::endl;
	outfile << "%COMMENT dimension = (3,2)" << std::endl;
	outfile << "%FORMAT(10I8)";
	linefactor = 0;
	for(i = 0; i < potentials.bondAtoms.size(); i++){
		if(potentials.bondCoeffType[potentials.bondCoeffNums[i]] != "REGULAR")
			continue;

		if(linefactor % 10 == 0)
			outfile << std::endl;

		outfile.width(8);
		outfile << potentials.bondAtoms[i][0];
		linefactor++;
		
		if(linefactor % 10 == 0)
			outfile << std::endl;

		outfile.width(8);
		outfile << potentials.bondAtoms[i][1];
		linefactor++;

		if(linefactor % 10 == 0)
			outfile << std::endl;

		outfile.width(8);
		outfile << potentials.bondCoeffNums[i]+1;
		linefactor++;
	}
	outfile << std::endl;
	*/
	/*
	linefactor = 0;
	std::vector<int> oxyStack, hydStack;
	int atomNum = 1;
	for(std::vector<Atom>::iterator it = geom.atoms.begin(); it != geom.atoms.end(); it++){
		if(it->symbol == "H"){
			hydStack.push_back(atomNum);
		} else if(it->symbol == "O"){
			oxyStack.push_back(atomNum);
			oxyStack.push_back(atomNum);
		}

		atomNum++;
	}

	while(!hydStack.empty() && !oxyStack.empty()){
		if(linefactor % 10 == 0) 
			outfile << std::endl;
		outfile.width(8);
		outfile << oxyStack.front();
		linefactor++;

		if(linefactor % 10 == 0) 
			outfile << std::endl;
		outfile.width(8);
		outfile << hydStack.front();
		linefactor++;

		if(linefactor % 10 == 0) 
			outfile << std::endl;
		outfile.width(8);
		outfile << 1;

		linefactor++;
		oxyStack.erase(oxyStack.begin());
		hydStack.erase(hydStack.begin());
	}
	if(linefactor % 10 != 1)
		outfile << std::endl;

	if(!hydStack.empty() || !oxyStack.empty())
		std::cout << "Warning: File is not in amoeba format." << std::endl;
	*/

	/* TODO	
	outfile << "%FLAG AMOEBA_REGULAR_BOND_NUM_PARAMS" << std::endl;
	outfile << "%FORMAT(I8)" << std::endl;
	int regular_bond_count = 0;
	for(i = 0; i < potentials.bondCoeffType.size(); i++){
		if(potentials.bondCoeffType[i] == "REGULAR")
			regular_bond_count++;
	}
	outfile << "     " << regular_bond_count << std::endl;
	
	outfile << "%FLAG AMOEBA_REGULAR_BOND_FORCE_CONSTANT" << std::endl;
	outfile << "%FORMAT(5E16.8)";
	linefactor = 0;
	//Iterate through all Coeff 1's for each atom-atom potential pair
	for(i = 0; i < potentials.bondCoeffs.size(); i++){
		if(potentials.bondCoeffType[i] != "REGULAR")
			continue;

		if(linefactor % 5 == 0)
			outfile << std::endl;

		outfile.width(16);
		outfile.precision(8);
		float value;
		std::istringstream iss(potentials.bondCoeffs[i][0]);
		if(iss >> value)
			outfile << value;
		else
			outfile << variables[potentials.bondCoeffs[i][0][0] - 65];

		linefactor++;
	}
	outfile << std::endl;
	
	outfile << "%FLAG AMOEBA_REGULAR_BOND_EQUIL_VALUE" << std::endl;
	outfile << "%FORMAT(5E16.8)";
	linefactor = 0;
	//Iterate through all Coeff 2's for each atom-atom potential pair
	for(i = 0; i < potentials.bondCoeffs.size(); i++){
		if(potentials.bondCoeffType[i] != "REGULAR")
			continue;

		if(linefactor % 5 == 0)
			outfile << std::endl;

		outfile.width(16);
		outfile.precision(8);
		float value;
		std::istringstream iss(potentials.bondCoeffs[i][1]);
		if(iss >> value)
			outfile << value;
		else
			outfile << variables[potentials.bondCoeffs[i][1][0] - 65];

		linefactor++;
	}
	outfile << std::endl;
	*/
	
	outfile << "%FLAG AMOEBA_REGULAR_BOND_FTAB_DEGREE" << std::endl;
	outfile << "%FORMAT(I8)" << std::endl;
	outfile << "       4" << std::endl;
	
	outfile << "%FLAG AMOEBA_REGULAR_BOND_FTAB_COEFFS" << std::endl;
	outfile << "%FORMAT(5E16.8)" << std::endl;
	outfile << "  0.00000000e+00  0.00000000e+00  1.00000000e+00 -2.55000000e+00  3.79312500e+00" << std::endl;
	
	/* TODO
	outfile << "%FLAG AMOEBA_UREY_BRADLEY_BOND_NUM_LIST" << std::endl;
	outfile << "%FORMAT(I8)" << std::endl;
	int urey_bradley_num = 0;
	for(i = 0; i < potentials.bondAtoms.size(); i++){
		if(potentials.bondCoeffType[potentials.bondCoeffNums[i]] == "UREY_BRADLEY")
			urey_bradley_num++;
	}
	outfile << "     " << urey_bradley_num << std::endl;

	outfile << "%FLAG AMOEBA_UREY_BRADLEY_BOND_LIST" << std::endl;
	outfile << "%COMMENT   dimension = (3,1)" << std::endl;
	outfile << "%FORMAT(10I8)";
	linefactor = 0;
	for(i = 0; i < potentials.bondAtoms.size(); i++){
		if(potentials.bondCoeffType[potentials.bondCoeffNums[i]] != "UREY_BRADLEY")
			continue;

		if(linefactor % 10 == 0)
			outfile << std::endl;

		outfile.width(8);
		outfile << potentials.bondAtoms[i][0];
		linefactor++;
		
		if(linefactor % 10 == 0)
			outfile << std::endl;

		outfile.width(8);
		outfile << potentials.bondAtoms[i][1];
		linefactor++;

		if(linefactor % 10 == 0)
			outfile << std::endl;

		outfile.width(8);
		outfile << potentials.bondCoeffNums[i];
		linefactor++;
	}
	outfile << std::endl;
	*/
	/*
	linefactor = 0;
	atomNum = 1;
	int hydNum = 0;
	for(std::vector<Atom>::iterator it = geom.atoms.begin(); it != geom.atoms.end(); it++){
		if(it->symbol == "H"){
			if(hydNum == 0){
				hydNum = atomNum;
			} else {
				if(linefactor % 10 == 0) 
					outfile << std::endl;
				outfile.width(8);
				outfile << hydNum;
				linefactor++;

				if(linefactor % 10 == 0) 
					outfile << std::endl;
				outfile.width(8);
				outfile << atomNum;
				linefactor++;

				if(linefactor % 10 == 0) 
					outfile << std::endl;
				outfile.width(8);
				outfile << 1;

				linefactor++;
				hydNum = 0;
			}
		}
		atomNum++;
	}
	if(linefactor % 10 != 1)
		outfile << std::endl;
	*/
	
	/* TODO
	outfile << "%FLAG AMOEBA_UREY_BRADLEY_BOND_NUM_PARAMS" << std::endl;
	outfile << "%FORMAT(I8)" << std::endl;
	int urey_bond_count = 0;
	for(i = 0; i < potentials.bondCoeffType.size(); i++){
		if(potentials.bondCoeffType[i] == "UREY_BRADLEY")
			urey_bond_count++;
	}
	outfile << "    " << urey_bond_count << std::endl;
	
	outfile << "%FLAG AMOEBA_UREY_BRADLEY_BOND_FORCE_CONSTANT" << std::endl;
	outfile << "%FORMAT(5E16.8)";
	linefactor = 0;
	//Iterate through all Coeff 1's for each atom-atom potential pair
	for(i = 0; i < potentials.bondCoeffs.size(); i++){
		if(potentials.bondCoeffType[i] != "UREY_BRADLEY")
			continue;

		if(linefactor % 5 == 0)
			outfile << std::endl;

		outfile.width(16);
		outfile.precision(8);
		float value;
		std::istringstream iss(potentials.bondCoeffs[i][0]);
		if(iss >> value)
			outfile << value;
		else
			outfile << variables[potentials.bondCoeffs[i][0][0] - 65];

		linefactor++;
	}
	outfile << std::endl;
	
	outfile << "%FLAG AMOEBA_UREY_BRADLEY_BOND_EQUIL_VALUE" << std::endl;
	outfile << "%FORMAT(5E16.8)";
	linefactor = 0;
	//Iterate through all Coeff 1's for each atom-atom potential pair
	for(i = 0; i < potentials.bondCoeffs.size(); i++){
		if(potentials.bondCoeffType[i] != "UREY_BRADLEY")
			continue;

		if(linefactor % 5 == 0)
			outfile << std::endl;

		outfile.width(16);
		outfile.precision(8);
		float value;
		std::istringstream iss(potentials.bondCoeffs[i][1]);
		if(iss >> value)
			outfile << value;
		else
			outfile << variables[potentials.bondCoeffs[i][1][0] - 65];

		linefactor++;
	}
	outfile << std::endl;
	*/
	
	outfile << "%FLAG AMOEBA_UREY_BRADLEY_BOND_FTAB_DEGREE" << std::endl;
	outfile << "%FORMAT(I8)" << std::endl;
	outfile << "       2" << std::endl;
	
	outfile << "%FLAG AMOEBA_UREY_BRADLEY_BOND_FTAB_COEFFS" << std::endl;
	outfile << "%FORMAT(5E16.8)" << std::endl;
	outfile << "  0.00000000e+00  0.00000000e+00  1.00000000e+00" << std::endl;

	/*	
	outfile << "%FLAG AMOEBA_REGULAR_ANGLE_NUM_LIST" << std::endl;
	outfile << "%FORMAT(I8)" << std::endl;
	int regular_angle_num = 0;
	for(i = 0; i < potentials.angleAtoms.size(); i++){
		if(potentials.angleCoeffType[potentials.angleCoeffNums[i]] == "REGULAR")
			regular_angle_num++;
	}
	outfile << "     " << regular_angle_num << std::endl;
	
	outfile << "%FLAG AMOEBA_REGULAR_ANGLE_LIST" << std::endl;
	outfile << "%COMMENT dimension = (4,1)" << std::endl;
	outfile << "%FORMAT(10I8)";

	linefactor = 0;
	for(i = 0; i < potentials.angleAtoms.size(); i++){
		if(potentials.angleCoeffType[potentials.angleCoeffNums[i]] != "REGULAR")
			continue;

		if(linefactor % 10 == 0)
			outfile << std::endl;

		outfile.width(8);
		outfile << potentials.angleAtoms[i][0];
		linefactor++;
		
		if(linefactor % 10 == 0)
			outfile << std::endl;

		outfile.width(8);
		outfile << potentials.angleAtoms[i][1];
		linefactor++;
		
		if(linefactor % 10 == 0)
			outfile << std::endl;

		outfile.width(8);
		outfile << potentials.angleAtoms[i][2];
		linefactor++;

		if(linefactor % 10 == 0)
			outfile << std::endl;

		outfile.width(8);
		outfile << potentials.angleCoeffNums[i]+1;
		linefactor++;
	}
	outfile << std::endl;
	*/

	/*
	linefactor = 0;
	i = 1;
	for(std::vector<Atom>::iterator it = geom.atoms.begin(); it != geom.atoms.end(); it++){
		if(it->symbol == "H" || it->symbol == "O"){
			if(linefactor % 10 == 0) 
				outfile << std::endl;
			outfile.width(8);
			outfile << i;
			linefactor++;

			if(linefactor % 3 == 0){
				if(linefactor % 10 == 0) 
					outfile << std::endl;
				outfile.width(8);
				outfile << 1;
				linefactor++;
			}
		}
		i++;
	}
	if(linefactor % 10 != 1)
		outfile << std::endl;
	*/
	
	/* TODO
	outfile << "%FLAG AMOEBA_REGULAR_ANGLE_NUM_PARAMS" << std::endl;
	outfile << "%FORMAT(I8)" << std::endl;
	int regular_angle_count = 0;
	for(i = 0; i < potentials.angleCoeffType.size(); i++){
		if(potentials.angleCoeffType[i] == "REGULAR")
			regular_angle_count++;
	}
	outfile << "     " << regular_angle_count << std::endl;
	
	outfile << "%FLAG AMOEBA_REGULAR_ANGLE_FORCE_CONSTANT" << std::endl;
	outfile << "%FORMAT(5E16.8)";
	linefactor = 0;
	//Iterate through all Coeff 1's for each atom-atom-atom potential set
	for(i = 0; i < potentials.angleCoeffs.size(); i++){
		if(potentials.angleCoeffType[i] != "REGULAR")
			continue;

		if(linefactor % 5 == 0)
			outfile << std::endl;

		outfile.width(16);
		outfile.precision(8);
		float value;
		std::istringstream iss(potentials.angleCoeffs[i][0]);
		if(iss >> value)
			outfile << value;
		else
			outfile << variables[potentials.angleCoeffs[i][0][0] - 65];

		linefactor++;
	}
	outfile << std::endl;
	
	outfile << "%FLAG AMOEBA_REGULAR_ANGLE_EQUIL_VALUE" << std::endl;
	outfile << "%FORMAT(5E16.8)";
	linefactor = 0;
	//Iterate through all Coeff 2's for each atom-atom-atom potential set
	for(i = 0; i < potentials.angleCoeffs.size(); i++){
		if(potentials.angleCoeffType[i] != "REGULAR")
			continue;

		if(linefactor % 5 == 0)
			outfile << std::endl;

		outfile.width(16);
		outfile.precision(8);
		float value;
		std::istringstream iss(potentials.angleCoeffs[i][1]);
		if(iss >> value)
			outfile << value;
		else
			outfile << variables[potentials.angleCoeffs[i][1][0] - 65];

		linefactor++;
	}
	outfile << std::endl;

	outfile << "%FLAG AMOEBA_REGULAR_ANGLE_FTAB_DEGREE" << std::endl;
	outfile << "%FORMAT(I8)" << std::endl;
	outfile << "       6" << std::endl;
	outfile << "%FLAG AMOEBA_REGULAR_ANGLE_FTAB_COEFFS" << std::endl;
	outfile << "%FORMAT(5E16.8)" << std::endl;
	outfile << "  0.00000000e+00  0.00000000e+00  1.00000000e+00 -1.40000000e-02  5.60000000e-05" << std::endl;
	outfile << " -7.00000000e-07  2.20000000e-08" << std::endl;
	outfile << "%FLAG AMOEBA_TRIGONAL_ANGLE_NUM_LIST" << std::endl;
	outfile << "%FORMAT(I8)" << std::endl;
	int trigonal_angle_num = 0;
	for(i = 0; i < potentials.angleAtoms.size(); i++){
		if(potentials.angleCoeffType[potentials.angleCoeffNums[i]] == "TRIGONAL")
			trigonal_angle_num++;
	}
	outfile << "     " << trigonal_angle_num << std::endl;

	outfile << "%FLAG AMOEBA_TRIGONAL_ANGLE_LIST" << std::endl;
	outfile << "%COMMENT dimension = (5,0)" << std::endl;
	outfile << "%FORMAT(10I8)";
	linefactor = 0;
	for(i = 0; i < potentials.angleAtoms.size(); i++){
		if(potentials.angleCoeffType[potentials.angleCoeffNums[i]] != "TRIGONAL")
			continue;

		if(linefactor % 10 == 0)
			outfile << std::endl;

		outfile.width(8);
		outfile << potentials.angleAtoms[i][0];
		linefactor++;
		
		if(linefactor % 10 == 0)
			outfile << std::endl;

		outfile.width(8);
		outfile << potentials.angleAtoms[i][1];
		linefactor++;
		
		if(linefactor % 10 == 0)
			outfile << std::endl;

		outfile.width(8);
		outfile << potentials.angleAtoms[i][2];
		linefactor++;

		if(linefactor % 10 == 0)
			outfile << std::endl;

		outfile.width(8);
		outfile << potentials.angleCoeffNums[i];
		linefactor++;
	}
	outfile << std::endl;

	outfile << "%FLAG AMOEBA_TRIGONAL_ANGLE_NUM_PARAMS" << std::endl;
	outfile << "%FORMAT(I8)" << std::endl;
	int trigonal_angle_count = 0;
	for(i = 0; i < potentials.angleCoeffType.size(); i++){
		if(potentials.angleCoeffType[i] == "TRIGONAL")
			trigonal_angle_count++;
	}
	outfile << "    " << trigonal_angle_count << std::endl;

	outfile << "%FLAG AMOEBA_TRIGONAL_ANGLE_FORCE_CONSTANT" << std::endl;
	outfile << "%FORMAT(5E16.8)";
	linefactor = 0;
	//Iterate through all Coeff 1's for each atom-atom-atom potential set
	for(i = 0; i < potentials.angleCoeffs.size(); i++){
		if(potentials.angleCoeffType[i] != "TRIGONAL")
			continue;

		if(linefactor % 5 == 0)
			outfile << std::endl;

		outfile.width(16);
		outfile.precision(8);
		float value;
		std::istringstream iss(potentials.angleCoeffs[i][0]);
		if(iss >> value)
			outfile << value;
		else
			outfile << variables[potentials.angleCoeffs[i][0][0] - 65];

		linefactor++;
	}
	outfile << std::endl;

	outfile << "%FLAG AMOEBA_TRIGONAL_ANGLE_EQUIL_VALUE" << std::endl;
	outfile << "%FORMAT(5E16.8)";
	linefactor = 0;
	//Iterate through all Coeff 2's for each atom-atom-atom potential set
	for(i = 0; i < potentials.angleCoeffs.size(); i++){
		if(potentials.angleCoeffType[i] != "TRIGONAL")
			continue;

		if(linefactor % 5 == 0)
			outfile << std::endl;

		outfile.width(16);
		outfile.precision(8);
		float value;
		std::istringstream iss(potentials.angleCoeffs[i][1]);
		if(iss >> value)
			outfile << value;
		else
			outfile << variables[potentials.angleCoeffs[i][1][0] - 65];

		linefactor++;
	}
	outfile << std::endl;
	*/

	outfile << "%FLAG AMOEBA_TRIGONAL_ANGLE_FTAB_DEGREE" << std::endl;
	outfile << "%FORMAT(I8)" << std::endl;
	outfile << "       6" << std::endl;
	outfile << "%FLAG AMOEBA_TRIGONAL_ANGLE_FTAB_COEFFS" << std::endl;
	outfile << "%FORMAT(5E16.8)" << std::endl;
	outfile << "  0.00000000e+00  0.00000000e+00  1.00000000e+00 -1.40000000e-02  5.60000000e-05" << std::endl;
	outfile << " -7.00000000e-07  2.20000000e-08" << std::endl;
	outfile << "%FLAG AMOEBA_OPBEND_ANGLE_NUM_LIST" << std::endl;
	outfile << "%FORMAT(I8)" << std::endl;
	outfile << "       0" << std::endl;
	outfile << "%FLAG AMOEBA_OPBEND_ANGLE_LIST" << std::endl;
	outfile << "%COMMENT dimension = (5,0)" << std::endl;
	outfile << "%FORMAT(10I8)" << std::endl;
	outfile << "" << std::endl;
	outfile << "%FLAG AMOEBA_OPBEND_ANGLE_NUM_PARAMS" << std::endl;
	outfile << "%FORMAT(I8)" << std::endl;
	outfile << "       0" << std::endl;
	outfile << "%FLAG AMOEBA_OPBEND_ANGLE_FORCE_CONSTANT" << std::endl;
	outfile << "%FORMAT(5E16.8)" << std::endl;
	outfile << "" << std::endl;
	outfile << "%FLAG AMOEBA_OPBEND_ANGLE_FTAB_DEGREE" << std::endl;
	outfile << "%FORMAT(I8)" << std::endl;
	outfile << "       6" << std::endl;
	outfile << "%FLAG AMOEBA_OPBEND_ANGLE_FTAB_COEFFS" << std::endl;
	outfile << "%FORMAT(5E16.8)" << std::endl;
	outfile << "  0.00000000e+00  0.00000000e+00  1.00000000e+00 -1.40000000e-02  5.60000000e-05" << std::endl;
	outfile << " -7.00000000e-07  2.20000000e-08" << std::endl;
	outfile << "%FLAG AMOEBA_TORSION_NUM_LIST" << std::endl;
	outfile << "%FORMAT(I8)" << std::endl;
	outfile << "       0" << std::endl;
	outfile << "%FLAG AMOEBA_TORSION_LIST" << std::endl;
	outfile << "%COMMENT   dimension = (5,0)" << std::endl;
	outfile << "%FORMAT(10I8)" << std::endl;
	outfile << "" << std::endl;
	outfile << "%FLAG AMOEBA_TORSION_NUM_PARAMS" << std::endl;
	outfile << "%FORMAT(I8)" << std::endl;
	outfile << "       0" << std::endl;
	outfile << "%FLAG AMOEBA_TORSION_FORCE_CONSTANT" << std::endl;
	outfile << "%FORMAT(5E16.8)" << std::endl;
	outfile << "" << std::endl;
	outfile << "%FLAG AMOEBA_TORSION_PERIODICITY" << std::endl;
	outfile << "%FORMAT(5E16.8)" << std::endl;
	outfile << "" << std::endl;
	outfile << "%FLAG AMOEBA_TORSION_PHASE" << std::endl;
	outfile << "%FORMAT(5E16.8)" << std::endl;
	outfile << "" << std::endl;
	outfile << "%FLAG AMOEBA_PI_TORSION_NUM_LIST" << std::endl;
	outfile << "%FORMAT(I8)" << std::endl;
	outfile << "       0" << std::endl;
	outfile << "%FLAG AMOEBA_PI_TORSION_LIST" << std::endl;
	outfile << "%COMMENT   dimension = (7,0)" << std::endl;
	outfile << "%FORMAT(10I8)" << std::endl;
	outfile << "" << std::endl;
	outfile << "%FLAG AMOEBA_PI_TORSION_NUM_PARAMS" << std::endl;
	outfile << "%FORMAT(I8)" << std::endl;
	outfile << "       0" << std::endl;
	outfile << "%FLAG AMOEBA_PI_TORSION_FORCE_CONSTANT" << std::endl;
	outfile << "%FORMAT(5E16.8)" << std::endl;
	outfile << "" << std::endl;
	outfile << "%FLAG AMOEBA_PI_TORSION_PERIODICITY" << std::endl;
	outfile << "%FORMAT(5E16.8)" << std::endl;
	outfile << "" << std::endl;
	outfile << "%FLAG AMOEBA_PI_TORSION_PHASE" << std::endl;
	outfile << "%FORMAT(5E16.8)" << std::endl;
	outfile << "" << std::endl;
	outfile << "%FLAG AMOEBA_STRETCH_BEND_NUM_LIST" << std::endl;
	outfile << "%FORMAT(I8)" << std::endl;
	outfile << "       0" << std::endl;
	outfile << "%FLAG AMOEBA_STRETCH_BEND_LIST" << std::endl;
	outfile << "%COMMENT   dimension = (4,0)" << std::endl;
	outfile << "%FORMAT(10I8)" << std::endl;
	outfile << "" << std::endl;
	outfile << "%FLAG AMOEBA_STRETCH_BEND_NUM_PARAMS" << std::endl;
	outfile << "%FORMAT(I8)" << std::endl;
	outfile << "       0" << std::endl;
	outfile << "%FLAG AMOEBA_STRETCH_BEND_FORCE_CONSTANT" << std::endl;
	outfile << "%FORMAT(5E16.8)" << std::endl;
	outfile << "" << std::endl;
	outfile << "%FLAG AMOEBA_STRETCH_BEND_ANGLE_EQUIL_VALUE" << std::endl;
	outfile << "%FORMAT(5E16.8)" << std::endl;
	outfile << "" << std::endl;
	outfile << "%FLAG AMOEBA_STRETCH_BEND_BOND1_EQUIL_VALUE" << std::endl;
	outfile << "%FORMAT(5E16.8)" << std::endl;
	outfile << "" << std::endl;
	outfile << "%FLAG AMOEBA_STRETCH_BEND_BOND2_EQUIL_VALUE" << std::endl;
	outfile << "%FORMAT(5E16.8)" << std::endl;
	outfile << "" << std::endl;
	outfile << "%FLAG AMOEBA_TORSION_TORSION_NUM_LIST" << std::endl;
	outfile << "%FORMAT(I8)" << std::endl;
	outfile << "0" << std::endl;
	outfile << "%FLAG AMOEBA_TORSION_TORSION_LIST" << std::endl;
	outfile << "%COMMENT   dimension = (6,147)" << std::endl;
	outfile << "%FORMAT(10I8)" << std::endl;
	outfile << "" << std::endl;
	outfile << "%FLAG AMOEBA_TORSION_TORSION_NUM_PARAMS" << std::endl;
	outfile << "%FORMAT(I8)" << std::endl;
	outfile << "       0" << std::endl;
	
	outfile << "%FLAG AMOEBA_VDW_ATOM_TYPES_NUM_LIST" << std::endl;
	outfile << "%FORMAT(I8)" << std::endl;
	outfile.width(8);
	outfile << geom.atoms.size() << std::endl;
	
	outfile << "%FLAG AMOEBA_VDW_ATOM_TYPES_LIST" << std::endl;
	outfile << "%COMMENT   dimension = (1,4)" << std::endl;
	outfile << "%FORMAT(10I8)";
	linefactor = 0;
	for(std::vector<Atom>::iterator atomit = geom.atoms.begin(); atomit != geom.atoms.end(); atomit++){
		int typeNum = 1;
		for(std::vector<std::string>::iterator it = types.begin(); it != types.end(); it++){
			if(atomit->symbol == *it){
				if(linefactor % 10 == 0)
					outfile << std::endl;
				outfile.width(8);
				outfile << typeNum;
				linefactor++;
			}
			typeNum++;
		}
	}
	outfile << std::endl;
	
	outfile << "%FLAG AMOEBA_VDW_ATOM_PARENT_LIST" << std::endl;
	outfile << "%COMMENT   dimension = (1,4)" << std::endl;
	outfile << "%FORMAT(10I8)";
	std::vector<int> oxyList;
	atomNum = 1;
	for(std::vector<Atom>::iterator atomit = geom.atoms.begin(); atomit != geom.atoms.end(); atomit++){
		if(atomit->symbol == "O") oxyList.push_back(atomNum);
		atomNum++;
	}

	int hCnt = 0;
	atomNum = 1;
	linefactor = 0;
	for(std::vector<Atom>::iterator atomit = geom.atoms.begin(); atomit != geom.atoms.end(); atomit++){
		if(linefactor % 10 == 0)
			outfile << std::endl;
		outfile.width(8);
		if(atomit->symbol == "H"){
			outfile << oxyList[hCnt/2];
			hCnt++;  // hack to fix the molecule number bug
		} else {
			outfile << atomNum;
		}
		atomNum++;
		linefactor++;
	}
	outfile << std::endl;
	
	outfile << "%FLAG AMOEBA_VDW_PARENT_COORD_WEIGHT_LIST" << std::endl;
	outfile << "%COMMENT   dimension = (1,4)" << std::endl;
	outfile << "%FORMAT(5E16.8)";
	linefactor = 0;
	for(std::vector<Atom>::iterator atomit = geom.atoms.begin(); atomit != geom.atoms.end(); atomit++){
			if(linefactor % 5 == 0)
				outfile << std::endl;
			outfile.width(16);
			outfile.precision(8);
			if(atomit->symbol == "H")
				outfile << 0.09;
			else
				outfile << 1.0;
			linefactor++;
	}
	outfile << std::endl;

	outfile << "%FLAG AMOEBA_VDW_BUFFER_DELTA" << std::endl;
	outfile << "%FORMAT(E16.8)" << std::endl;
	outfile << "  0.70000000E-01" << std::endl;
	outfile << "%FLAG AMOEBA_VDW_BUFFER_GAMMA" << std::endl;
	outfile << "%FORMAT(E16.8)" << std::endl;
	outfile << "  0.12000000E+00" << std::endl;

	outfile << "%FLAG AMOEBA_VDW_PARAMS_NUM_LIST" << std::endl;
	outfile << "%FORMAT(I8)" << std::endl;
	outfile << "       3" << std::endl;

	/*	
	float radii[3][3] = {{0.0,0.0,0.0},{0.0,0.0,0.0},{0.0,0.0,0.0}};
	float epsilons[3][3] = {{0.0,0.0,0.0},{0.0,0.0,0.0},{0.0,0.0,0.0}};
	int varindex;
	for(i = 0; i < potentials.pairAtoms.size(); i++){
	  int x, y, coeffNum;

	  coeffNum=potentials.pairCoeffNums[i];  // get the pair potential number
	  if ( coeffNum < 0 || coeffNum > potentials.pairCoeffType.size() ) {
	    MolecularDynamicsException excep("Coefficient number for atom types in pair potentials is out of range.");
	    throw excep;
	  } // Trap invalid coefficient numbers
	  
	  std::istringstream issx(potentials.pairAtoms[coeffNum][0]), issy(potentials.pairAtoms[coeffNum][1]);
	  if(!(issx >> x))
	    x = std::find(types.begin(), types.end(), potentials.pairAtoms[coeffNum][0]) - types.begin();
	  if(!(issy >> y))
	    y = std::find(types.begin(), types.end(), potentials.pairAtoms[coeffNum][1]) - types.begin();
	  x--;
	  y--;
	  std::istringstream issr(potentials.pairCoeffs[coeffNum][0]);
	  if(!(issr >> radii[x][y])) {
	    varindex=potentials.pairCoeffs[coeffNum][0][0]-65;
	    if ( varindex < 0 || varindex > 1 ) {
	      MolecularDynamicsException excep("Variable index is not 0 or 1 for radii. Fit cannot continue");
	      throw excep;
	    }
	    radii[x][y] = variables[varindex];
	  }
	  if(radii[x][y] < 0.00001){
	    MolecularDynamicsException excep("Detected zero value in 14-7 radius parameter. Fit cannot continue");
	    throw excep;
	  }
	  radii[y][x] = radii[x][y];
	  std::istringstream isse(potentials.pairCoeffs[coeffNum][1]);		
	  if(!(isse >> epsilons[x][y])) {
	    varindex=potentials.pairCoeffs[coeffNum][1][0]-65;
	    if ( varindex < 0 || varindex > 1 ) {
	      MolecularDynamicsException excep("Variable index is not 0 or 1 for epsilon. Fit cannot continue");
	      throw excep;
	    }		  
	    epsilons[x][y] = variables[varindex];
	    epsilons[y][x] = epsilons[x][y];
	  }

	}

	outfile << "%FLAG AMOEBA_VDW_MIXED_RADII_LIST" << std::endl;
	outfile << "%COMMENT   dimension = (3,3)" << std::endl;
	outfile << "%FORMAT(5E16.8)";
	linefactor = 0;
	for(i = 0; i < 3; i++){
		for(int j = 0; j < 3; j++){
			if(linefactor % 5 == 0)
				outfile << std::endl;
			outfile.width(16);
			outfile.precision(8);
			outfile << radii[i][j];
			linefactor++;
		}
	}
	outfile << std::endl;

	outfile << "%FLAG AMOEBA_VDW_MIXED_EPSILONS_LIST" << std::endl;
	outfile << "%COMMENT   dimension = (3,3)" << std::endl;
	outfile << "%FORMAT(5E16.8)";
	linefactor = 0;
	for(i = 0; i < 3; i++){
		for(int j = 0; j < 3; j++){
			if(linefactor % 5 == 0)
				outfile << std::endl;
			outfile.width(16);
			outfile.precision(8);
			outfile << epsilons[i][j];
			linefactor++;
		}
	}
	outfile << std::endl;
	*/
	
	outfile << "%FLAG AMOEBA_LOCAL_FRAME_MULTIPOLES_NUM_LIST" << std::endl;
	outfile << "%FORMAT(I8)" << std::endl;
	outfile.width(8);
	outfile << geom.atoms.size() << std::endl;
	
	outfile << "%FLAG AMOEBA_LOCAL_FRAME_MULTIPOLES_LIST" << std::endl;
	outfile << "%COMMENT dimension = (10,4)" << std::endl;
	outfile << "%FORMAT(5E16.8)" << std::endl;
	for(std::vector<Atom>::iterator atomit = geom.atoms.begin(); atomit != geom.atoms.end(); atomit++){
		if(atomit->symbol == "H"){
			outfile << "  2.59830000e-01 -3.85900000e-02  0.00000000e+00 -5.81800000e-02 -1.83650000e-02" << std::endl;
			outfile << " -5.36950000e-02  7.20600000e-02  0.00000000e+00 -2.03000000e-03  0.00000000e+00" << std::endl;
		} else if(atomit->symbol == "O"){
			outfile << " -5.19660000e-01  0.00000000e+00  0.00000000e+00  1.42790000e-01  1.89640000e-01" << std::endl;
			outfile << " -2.09045000e-01  1.94050000e-02  0.00000000e+00  0.00000000e+00  0.00000000e+00" << std::endl;
		} else {
			outfile.width(16);
			outfile.precision(8);
			outfile << metalCharge;
			outfile << "  0.00000000e+00  0.00000000e+00  0.00000000e+00  0.00000000e+00" << std::endl;
			outfile << "  0.00000000e+00  0.00000000e+00  0.00000000e+00  0.00000000e+00  0.00000000e+00" << std::endl;
		}
	}	
	
	outfile << "%FLAG AMOEBA_CHIRAL_FRAME_NUM_LIST" << std::endl;
	outfile << "%FORMAT(I8)" << std::endl;
	outfile << "       0" << std::endl;
	outfile << "%FLAG AMOEBA_CHIRAL_FRAME_LIST" << std::endl;
	outfile << "%COMMENT   dimension = (3,0)" << std::endl;
	outfile << "%FORMAT(10I8)" << std::endl;
	outfile << "" << std::endl;
	outfile << "%FLAG AMOEBA_FRAME_DEF_NUM_LIST" << std::endl;
	outfile << "%FORMAT(I8)" << std::endl;
	//	outfile << "       7" << std::endl; //
	outfile <<  7*((geom.atoms.size()-1)/3); // This entry appears to the number of water molecules * 7
	outfile << "" << std::endl;

	outfile << "%FLAG AMOEBA_FRAME_DEF_LIST" << std::endl;
	outfile << "%COMMENT   dimension = (5,x)" << std::endl;
	outfile << "%FORMAT(10I8)";
	int h1 = 0, h2 = 0, o1 = 0;
	atomNum = 1;
	linefactor = 0;
	for(std::vector<Atom>::iterator atomit = geom.atoms.begin(); atomit != geom.atoms.end(); atomit++){
		if(atomit->symbol == "H") {
			if(h1 == 0) h1 = atomNum;
			else if(h2 == 0) h2 = atomNum;
			else std::cout << "Warning: File is not in amoeba format. 3" << std::endl;
		} else if(atomit->symbol == "O") {
			if(o1 == 0) o1 = atomNum;
			else std::cout << "Warning: File is not in amoeba format. 4" << std::endl;
		}

		if(h1 != 0 && h2 != 0 && o1 != 0) {
			if(linefactor % 2 == 0)
				outfile << std::endl;
			outfile << std::setw(8) << o1 << std::setw(8) << 1 << std::setw(8) << o1 << std::setw(8) << h1 << std::setw(8) << 2;
			linefactor++;

			if(linefactor % 2 == 0)
				outfile << std::endl;
			outfile << std::setw(8) << o1 << std::setw(8) << 1 << std::setw(8) << o1 << std::setw(8) << h2 << std::setw(8) << 2;
			linefactor++;

			if(linefactor % 2 == 0)
				outfile << std::endl;
			outfile << std::setw(8) << o1 << std::setw(8) << 2 << std::setw(8) << o1 << std::setw(8) << h2 << std::setw(8) << 1;
			linefactor++;

			if(linefactor % 2 == 0)
				outfile << std::endl;
			outfile << std::setw(8) << h1 << std::setw(8) << 1 << std::setw(8) << h1 << std::setw(8) << o1 << std::setw(8) << 1;
			linefactor++;

			if(linefactor % 2 == 0)
				outfile << std::endl;
			outfile << std::setw(8) << h1 << std::setw(8) << 2 << std::setw(8) << h1 << std::setw(8) << h2 << std::setw(8) << 1;
			linefactor++;

			if(linefactor % 2 == 0)
				outfile << std::endl;
			outfile << std::setw(8) << h2 << std::setw(8) << 1 << std::setw(8) << h2 << std::setw(8) << o1 << std::setw(8) << 1;
			linefactor++;

			if(linefactor % 2 == 0)
				outfile << std::endl;
			outfile << std::setw(8) << h2 << std::setw(8) << 2 << std::setw(8) << h2 << std::setw(8) << h1 << std::setw(8) << 1;
			linefactor++;

			h1 = 0;
			h2 = 0;
			o1 = 0;
		}
		atomNum++;
	}
	//	if(linefactor % 2 == 1) //
	outfile << std::endl;

	outfile << "%FLAG AMOEBA_ADJUST_NUM_LIST" << std::endl;
	outfile << "%FORMAT(I8)" << std::endl;
	outfile.width(8);
	outfile << geom.atoms.size()-1 << std::endl;

	outfile << "%FLAG AMOEBA_ADJUST_LIST" << std::endl;
	outfile << "%COMMENT dimension = (3,3)" << std::endl;
	outfile << "%FORMAT(10I8)";
	h1 = 0, h2 = 0, o1 = 0;
	atomNum = 1;
	linefactor = 0;
	for(std::vector<Atom>::iterator atomit = geom.atoms.begin(); atomit != geom.atoms.end(); atomit++){
		if(atomit->symbol == "H") {
			if(h1 == 0) h1 = atomNum;
			else if(h2 == 0) h2 = atomNum;
			else std::cout << "Warning: File is not in amoeba format. 5" << std::endl;
		} else if(atomit->symbol == "O") {
			if(o1 == 0) o1 = atomNum;
			else std::cout << "Warning: File is not in amoeba format. 6" << std::endl;
		}

		if(h1 != 0 && h2 != 0 && o1 != 0) {
			if(linefactor % 10 == 0)
				outfile << std::endl;
			outfile.width(8);
			outfile << o1;
			linefactor++;

			if(linefactor % 10 == 0)
				outfile << std::endl;
			outfile.width(8);
			outfile << h1;
			linefactor++;

			if(linefactor % 10 == 0)
				outfile << std::endl;
			outfile.width(8);
			outfile << 1;
			linefactor++;

			if(linefactor % 10 == 0)
				outfile << std::endl;
			outfile.width(8);
			outfile << o1;
			linefactor++;

			if(linefactor % 10 == 0)
				outfile << std::endl;
			outfile.width(8);
			outfile << h2;
			linefactor++;

			if(linefactor % 10 == 0)
				outfile << std::endl;
			outfile.width(8);
			outfile << 1;
			linefactor++;

			if(linefactor % 10 == 0)
				outfile << std::endl;
			outfile.width(8);
			outfile << h1;
			linefactor++;

			if(linefactor % 10 == 0)
				outfile << std::endl;
			outfile.width(8);
			outfile << h2;
			linefactor++;

			if(linefactor % 10 == 0)
				outfile << std::endl;
			outfile.width(8);
			outfile << 2;
			linefactor++;

			o1 = h1 = h2 = 0;
		}
		atomNum++;
	}
	outfile << std::endl;

	outfile << "%FLAG AMOEBA_ADJUST_VDW_WEIGHTS_LIST" << std::endl;
	outfile << "%COMMENT   dimension = (9)" << std::endl;
	outfile << "%FORMAT(5E16.8)" << std::endl;
	outfile << "  0.00000000e+00  0.00000000e+00  1.00000000e+00  1.00000000e+00  0.00000000e+00" << std::endl;
	outfile << "  0.00000000e+00  1.00000000e+00  1.00000000e+00  1.00000000e+00" << std::endl;
	outfile << "%FLAG AMOEBA_ADJUST_MPOLE_WEIGHTS_LIST" << std::endl;
	outfile << "%COMMENT   dimension = (9)" << std::endl;
	outfile << "%FORMAT(5E16.8)" << std::endl;
	outfile << "  0.00000000e+00  0.00000000e+00  4.00000000e-01  8.00000000e-01  0.00000000e+00" << std::endl;
	outfile << "  0.00000000e+00  4.00000000e-01  8.00000000e-01  1.00000000e+00" << std::endl;
	outfile << "%FLAG AMOEBA_ADJUST_DIRECT_WEIGHTS_LIST" << std::endl;
	outfile << "%COMMENT   dimension = (9)" << std::endl;
	outfile << "%FORMAT(5E16.8)" << std::endl;
	outfile << "  1.00000000e+00  1.00000000e+00  1.00000000e+00  1.00000000e+00  0.00000000e+00" << std::endl;
	outfile << "  0.00000000e+00  0.00000000e+00  0.00000000e+00  0.00000000e+00" << std::endl;
	outfile << "%FLAG AMOEBA_ADJUST_POLAR_WEIGHTS_LIST" << std::endl;
	outfile << "%COMMENT   dimension = (9)" << std::endl;
	outfile << "%FORMAT(5E16.8)" << std::endl;
	outfile << "  0.00000000e+00  0.00000000e+00  1.00000000e+00  1.00000000e+00  0.00000000e+00" << std::endl;
	outfile << "  0.00000000e+00  5.00000000e-01  1.00000000e+00  1.00000000e+00" << std::endl;
	outfile << "%FLAG AMOEBA_ADJUST_MUTUAL_WEIGHTS_LIST" << std::endl;
	outfile << "%COMMENT   dimension = (9)" << std::endl;
	outfile << "%FORMAT(5E16.8)" << std::endl;
	outfile << "  1.00000000e+00  1.00000000e+00  1.00000000e+00  1.00000000e+00  1.00000000e+00" << std::endl;
	outfile << "  1.00000000e+00  1.00000000e+00  1.00000000e+00  1.00000000e+00" << std::endl;
	outfile << "%FLAG AMOEBA_POLARIZABILITY_NUM_LIST" << std::endl;
	outfile << "%FORMAT(I8)" << std::endl;
	outfile.width(8);
	outfile << geom.atoms.size() << std::endl;
	outfile << "%FLAG AMOEBA_POLARIZABILITY_LIST" << std::endl;
	outfile << "%COMMENT   dimension = (4)" << std::endl;
	outfile << "%FORMAT(5E16.8)";
	linefactor = 0;
	for(std::vector<Atom>::iterator atomit = geom.atoms.begin(); atomit != geom.atoms.end(); atomit++){
			if(linefactor % 5 == 0)
				outfile << std::endl;
			outfile.width(16);
			outfile.precision(8);
			if(atomit->symbol == "H")
				outfile << 0.496;
			else if(atomit->symbol == "O")
				outfile << 0.837;
			else
				outfile << polarizability;
			linefactor++;
	}
	outfile << std::endl;
	
	outfile.close();

}

void Amber::generateInpcrd(Geometry geom, std::string directory){
	mkdir(directory.c_str(),0700);
	std::ofstream outfile((directory + "/inpcrd").c_str());
	outfile << "newffm Generated Geometry" << std::endl;
	outfile.width(6);
	outfile << geom.atoms.size() << std::endl;
	
	outfile.precision(7);
	outfile.setf(std::ios_base::fixed,std::ios::floatfield);
	int origAtomNum = 1;
	for(std::vector<Atom>::iterator atomit = geom.atoms.begin(); atomit != geom.atoms.end(); atomit++){
		outfile << "  ";
		outfile.width(10);
		outfile << atomit->x;
		outfile << "  ";
		outfile.width(10);
		outfile << atomit->y;
		outfile << "  ";
		outfile.width(10);
		outfile << atomit->z;

		if(origAtomNum % 2 == 0){
			outfile << std::endl;
		}
		origAtomNum++;
	}
	if(origAtomNum % 2 != 1)
		outfile << std::endl;
	outfile << "  ";
	outfile.width(10);
	outfile << 40.0;
	outfile << "  ";
	outfile.width(10);
	outfile << 40.0;
	outfile << "  ";
	outfile.width(10);
	outfile << 40.0;
	outfile << "  ";
	outfile.width(10);
	outfile << 90.0;
	outfile << "  ";
	outfile.width(10);
	outfile << 90.0;
	outfile << "  ";
	outfile.width(10);
	outfile << 90.0 << std::endl;
	outfile.close();
}

void Amber::readForces(std::string directory){
	std::ifstream inFile((directory + "/mdout").c_str());
	std::string line;
	float energy, temp, x, y, z;
	int step, atom, i, totatoms=0;
	bool finished = false;
	std::vector<float> atomForce;
	float number;
	
	while(std::getline(inFile, line)){
		if(line.find("FINAL RESULTS") == std::string::npos){		
			finished = true;
			break;
		}
	}

	if(!finished){
		MolecularDynamicsException excep("Sander failed. Check amberRun/ for errors and try again.");
		throw excep;
	}

	while(std::getline(inFile, line) && line.find("NSTEP") == std::string::npos);

	std::getline(inFile, line);
	std::istringstream iss(line);
	iss >> step >> energy;
	energies.push_back(energy);

	std::cout << "    Energy: " << energy << std::endl;

	if(step != 1){
		MolecularDynamicsException excep("Couldn't find step 1 in mdout. Check amberRun/ for errors and try again.");
		throw excep;
	}

	inFile.close();

	std::ifstream velFile((directory + "/mdvel").c_str());

	//Skip title
	std::getline(velFile, line);

	std::vector< std::vector<float> > geom;
	while(std::getline(velFile, line)){
		std::istringstream iss2(line);
		int atomnum=0;
		while(iss2 >> number ){
		    atomForce.push_back(number);
		    if (atomForce.size() == 3 ) {
			geom.push_back(atomForce);
			atomnum++;
			std::cout << "Forces:   Atom " << atomnum << " " << atomForce[0] << " " << atomForce[1] << " " << atomForce[2] << std::endl;
			atomForce.clear();
		    }
		}
	}
	if(!atomForce.empty()){
		MolecularDynamicsException excep("No forces were read. Check amberRun/ for errors and try again.");
		throw excep;
	}
	forces.push_back(geom);
	velFile.close();
}

std::vector<struct mdQuestion> Amber::getQuestions(){
	std::vector<struct mdQuestion> result;

	struct mdQuestion qMetalCharge;
	qMetalCharge.phrase = "Metal Charge";
	qMetalCharge.type = MDSTRING;
	qMetalCharge.stringAnswer = "3.0";
	
	result.push_back(qMetalCharge);

	struct mdQuestion qPolarizability;
	qPolarizability.phrase = "Metal Polarizability";
	qPolarizability.type = MDSTRING;
	qPolarizability.stringAnswer = "1.134";

	result.push_back(qPolarizability);

	struct mdQuestion qSanderExe;

	qSanderExe.phrase = "Sander executable";
	qSanderExe.type = MDFILENAME;

	result.push_back(qSanderExe);

	return result;
}

void Amber::prepareMD(std::vector<GeometrySet> geoSets, std::vector<float> variables, std::vector<struct mdQuestion> questions){
	this->geoSets = geoSets;
	this->variables = variables;

	metalCharge = atof(questions[0].stringAnswer.c_str());
	polarizability = atof(questions[1].stringAnswer.c_str());
	amberPath = questions[2].stringAnswer;

	mkdir("amberRun",0700);
}

void Amber::runMD(){
	int i;
	std::string line;
	std::ofstream outfile;
	
	resetRun();

	i++;

	int set = 1;
	for(std::vector<GeometrySet>::iterator setit = geoSets.begin(); setit != geoSets.end(); setit++){
		i = 1;
		std::string filename;
		for(std::vector<Geometry>::iterator geomit = setit->geometries.begin(); geomit != setit->geometries.end(); geomit++){
			int atomNum = 1;
			std::ostringstream oss(std::ostringstream::out);
			oss << "amberRun/set_" << set << "_geometry_" << i;
			generateInpcrd(*geomit, oss.str());
			i++;
		}
	
		i = 1;
		for(std::vector<Geometry>::iterator geomit = setit->geometries.begin(); geomit != setit->geometries.end(); geomit++){
			std::ifstream gInControl("amberRun/prmtop");
			std::ostringstream oss(std::ostringstream::out);
			
			oss << "amberRun/set_" << set << "_geometry_" << i;

			generateMdin(oss.str());
			generatePrmtop(*geomit, oss.str());
			
			i++;
		}
		set++;
	}

	struct stat fileInfo;
	if(stat(amberPath.c_str(), &fileInfo)){
		std::cout << "Error: Couldn't stat " << amberPath << std::endl;
		return;
	}

	std::cout << "Running sander..." << std::endl;
	
	set = 1;
	for(std::vector<GeometrySet>::iterator setit = geoSets.begin(); setit != geoSets.end(); setit++){
		i = 1;
		for(std::vector<Geometry>::iterator geomit = setit->geometries.begin(); geomit != setit->geometries.end(); geomit++){
			std::ostringstream oss(std::ostringstream::out);
			
			oss << "amberRun/set_" << set << "_geometry_" << i;
			chdir(oss.str().c_str());
			// Removes moved here to allow the files to remain after a crash
			remove("mdcrd");
			remove("mdinfo");
			remove("mdout");
			remove("mdvel");
			remove("restrt");
			system(amberPath.c_str());

			try {
				readForces(".");
			} catch(const MolecularDynamicsException & e){
				chdir("..");
				chdir("..");
				throw e;
			}
			// Removes were here originally
			
			chdir("..");
			chdir("..");
			i++;
		}
	}

	std::cout << "Sander complete." << std::endl;
}
