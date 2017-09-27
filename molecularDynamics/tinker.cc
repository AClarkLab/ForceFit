#include <iostream>
#include <sstream>
#include <fstream>
#include <sys/stat.h>
#include <cstdlib>

#include "tinker.h"
	
std::string const Tinker::name = "Tinker";

std::string const Tinker::getName(){
	return "Tinker";
}

std::vector<struct mdQuestion> Tinker::getQuestions(){
	std::vector<struct mdQuestion> result;

	struct mdQuestion qXyzFileTemplate;
	qXyzFileTemplate.phrase = "XYZ File Template";
	qXyzFileTemplate.type = MDFILENAME;
	
	result.push_back(qXyzFileTemplate);

	struct mdQuestion qPrmFile;
	qPrmFile.phrase = "PRM File";
	qPrmFile.type = MDFILENAME;
	
	result.push_back(qPrmFile);

	struct mdQuestion qKeyFile;
	qKeyFile.phrase = "Key File";
	qKeyFile.type = MDFILENAME;
	
	result.push_back(qKeyFile);

	struct mdQuestion qTinkerExe;
	qTinkerExe.phrase = "Tinker executable";
	qTinkerExe.type = MDFILENAME;

	result.push_back(qTinkerExe);

	return result;
}

void Tinker::generateXyz(std::string directory, Geometry geom){
	std::ifstream inXyzFile(xyzFile.c_str());
	std::ofstream outXyzFile((directory + "/ForceFit.xyz").c_str());
	std::string line;

	//Copy first line
	std::getline(inXyzFile, line);
	outXyzFile << line << std::endl;

	int i = 0;
	while(std::getline(inXyzFile, line)){
		if(i >= geom.atoms.size()){
			MolecularDynamicsException excep("Tinker XYZ file has incorrect number of atoms.");
			throw excep;
		}

		std::istringstream iss(line);
		int num;
		std::string symbol;
		float coordinate;
		int type;
		int bond;

		iss >> num >> symbol >> coordinate >> coordinate >> coordinate >> type;

		outXyzFile << " " << num << " " << symbol << " " << geom.atoms[i].x << " " << geom.atoms[i].y << " " << geom.atoms[i].z << " " << type;

		while(iss >> bond){
			outXyzFile << " " << bond;
		}
		outXyzFile << std::endl;
		i++;
	}
}

void Tinker::generatePrm(std::string directory, Geometry geom){
	std::ifstream inPrmFile(prmFile.c_str());
	std::ofstream outPrmFile((directory + "/ForceFit.prm").c_str());
	std::string line;

	while(std::getline(inPrmFile, line)) {
		for(std::vector<float>::iterator it = variables.begin(); it != variables.end(); it++){
			std::string needle = "#_#";
			needle[1] = (it - variables.begin()) + 65;
			std::string::size_type pos = 0;
			while((pos = line.find(needle, pos)) != std::string::npos){
				std::ostringstream variable(std::ostringstream::out);
				variable << (*it);
				line.replace(pos, needle.size(), variable.str());
			}
		}
		outPrmFile << line << std::endl;
	}

	inPrmFile.close();
	outPrmFile.close();
}

void Tinker::readForces(std::string directory){
	std::ifstream stdoutFile((directory + "/stdout").c_str());
	std::string line;
	std::string temp;
	std::string forcesFileName;
	
	std::getline(stdoutFile, line);
	
	while(line.find("E Total") == std::string::npos)
		std::getline(stdoutFile, line);
	
	std::getline(stdoutFile, line);
	std::getline(stdoutFile, line);
	
	std::istringstream stdout2Iss(line);
	float energy;

	stdout2Iss >> temp >> energy;

	energies.push_back(energy);
	while(line.find("Force Vector File") == std::string::npos)
		std::getline(stdoutFile, line);

	std::istringstream stdoutIss(line);
	stdoutIss >> temp >> temp >> temp >> forcesFileName;

	stdoutFile.close();

	std::cout << (directory + "/" + forcesFileName) << std::endl;

	std::ifstream forcesFile((directory + "/" + forcesFileName).c_str());
	int n;
	std::string symbol;
	float x, y, z;

	std::getline(forcesFile, line);

	std::vector< std::vector<float> > geom;
	while(std::getline(forcesFile, line)){
		std::istringstream iss(line);
		iss >> n >> symbol >> x >> y >> z;
		std::vector<float> forcesRow;
		forcesRow.push_back(x);
		forcesRow.push_back(y);
		forcesRow.push_back(z);

		geom.push_back(forcesRow);
	}
	forces.push_back(geom);

	forcesFile.close();
}

void Tinker::prepareMD(std::vector<GeometrySet> geoSets, std::vector<float> variables, std::vector<struct mdQuestion> questions){
	this->geoSets = geoSets;
	this->variables = variables;

	xyzFile = questions[0].stringAnswer;
	prmFile = questions[1].stringAnswer;
	keyFile = questions[2].stringAnswer;
	tinkerExe = questions[3].stringAnswer;

	mkdir("tinkerRun",0700);
}

void Tinker::runMD(){
	std::string line;
	resetRun();
	
	int set = 1;
	for(std::vector<GeometrySet>::iterator setit = geoSets.begin(); setit != geoSets.end(); setit++){
		int i = 1;
		for(std::vector<Geometry>::iterator geomit = setit->geometries.begin(); geomit != setit->geometries.end(); geomit++){
			std::ostringstream oss(std::ostringstream::out);
			oss << "tinkerRun/set_" << set << "_geometry_" << i;
			
			{
				mkdir(oss.str().c_str(),0700);
			}

			{
				std::ifstream inKeyFile(keyFile.c_str());
				std::ostringstream oss2(std::ostringstream::out);
				
				oss2 << "tinkerRun/set_" << set << "_geometry_" << i << "/ForceFit.key";
				std::ofstream outKeyFile(oss2.str().c_str());
				while(std::getline(inKeyFile, line))
					outKeyFile << line << std::endl;
				inKeyFile.close();
				outKeyFile.close();
			}

			generateXyz(oss.str(), *geomit);

			generatePrm(oss.str(), *geomit);
			
			i++;
		}
		set++;
	}

	set = 1;
	for(std::vector<GeometrySet>::iterator setit = geoSets.begin(); setit != geoSets.end(); setit++){
		int i = 1;
		for(std::vector<Geometry>::iterator geomit = setit->geometries.begin(); geomit != setit->geometries.end(); geomit++){
			std::ostringstream oss(std::ostringstream::out);
			
			oss << "tinkerRun/set_" << set << "_geometry_" << i;
			chdir(oss.str().c_str());
			std::ostringstream execString;
			execString << tinkerExe << " ForceFit 1 1.0 0.001 1 >& stdout";
			system(execString.str().c_str());
			chdir("..");
			chdir("..");
			readForces(oss.str());

			i++;
		}
		set++;
	}

}
