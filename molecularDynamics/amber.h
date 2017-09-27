#ifndef _AMBER_H
#define _AMBER_H

#include "molecularDynamics.h"

class Amber : public MolecularDynamics {
public:
	Amber();
	virtual std::vector<struct mdQuestion> getQuestions();
	virtual void prepareMD(std::vector<GeometrySet> geoSets, std::vector<float> variables, std::vector<struct mdQuestion> questions);
	virtual void runMD();
	static std::string const name;
	virtual std::string const getName();
private:
	void generatePrmtop(Geometry geom, std::string directory);
	void generateMdin(std::string directory);
	void generateInpcrd(Geometry geom, std::string directory);
	void readForces(std::string);

	float metalCharge;
	float polarizability;
	std::string amberPath;
};

#endif
