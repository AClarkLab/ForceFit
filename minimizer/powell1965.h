#ifndef _POWELL1965_H
#define _POWELL1965_H

#include <fstream>

#include "minimizer.h"
#include "../molecularDynamics/molecularDynamics.h"

class Powell1965 : public Minimizer {
public:
	virtual std::vector<struct minQuestion> getQuestions();
	virtual void prepareMin(const std::vector<GeometrySet> & geoSets, std::vector<struct minQuestion> questions);
	virtual void runMin(MolecularDynamics * md);
	static std::string const name;
	virtual std::string const getName();
private:
	std::vector<float> calculateUForce(struct dataPoint data);
	std::vector<float> calculateUEnergy(struct dataPoint data);
	struct dataPoint runValues(MolecularDynamics * md, std::vector<float> values);

	MolecularDynamics * md;
	
	float linmin(std::vector<float> & p, std::vector<float> & xi);
	float fidim(std::vector<float> & p, std::vector<float> & xi, float x);
	float brent(float ax, float bx, float cx, std::vector<float> & p, std::vector<float> & xi, float tol, float *xmin);
	void mnbrak(float * ax, float * bx, float * cx, float * fa, float * fb, float *fc, std::vector<float> & p, std::vector<float> & xi);

	float epsilon;
	bool useForce, useEnergy;
	std::vector<float> finalVariables;

	std::ofstream logfile;
};

#endif
