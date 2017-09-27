#include <cmath>
#include <ctime>
#include <gmm/gmm.h>

#include "powell1965.h"
#include "../molecularDynamics/molecularDynamics.h"

#define TOL 2.0e-4
#define ITMAX 100
#define CGOLD 0.3819660
#define ZEPS 1.0e-10
#define GOLD 1.618034
#define GLIMIT 100.0
#define TINY 1.0e-20
#define MAX(a,b) ((a) > (b) ? (a) : (b))
#define SIGN(a,b) ((b) > 0.0 ? fabs(a) : -fabs(a))
#define SHFT(a,b,c,d) (a)=(b);(b)=(c);(c)=(d);

std::string const Powell1965::name = "Powell 1965";

std::string const Powell1965::getName(){
	return "Powell 1965";
}

std::vector<struct minQuestion> Powell1965::getQuestions(){
	std::vector<struct minQuestion> questions;

	struct minQuestion one;
	one.type = MINSTRING;
	one.phrase = "Step size for minimization in percent";
	one.stringAnswer = std::string("10");
	questions.push_back(one);

	struct minQuestion two;
	two.type = MINBOOL;
	two.phrase = "Use forces in force field fitting";
	two.boolAnswer = true;
	questions.push_back(two);

	struct minQuestion three;
	three.type = MINBOOL;
	three.phrase = "Use energies in force field fitting";
	three.boolAnswer = true;
	questions.push_back(three);

	return questions;
}

void Powell1965::prepareMin(const std::vector<GeometrySet> & geoSets, std::vector<struct minQuestion> questions){
	this->geoSets = geoSets;

	epsilon = atof(questions[0].stringAnswer.c_str())/100;
	useForce = questions[1].boolAnswer;
	useEnergy = questions[2].boolAnswer;
}

struct dataPoint {
	std::vector<std::vector<std::vector<float> > > forces;
	std::vector<float> energy;
	std::vector<float> values;
};

std::vector<float> Powell1965::calculateUForce(struct dataPoint data){
	std::vector<float> u;
	float q = 2.0*23.0605423, c = 1.0;

	std::vector<Geometry> geoms;
	for(std::vector<GeometrySet>::iterator setit = geoSets.begin(); setit != geoSets.end(); setit++){
		for(std::vector<Geometry>::iterator geomit = setit->geometries.begin(); geomit != setit->geometries.end(); geomit++){
			geoms.push_back(*geomit);
		}
	}
	
	if(geoms.size() != data.forces.size()){
		std::cout << "Incorrect number of geometries in forces from MD." << std::endl;
		MolecularDynamicsException excep("Incorrect number of geometries in forces from MD.");
		throw excep;
	}
	float rmsdx = 0.0, rmsdy = 0.0, rmsdz = 0.0;
	int count = 0;

	for(int geom = 0; geom < geoms.size(); geom++){
		if(geoms[geom].atoms.size() != data.forces[geom].size()){
			std::cout << "Incorrect number of atoms in forces from MD." << std::endl;
			MolecularDynamicsException excep("Incorrect number of atoms in forces from MD.");
			throw excep;
		}
		
		float contribution = 0.0;
		if(useForce){
			for(int atom = 0; atom < geoms[geom].atoms.size(); atom++){
				float w = c / (1 + (pow(geoms[geom].atoms[atom].forcex,2) + pow(geoms[geom].atoms[atom].forcey,2) + pow(geoms[geom].atoms[atom].forcez,2))/pow(q,2));

				contribution += ((sqrt(pow(data.forces[geom][atom][0],2)*w) - sqrt(pow(geoms[geom].atoms[atom].forcex,2)*w)) + 
						(sqrt(pow(data.forces[geom][atom][1],2)*w) - sqrt(pow(geoms[geom].atoms[atom].forcey,2)*w)) +
						(sqrt(pow(data.forces[geom][atom][2],2)*w) - sqrt(pow(geoms[geom].atoms[atom].forcez,2)*w)));

				rmsdx += pow(data.forces[geom][atom][0] - geoms[geom].atoms[atom].forcex,2);
				rmsdy += pow(data.forces[geom][atom][1] - geoms[geom].atoms[atom].forcey,2);
				rmsdz += pow(data.forces[geom][atom][2] - geoms[geom].atoms[atom].forcez,2);

				count++;
			}
		}
		u.push_back(contribution);
	}
	
	float rmsdtotal = sqrt((rmsdx + rmsdy + rmsdz)/count);

	rmsdx = sqrt(rmsdx/count);
	rmsdy = sqrt(rmsdy/count);
	rmsdz = sqrt(rmsdz/count);

	logfile << "    RMSD: " << rmsdtotal << std::endl;
	logfile << "    X: " << rmsdx << " Y: " << rmsdy << " Z: " << rmsdz << std::endl;
	logfile.flush();

	return u;
}

std::vector<float> Powell1965::calculateUEnergy(struct dataPoint data){
	std::vector<float> u;
	float q = 2.0*23.0605423, c = 1.0;

	std::vector<Geometry> geoms;
	for(std::vector<GeometrySet>::iterator setit = geoSets.begin(); setit != geoSets.end(); setit++){
		for(std::vector<Geometry>::iterator geomit = setit->geometries.begin(); geomit != setit->geometries.end(); geomit++){
			geoms.push_back(*geomit);
		}
	}
	
	if(geoms.size() != data.energy.size()){
		std::cout << "Incorrect number of atoms in energies from MD. QM: " << geoms.size() << " MD:" << data.energy.size() << std::endl;
		return u;
	}
	
	for(int geom = 0; geom < geoms.size(); geom++){
		if(useEnergy){
			float w = c / (1 + pow(geoms[geom].energy,2) / pow(q,2));

			u.push_back(sqrt(pow(data.energy[geom],2)*w) - sqrt(pow(geoms[geom].energy,2)*w));
		} else {
			u.push_back(0.0);
		}
	}

	return u;
}

float Powell1965::linmin(std::vector<float> & p, std::vector<float> & xi){
	float xx, xmin, fx, fb, fa, bx, ax, fret;
	std::vector<float> pcom, xicom;

	for(int i = 0; i < p.size(); i++) pcom.push_back(p[i]);
	for(int i = 0; i < xi.size(); i++) xicom.push_back(xi[i]);

	ax = 0.0;
	xx = 1.0;

	mnbrak(&ax, &xx, &bx, &fa, &fx, &fb, p, xi);
	fret = brent(ax,xx,bx,p, xi, TOL, &xmin);

	return fret;
}

float Powell1965::fidim(std::vector<float> & p, std::vector<float> & xi, float x){
	std::vector<float> xt(p.size());
	float f = 0.0;

	for(int j = 0; j < p.size(); j++){
		xt[j] = p[j] + x * xi[j];
	}
	
	struct dataPoint pos = runValues(md, xt);
	std::vector<float> uforce = calculateUForce(pos);
	std::vector<float> uenergy = calculateUEnergy(pos);

	for(int i = 0; i < uforce.size() || i < uenergy.size(); i++){
		if(useForce)
			f += uforce[i]*uforce[i];
		if(useEnergy)
			f += uenergy[i]*uenergy[i];
	}
	finalVariables = xt;
	return f;
}

float Powell1965::brent(float ax, float bx, float cx, std::vector<float> & p, std::vector<float> & xi, float tol, float *xmin){
	float a,b,d,etemp,fu,fv,fw,fx,o,q,r,tol1,tol2,u,v,w,x,xm;
	float e=0.0;

	a = ((ax < cx) ? ax : cx);
	b = ((ax > cx) ? ax : cx);
	x=w=v=bx;
	fw=fv=fx=fidim(p,xi,x);

	for(int iter = 1; iter <= ITMAX; iter++){
		xm = 0.5 * (a + b);
		tol2 = 2.0 * (tol1 = tol * fabs(x) + ZEPS);

		if(fabs(x - xm) <= (tol2 - 0.5 * (b - a))){
			*xmin = x;
			return fx;
		}

		if(fabs(e) > tol1){
			r = (x - w) * (fx - fv);
			q = (x - v) * (fx - fw);
			o = (x - v) * q - (x - w) * r;
			q = 2.0 * (q - r);
			
			if(q > 0.0)
				o = -o;

			q = fabs(q);
			etemp = e;
			e = d;

			if(fabs(o) >= fabs(0.5 * q * etemp) || o <= q * (a - x) || o >= q * (b - x)) {
				d = CGOLD * (e = (x >= xm ? a - x : b - x));
			} else {
				d = o / q;
				u = x + d;
				if(u - a < tol2 || b - u < tol2)
					d = SIGN(tol1, xm - x);
			}
		} else {
			d = CGOLD * (e = (x >= xm ? a - x : b -x ));
		}
		u = (fabs(d) >= tol1 ? x + d : x + SIGN(tol1, d));
		fu = fidim(p,xi,u);

		if(fu <= fx) {
			if(u >= x) a = x;
			else b = x;
			SHFT(v,w,x,u)
			SHFT(fv,fw,fx,fu)
		} else {
			if(u < x) a = u;
			else b = u;
			if(fu < fw || w == x){
				v = w;
				w = u;
				fv = fw;
				fw = fu;
			} else if(fu <= fv || v == x || v == w) {
				v = u;
				fv = fu;
			}
		}
	}
	*xmin = x;

	return fx;
}

void Powell1965::mnbrak(float * ax, float * bx, float * cx, float * fa, float * fb, float *fc, std::vector<float> & p, std::vector<float> & xi){
	float ulim, u, r, q, fu, dum;

	*fa = fidim(p,xi,*ax);
	*fb = fidim(p,xi,*bx);

	if(*fb > *fa){
		SHFT(dum, *ax, *bx, dum)
		SHFT(dum, *fb, *fa, dum)
	}

	*cx = (*bx)+GOLD*(*bx-*ax);
	*fc = fidim(p,xi,*cx);

	while(*fb > *fc){
		r = (*bx - *ax) * (*fb - *fc);
		q = (*bx - *cx) * (*fb - *fa);
		u = (*bx) - ((*bx - *cx) * q - (*bx - *ax) * r) / (2.0 * SIGN(MAX(fabs(q-r), TINY), q-r));
		ulim = (*bx) + GLIMIT * (*cx - *bx);

		if((*bx - u) * (u - *cx) > 0.0){
			fu = fidim(p,xi,u);
			
			if(fu < *fc){
				*ax = *bx;
				*bx = u;
				*fa = *fb;
				*fb = fu;
				
				return;
			} else if(fu > *fb) {
				*cx = u;
				*fc = fu;

				return;
			}

			u = (*cx) + GOLD * (*cx - *bx);
			fu = fidim(p,xi,u);
		} else if((*cx - u) * (u - ulim) > 0.0) {
			fu = fidim(p,xi,u);

			if(fu < *fc) {
				SHFT(*bx, *cx, u, *cx+GOLD*(*cx - *bx))
				SHFT(*fb, *fc, fu, fidim(p,xi,u))
			}
		} else if((u - ulim) * (ulim - *cx) >= 0.0) {
			u = ulim;
			fu = fidim(p,xi,u);
		} else {
			u = (*cx) + GOLD * (*cx - *bx);
			fu = fidim(p,xi,u);
		}
		SHFT(*ax, *bx, *cx, u)
			SHFT(*fa, *fb, *fc, fu);
	}
}

struct dataPoint Powell1965::runValues(MolecularDynamics * md, std::vector<float> values){
	struct dataPoint result;
	
	time_t now;
	logfile.open("powell.log");
	time(&now);
	logfile << "*** Running MD at " << ctime(&now);
	logfile << "    Starting variables: ";
	for(std::vector<float>::iterator it = values.begin(); it != values.end(); it++){
		logfile << (*it) << " ";
	}
	logfile << std::endl;

	md->setVariables(values);
	md->runMD();

	result.forces = md->getForces();
	logfile << "    Energies: ";
	for(std::vector<std::vector<std::vector<float> > >::iterator it = result.forces.begin(); it != result.forces.end(); it++){
		logfile << "    Geometry " << (it - result.forces.begin() + 1) << std::endl;
		
		for(std::vector<std::vector<float> >::iterator geomit = it->begin(); geomit != it->end(); geomit++){
			logfile << "       ";
			for(std::vector<float>::iterator valit = geomit->begin(); valit != geomit->end(); valit++){
				logfile << (*valit) << " ";
			}
			logfile << std::endl;
		}
	}
	logfile << std::endl;
	result.energy = md->getEnergies();
	logfile << "    Energies: ";
	for(std::vector<float>::iterator it = result.energy.begin(); it != result.energy.end(); it++){
		logfile << (*it) << " ";
	}
	logfile << std::endl;
	result.values = md->getVariables();
	logfile << "    Variables: ";
	for(std::vector<float>::iterator it = values.begin(); it != values.end(); it++){
		logfile << (*it) << " ";
	}
	logfile << std::endl;
	time(&now);
	logfile << "*** Ended MD at " << ctime(&now);

	return result;
}

void Powell1965::runMin(MolecularDynamics * md){
	this->md = md;
	//Center point and steps in each direction
	struct dataPoint center;
	std::vector<struct dataPoint> positive;

	std::vector<float> centerForce, centerEnergy;

	std::vector<float> variableValues = md->getVariables();
	
	center = runValues(md, variableValues);

	centerForce = calculateUForce(center);

	logfile << "    Center Forces Us:" << std::endl;
	for(int i = 0; i < centerForce.size(); i++){
		logfile << "\t" << centerForce[i] << std::endl;
	}
	
	centerEnergy = calculateUEnergy(center);
	
	logfile << "    Center Energy Us:" << std::endl;
	for(int i = 0; i < centerEnergy.size(); i++){
		logfile << "\t" << centerEnergy[i] << std::endl;
	}
	
	std::vector< std::vector<float> > gammaPosForce, gammaPosEnergy;
	std::vector< std::vector<float> > tPosForceMatrix, tPosEnergyMatrix;
	std::vector<float> sPosForce, sPosEnergy;
	std::vector<float> pForce, pEnergy;
	
	std::vector<float> newValues(center.values);
	for(int i = 0; i < newValues.size(); i++){
		newValues[i] *= 1.0 + epsilon;

		struct dataPoint pos = runValues(md, newValues);
		positive.push_back(pos);
		if(useForce)
			tPosForceMatrix.push_back(calculateUForce(pos));
		if(useEnergy)
			tPosEnergyMatrix.push_back(calculateUEnergy(pos));

		newValues[i] /= 1.0 + epsilon;
	}

	for(int i = 0; i < tPosForceMatrix.size(); i++){
		float sum = 0.0;

		for(int j = 0; j < tPosForceMatrix[i].size(); j++){
			sum += pow((tPosForceMatrix[i][j] - centerForce[j])/epsilon,2);
		}
		
		float s = (1/sqrt(sum));
		
		sPosForce.push_back(s);

		std::vector<float> gammaList;
		float pContrib = 0.0;
		for(int j = 0; j < tPosForceMatrix[i].size(); j++){
			float gamma = s * (tPosForceMatrix[i][j] - centerForce[j])/epsilon;
			
			gammaList.push_back(gamma);

			pContrib += gamma * centerForce[j];
		}
		gammaPosForce.push_back(gammaList);

		pForce.push_back(0.0-pContrib);
	}

	for(int i = 0; i < tPosEnergyMatrix.size(); i++){
		float sum = 0.0;

		for(int j = 0; j < tPosEnergyMatrix[i].size(); j++){
			sum += pow((tPosEnergyMatrix[i][j] - centerEnergy[j])/epsilon,2);
		}

		float s = (1/sqrt(sum));
		
		sPosEnergy.push_back(s);

		std::vector<float> gammaList;
		float pContrib = 0.0;
		for(int j = 0; j < tPosEnergyMatrix[i].size(); j++){
			float gamma = s * (tPosEnergyMatrix[i][j] - centerEnergy[j])/epsilon;
			
			gammaList.push_back(gamma);

			pContrib += gamma * centerEnergy[j];
		}
		gammaPosEnergy.push_back(gammaList);

		//Add to the p already created by the force
		pEnergy.push_back(0.0-pContrib);
	}

	gmm::dense_matrix<float> aForce(pForce.size(), pForce.size());
	gmm::dense_matrix<float> aEnergy(pEnergy.size(), pEnergy.size());

	for(int i = 0; i < pForce.size() || i < pEnergy.size(); i++){
		for(int j = 0; j < pForce.size() || j < pEnergy.size(); j++){
			float sum = 0.0;
			
			if(i < pForce.size() && j < pForce.size()){
				for(int k = 0; k < gammaPosForce[i].size(); k++)
					sum += gammaPosForce[i][k] * gammaPosForce[j][k];

				aForce(i, j) = sum;
			}
		
			if(i < pEnergy.size() && j < pEnergy.size()){
				sum = 0.0;
				for(int k = 0; k < gammaPosEnergy[i].size(); k++)
					sum += gammaPosEnergy[i][k] * gammaPosEnergy[j][k];

				aEnergy(i, j) = sum;
			}
		}
	}

	std::vector<float> qForce(pForce.size());
	gmm::lu_solve(aForce, qForce, pForce);
	
	std::vector<float> qEnergy(pEnergy.size());
	gmm::lu_solve(aEnergy, qEnergy, pEnergy);

	std::vector<float> delta(MAX(sPosForce.size(), sPosEnergy.size()));

	for(int i = 0; i < sPosForce.size() || i < sPosEnergy.size(); i++){
		if(i < sPosForce.size())
			delta[i] = qForce[i] * sPosForce[i];
		if(i < sPosEnergy.size())
			delta[i] += qEnergy[i] * sPosEnergy[i];
	}
	
	logfile << "    Force\tq\ts:" << std::endl;
	for(int i = 0; i < sPosForce.size(); i++){
		logfile << "\t" << qForce[i] << "\t" <<  sPosForce[i] << std::endl;
	}
	
	logfile << "    Energy\tq\ts:" << std::endl;
	for(int i = 0; i < sPosEnergy.size(); i++){
		logfile << "\t" << qEnergy[i] << "\t" <<  sPosEnergy[i] << std::endl;
	}

	float u = 0.0;
	for(int i = 0; i < centerForce.size() || i < centerEnergy.size(); i++){
		if(i < centerForce.size())
			u += centerForce[i]*centerForce[i];
		if(i < centerEnergy.size())
			u += centerEnergy[i]*centerEnergy[i];
	}
	logfile << "    U: " << u << " Variables:";
	for(int i = 0; i < variableValues.size(); i++){
		if(i != 0) logfile << ",";
		logfile << " " << variableValues[i];
	}
	logfile << std::endl;
	u = linmin(variableValues, delta);

	logfile << "    U: " << u << " Variables:";
	for(int i = 0; i < finalVariables.size(); i++){
		if(i != 0) logfile << ",";
		logfile << " " << finalVariables[i];
	}
	logfile << std::endl;
}
