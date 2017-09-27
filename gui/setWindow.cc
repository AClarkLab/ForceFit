#include "setWindow.h"

#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm>

#include "addSet.h"
#include "mdQuestionWindow.h"
#include "minimizerQuestionWindow.h"
#include "geometryWindow.h"

#include "../pugixml/pugixml.hpp"

#define MolecularDynamicsInclude
#define GradientCreatorInclude
#define MinimizerInclude
#include "../classes.h"
#undef MolecularDynamicsInclude
#undef GradientCreatorInclude
#undef MinimizerInclude

SetWindow::SetWindow(BaseObjectType *cobject, const Glib::RefPtr<Gnome::Glade::Xml>& refGlade) : Gtk::Window(cobject), refGlade(refGlade) {
	refGlade->get_widget("setView", setView);
	setModel = Gtk::TreeStore::create(setColumns);
	setView->set_model(setModel);
	setView->signal_row_activated().connect(sigc::mem_fun(*this, &SetWindow::on_setView_row_activated));
	refGlade->get_widget("mdCombo", mdCombo);
	refGlade->get_widget("mdButton", mdButton);
	mdButton->signal_clicked().connect(sigc::mem_fun(*this, &SetWindow::on_mdButton_clicked));
	refGlade->get_widget("gradientCombo", gradientCombo);
	refGlade->get_widget("gradientButton", gradientButton);

	refGlade->get_widget("saveButton", saveButton);
	saveButton->signal_clicked().connect(sigc::mem_fun(*this, &SetWindow::on_saveButton_clicked));
	refGlade->get_widget("openButton", openButton);
	openButton->signal_clicked().connect(sigc::mem_fun(*this, &SetWindow::on_openButton_clicked));

	refGlade->get_widget("variablesView", variablesView);

	variablesModel = Gtk::ListStore::create(variablesColumns);
	variablesView->set_model(variablesModel);
	variablesModel->append();
	variablesView->append_column("Name", variablesColumns.colVar);
	variablesView->append_column_editable("Value", variablesColumns.colValue);
	variablesModel->signal_row_changed().connect(sigc::mem_fun(*this, &SetWindow::on_variablesModel_changed));
	
	refGlade->get_widget("addSetButton", addSetButton);
	addSetButton->signal_clicked().connect(sigc::mem_fun(*this, &SetWindow::on_setAddButton_clicked));
	
	refGlade->get_widget("removeSetButton", removeSetButton);
	removeSetButton->signal_clicked().connect(sigc::mem_fun(*this, &SetWindow::on_setRemoveButton_clicked));
		
	setView->append_column("Name", setColumns.colName);

	Gtk::TreeModel::Row row;

	mdModel = Gtk::ListStore::create(setColumns);
	mdCombo->set_model(mdModel);
#define MolecularDynamicsClass(Class) \
	row = *(mdModel->append()); \
	row[setColumns.colName] = Class::name;
#include "../classes.h"
#undef MolecularDynamicsClass
	mdCombo->pack_start(setColumns.colName);
	mdCombo->set_active(0);

	Glib::RefPtr<Gtk::ListStore> gradientModel = Gtk::ListStore::create(setColumns);
	gradientCombo->set_model(gradientModel);
#define GradientCreatorClass(Class) \
	row = *(gradientModel->append()); \
	row[setColumns.colName] = Class::name;
#include "../classes.h"
#undef GradientCreatorClass
	gradientCombo->pack_start(setColumns.colName);
	gradientCombo->set_active(0);
}

SetWindow::~SetWindow(){
}

void SetWindow::addGeometrySet(GeometrySet * set){
	Gtk::TreeModel::Row row = *(setModel->append());
	Gtk::TreeModel::Row childrow;
	row[setColumns.colName] = set->name;
	row[setColumns.set] = set;

	int n = 1;
	for(std::vector<Geometry>::iterator geomit = set->geometries.begin(); geomit != set->geometries.end(); geomit++){
		childrow = *(setModel->append(row.children()));
		std::ostringstream oss;
		oss << "Geometry " << n++;
		childrow[setColumns.colName] = oss.str();
		childrow[setColumns.geomPtr] = &(*geomit);
	}
}

void SetWindow::on_setAddButton_clicked(){
	AddSet *addSetWindow;
	refGlade->get_widget_derived("addSetWindow", addSetWindow);
	Gtk::Main::run(*addSetWindow);

	if(addSetWindow->getSet().geometries.size() > 0){
		GeometrySet * set = new GeometrySet(addSetWindow->getSet());
		addGeometrySet(set);

	}
}

void SetWindow::on_setRemoveButton_clicked(){
	if(!setView->get_selection()->get_selected()->parent()) // Erasing a geometry set
		setModel->erase(setView->get_selection()->get_selected());
	else { // Erasing a geometry
		Gtk::TreeModel::iterator sel = setView->get_selection()->get_selected();
		GeometrySet * parent = (*sel->parent())[setColumns.set];
		//Geometry * geom = (*sel)[setColumns.geomPtr];
		//parent.geometries.erase(std::find(parent.geometries.begin(), parent.geometries.end(),	*geom));
	}
}

void SetWindow::on_setView_row_activated(const Gtk::TreeModel::Path& path, Gtk::TreeViewColumn* column){
	Gtk::TreeModel::iterator iter = setModel->get_iter(path);

	if(!iter)
		return;

	Gtk::TreeModel::Row row = *iter;

	if(!row->parent()){
	} else {
		Geometry * geom = row[setColumns.geomPtr];

		if(!geom){
			std::cout << "Invalid geometry" << std::endl;
			return;
		}

		/*
		GeometryWindow *geomWindow;
		refGlade->get_widget_derived("geometryEditWindow", geomWindow);

		geomWindow->setup(geom);

		Gtk::Main::run(*geomWindow);
		*/
	}
}

bool changingVariablesModel = false;

void SetWindow::on_variablesModel_changed(const Gtk::TreeModel::Path& path, const Gtk::TreeModel::iterator& iter){
	if(changingVariablesModel)
		return;

	changingVariablesModel = true;
	std::vector<Gtk::TreeRow> toErase;
	std::string var = "A";
	for(Gtk::TreeModel::Children::iterator chIt = variablesModel->children().begin(); chIt != variablesModel->children().end(); chIt++){
		chIt->set_value(variablesColumns.colVar, Glib::ustring(var));

		if((*chIt)[variablesColumns.colValue] == 0.0)
			toErase.push_back(*chIt);
		else
			var[0]++;
	}
	for(std::vector<Gtk::TreeRow>::iterator erIt = toErase.begin(); erIt != toErase.end(); erIt++)
		variablesModel->erase(*erIt);

	variablesModel->append();

	changingVariablesModel = false;
}

void SetWindow::on_mdButton_clicked(){
	std::vector<GeometrySet> geoSets;
	std::vector<GeometrySet *> setPtrs;
	std::vector<float> variables;
	
	for(Gtk::TreeModel::Children::iterator chIt = setModel->children().begin(); chIt != setModel->children().end(); chIt++){
		geoSets.push_back(*(*chIt)[setColumns.set]);
		setPtrs.push_back((*chIt)[setColumns.set]);
	}

	for(Gtk::TreeModel::Children::iterator chIt = variablesModel->children().begin(); chIt != variablesModel->children().end(); chIt++){
		if((*chIt)[variablesColumns.colValue] != 0.0)
			variables.push_back((*chIt)[variablesColumns.colValue]);
	}

	MolecularDynamics * md;
	Gtk::TreeRow row = *(mdCombo->get_active());
	if(0);
#define MolecularDynamicsClass(Class) \
	else if(row[setColumns.colName] == Class::name) md = new Class();
#include "../classes.h"
#undef MolecularDynamicsClass

	std::vector<struct mdQuestion> questions = md->getQuestions();
	MDQuestionWindow mdQuestionWindow(&questions);
	Gtk::Main::run(mdQuestionWindow);

	if(mdQuestionWindow.ok()){
		try {
			md->prepareMD(geoSets, variables, questions);
			Minimizer * min;
			if(0);
#define MinimizerClass(Class) \
			else if(mdQuestionWindow.getMin() == Class::name) min = new Class();
#include "../classes.h"
#undef MinimizerClass

			std::vector<struct minQuestion> minQuestions = min->getQuestions();
			MinimizerQuestionWindow minQuestionWindow(&minQuestions);
			Gtk::Main::run(minQuestionWindow);
			if(minQuestionWindow.ok()){
				min->prepareMin(geoSets, minQuestions);
				min->runMin(md);
			}

			delete min;
		} catch(const MolecularDynamicsException & e){
			Gtk::MessageDialog dialog(*this, "Molecular Dynamics Error!", false, Gtk::MESSAGE_ERROR);
			dialog.set_secondary_text(e.what());
			dialog.run();
		}
	}

	delete md;
}

void SetWindow::on_saveButton_clicked(){
	Gtk::FileChooserDialog dialog("Select a file", Gtk::FILE_CHOOSER_ACTION_SAVE);
	dialog.set_transient_for(*this);
	dialog.add_button(Gtk::Stock::CANCEL, Gtk::RESPONSE_CANCEL);
	dialog.add_button(Gtk::Stock::SAVE, Gtk::RESPONSE_OK);

	int result = dialog.run();
	if(result == Gtk::RESPONSE_OK){
		std::ofstream saveFile(dialog.get_filename().c_str());

		saveFile << "<forcefit>" << std::endl;

		for(Gtk::TreeModel::Children::iterator chIt = setModel->children().begin(); chIt != setModel->children().end(); chIt++){
			(*(*chIt)[setColumns.set]).write(saveFile, 1);
		}

		saveFile << "\t<variables>" << std::endl;
		for(Gtk::TreeModel::Children::iterator chIt = variablesModel->children().begin(); chIt != variablesModel->children().end(); chIt++){
			if((*chIt)[variablesColumns.colValue] != 0.0){
				saveFile << "\t\t<variable>" << (*chIt)[variablesColumns.colValue] << "</variable>" << std::endl;
			}
		}
		saveFile << "\t</variables>" << std::endl;
		
		saveFile << "</forcefit>" << std::endl;

		saveFile.close();
	}
}

void SetWindow::on_openButton_clicked(){
	Gtk::FileChooserDialog dialog("Select a file", Gtk::FILE_CHOOSER_ACTION_OPEN);
	dialog.set_transient_for(*this);
	dialog.add_button(Gtk::Stock::CANCEL, Gtk::RESPONSE_CANCEL);
	dialog.add_button(Gtk::Stock::OPEN, Gtk::RESPONSE_OK);

	int result = dialog.run();
	if(result == Gtk::RESPONSE_OK){
		pugi::xml_document xmlDoc;
		pugi::xml_parse_result result = xmlDoc.load_file(dialog.get_filename().c_str());

		for(pugi::xml_node geometrySetNode = xmlDoc.child("forcefit").child("geometrySet"); geometrySetNode; geometrySetNode = geometrySetNode.next_sibling("geometrySet")){
			GeometrySet * set = new GeometrySet();

			set->load(geometrySetNode);

			addGeometrySet(set);
		}

		{
			variablesModel = Gtk::ListStore::create(variablesColumns);
			std::string var = "A";
			for(pugi::xml_node variable = xmlDoc.child("forcefit").child("variables").child("variable"); variable; variable = variable.next_sibling("variable")){
				Gtk::TreeModel::Row row = *(variablesModel->append());
				row[variablesColumns.colVar] = var;
				row[variablesColumns.colValue] = atof(variable.child_value());
				var[0]++;
			}
			variablesModel->append();
			variablesModel->signal_row_changed().connect(sigc::mem_fun(*this, &SetWindow::on_variablesModel_changed));
			variablesView->set_model(variablesModel);
			variablesView->remove_all_columns();
			variablesView->append_column("Name", variablesColumns.colVar);
			variablesView->append_column_editable("Value", variablesColumns.colValue);
			
		}
	}
}
