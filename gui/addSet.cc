#include "addSet.h"

#include <iostream>

AddSet::AddSet(BaseObjectType *cobject, const Glib::RefPtr<Gnome::Glade::Xml>& refGlade) : Gtk::Window(cobject), refGlade(refGlade) {
	refGlade->get_widget("inputFileView", inputFileView);

	refGlade->get_widget("addInputButton", addInputButton);
	addInputButton->signal_clicked().connect(sigc::mem_fun(*this, &AddSet::on_addInputButton_clicked));

	refGlade->get_widget("removeInputButton", removeInputButton);
	removeInputButton->signal_clicked().connect(sigc::mem_fun(*this, &AddSet::on_removeInputButton_clicked));
	
	refGlade->get_widget("okScanButton", okScanButton);
	okScanButton->signal_clicked().connect(sigc::mem_fun(*this, &AddSet::on_okScanButton_clicked));

	refGlade->get_widget("cancelScanButton", cancelScanButton);
	cancelScanButton->signal_clicked().connect(sigc::mem_fun(*this, &AddSet::on_cancelScanButton_clicked));

	refGlade->get_widget("scanReaderCombo", scanReaderCombo);
	scanModel = Gtk::ListStore::create(scanColumns);
	scanReaderCombo->set_model(scanModel);

	//Populate the drop down box
	Gtk::TreeModel::Row row;
#define ScanReaderClass(Class) \
	row = *(scanModel->append()); \
	row[scanColumns.colName] = Class::name;
#include "../classes.h"
#undef ScanReaderClass
	
	scanReaderCombo->pack_start(scanColumns.colName);
	scanReaderCombo->set_active(0);
	
	inputFileView->append_column("File Name", scanColumns.colName);
	inputModel = Gtk::ListStore::create(scanColumns);
	inputFileView->set_model(inputModel);
	inputModel->clear();
}

AddSet::~AddSet(){
}

GeometrySet & AddSet::getSet(){
	return set;
}

void AddSet::on_addInputButton_clicked(){
	Gtk::FileChooserDialog dialog("Select a file", Gtk::FILE_CHOOSER_ACTION_OPEN);
	dialog.set_transient_for(*this);
	dialog.add_button(Gtk::Stock::CANCEL, Gtk::RESPONSE_CANCEL);
	dialog.add_button(Gtk::Stock::OPEN, Gtk::RESPONSE_OK);

	int result = dialog.run();
	if(result == Gtk::RESPONSE_OK){
		Gtk::TreeModel::Row row = *(inputModel->append());
		row[scanColumns.colName] = dialog.get_filename();
	}
}

void AddSet::on_removeInputButton_clicked(){
	inputModel->erase(inputFileView->get_selection()->get_selected());
}

void AddSet::on_okScanButton_clicked(){
	int i = -1;

	ScanReader * reader;

	if(scanReaderCombo->get_active_row_number() == -1);
#define ScanReaderClass(Class) \
	else if(scanReaderCombo->get_active_row_number() == ++i) reader = new Class();
#include "../classes.h"
#undef ScanReaderClass

	for(Gtk::TreeModel::Children::iterator it = inputModel->children().begin(); it != inputModel->children().end(); it++){
		Glib::ustring file = (*it)[scanColumns.colName];
		reader->inputFiles.push_back(file.raw());
	}

	set = reader->getGeometries();
	hide();
}

void AddSet::on_cancelScanButton_clicked(){
	hide();
}
