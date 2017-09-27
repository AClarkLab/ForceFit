#include "mdQuestionWindow.h"

#define MinimizerInclude
#include "../classes.h"
#undef MinimizerInclude

MDQuestionWindow::MDQuestionWindow(std::vector<struct mdQuestion> * qs) : table(3,qs->size()+2,false), minimizerLabel("Minimizer:"), okButton(Gtk::Stock::OK), cancelButton(Gtk::Stock::CANCEL), okClicked(false), questions(qs) {
	set_title("Molecular Dynamics Run");

	add(table);

	minimizerModel = Gtk::ListStore::create(minimizerColumns);
	minimizerCombo.set_model(minimizerModel);

	Gtk::TreeModel::Row row;
#define MinimizerClass(Class) \
	row = *(minimizerModel->append()); \
	row[minimizerColumns.colName] = Class::name;
#include "../classes.h"
#undef ScanReaderClass

	minimizerCombo.pack_start(minimizerColumns.colName);
	minimizerCombo.set_active(0);

	table.attach(minimizerLabel, 0, 1, 0, 1);
	table.attach(minimizerCombo, 1, 3, 0, 1);

	int i = 0;
	for(std::vector<struct mdQuestion>::iterator it = questions->begin(); it != questions->end(); it++){
		Gtk::Label * label = new Gtk::Label(it->phrase);
		questionLabels.push_back(label);
		table.attach(*label, 0, 1, i+1, i+2);
		Gtk::Entry * entry;

		if(it->type == MDBOOL) {
			Gtk::CheckButton * check = new Gtk::CheckButton();
			check->set_active(it->boolAnswer);
			questionFields.push_back((Gtk::Container *)check);
		} else {
			entry = new Gtk::Entry();
			entry->set_text(it->stringAnswer);
			questionFields.push_back((Gtk::Container *)entry);
		}
		if(it->type != MDFILENAME) {
			table.attach(*(questionFields[i]),1,3,i+1,i+2);
		} else {
			table.attach(*(questionFields[i]),1,2,i+1,i+2);
			Gtk::Button * browse = new Gtk::Button("Browse...");
			browse->signal_clicked().connect(sigc::bind(sigc::mem_fun(*this, &MDQuestionWindow::on_browseButton_clicked), entry));
			browseButtons.push_back(browse);
			table.attach(*browse,2,3,i+1,i+2);
		}
		i++;
	}

	buttonBox.add(cancelButton);
	buttonBox.add(okButton);
	table.attach(buttonBox, 0, 3, i+1, i+2);

	cancelButton.signal_clicked().connect(sigc::mem_fun(*this, &MDQuestionWindow::on_cancelButton_clicked));
	okButton.signal_clicked().connect(sigc::mem_fun(*this, &MDQuestionWindow::on_okButton_clicked));

	show_all_children();
}

MDQuestionWindow::~MDQuestionWindow(){
	for(std::vector<Gtk::Label *>::iterator it = questionLabels.begin(); it != questionLabels.end(); it++)
		delete (*it);

	for(std::vector<Gtk::Container *>::iterator it = questionFields.begin(); it != questionFields.end(); it++)
		delete (*it);

	for(std::vector<Gtk::Button *>::iterator it = browseButtons.begin(); it != browseButtons.end(); it++)
		delete (*it);
}

const bool & MDQuestionWindow::ok(){
	return okClicked;
}

std::string MDQuestionWindow::getMin(){
	Gtk::TreeRow row = *(minimizerCombo.get_active());
	Glib::ustring str = row[minimizerColumns.colName];
	return str.raw();
}

void MDQuestionWindow::on_cancelButton_clicked(){
	okClicked = false;
	hide();
}

void MDQuestionWindow::on_okButton_clicked(){
	for(int i = 0; i < questions->size() && i < questionFields.size(); i++){
		if((*questions)[i].type == MDBOOL)
			(*questions)[i].boolAnswer = ((Gtk::CheckButton *)questionFields[i])->get_active();
		else
			(*questions)[i].stringAnswer = ((Gtk::Entry *)questionFields[i])->get_text();
	}

	okClicked = true;
	hide();
}

void MDQuestionWindow::on_browseButton_clicked(Gtk::Entry *field){
	Gtk::FileChooserDialog dialog("Select a file", Gtk::FILE_CHOOSER_ACTION_OPEN);
	dialog.set_transient_for(*this);
	dialog.add_button(Gtk::Stock::CANCEL, Gtk::RESPONSE_CANCEL);
	dialog.add_button(Gtk::Stock::OPEN, Gtk::RESPONSE_OK);

	int result = dialog.run();
	if(result == Gtk::RESPONSE_OK){
		field->set_text(dialog.get_filename());
	}
}
