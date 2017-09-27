#include "minimizerQuestionWindow.h"

MinimizerQuestionWindow::MinimizerQuestionWindow(std::vector<struct minQuestion> * qs) : table(3,qs->size()+1,false), okButton(Gtk::Stock::OK), cancelButton(Gtk::Stock::CANCEL), okClicked(false), questions(qs) {
	set_title("Minimizer Run");

	add(table);

	int i = 0;
	for(std::vector<struct minQuestion>::iterator it = questions->begin(); it != questions->end(); it++){
		Gtk::Label * label = new Gtk::Label(it->phrase);
		questionLabels.push_back(label);
		table.attach(*label, 0, 1, i+1, i+2);
		Gtk::Entry * entry;

		if(it->type == MINBOOL) {
			Gtk::CheckButton * check = new Gtk::CheckButton();
			check->set_active(it->boolAnswer);
			questionFields.push_back((Gtk::Container *)check);
		} else {
			entry = new Gtk::Entry();
			entry->set_text(it->stringAnswer);
			questionFields.push_back((Gtk::Container *)entry);
		}
		if(it->type != MINFILENAME) {
			table.attach(*(questionFields[i]),1,3,i+1,i+2);
		} else {
			table.attach(*(questionFields[i]),1,2,i+1,i+2);
			Gtk::Button * browse = new Gtk::Button("Browse...");
			browse->signal_clicked().connect(sigc::bind(sigc::mem_fun(*this, &MinimizerQuestionWindow::on_browseButton_clicked), entry));
			browseButtons.push_back(browse);
			table.attach(*browse,2,3,i+1,i+2);
		}
		i++;
	}

	buttonBox.add(cancelButton);
	buttonBox.add(okButton);
	table.attach(buttonBox, 0, 3, i+1, i+2);

	cancelButton.signal_clicked().connect(sigc::mem_fun(*this, &MinimizerQuestionWindow::on_cancelButton_clicked));
	okButton.signal_clicked().connect(sigc::mem_fun(*this, &MinimizerQuestionWindow::on_okButton_clicked));

	show_all_children();
}

MinimizerQuestionWindow::~MinimizerQuestionWindow(){
	for(std::vector<Gtk::Label *>::iterator it = questionLabels.begin(); it != questionLabels.end(); it++)
		delete (*it);

	for(std::vector<Gtk::Container *>::iterator it = questionFields.begin(); it != questionFields.end(); it++)
		delete (*it);

	for(std::vector<Gtk::Button *>::iterator it = browseButtons.begin(); it != browseButtons.end(); it++)
		delete (*it);
}

const bool & MinimizerQuestionWindow::ok(){
	return okClicked;
}

void MinimizerQuestionWindow::on_cancelButton_clicked(){
	okClicked = false;
	hide();
}

void MinimizerQuestionWindow::on_okButton_clicked(){
	for(int i = 0; i < questions->size() && i < questionFields.size(); i++){
		if((*questions)[i].type == MINBOOL)
			(*questions)[i].boolAnswer = ((Gtk::CheckButton *)questionFields[i])->get_active();
		else
			(*questions)[i].stringAnswer = ((Gtk::Entry *)questionFields[i])->get_text();
	}

	okClicked = true;
	hide();
}

void MinimizerQuestionWindow::on_browseButton_clicked(Gtk::Entry *field){
	Gtk::FileChooserDialog dialog("Select a file", Gtk::FILE_CHOOSER_ACTION_OPEN);
	dialog.set_transient_for(*this);
	dialog.add_button(Gtk::Stock::CANCEL, Gtk::RESPONSE_CANCEL);
	dialog.add_button(Gtk::Stock::OPEN, Gtk::RESPONSE_OK);

	int result = dialog.run();
	if(result == Gtk::RESPONSE_OK){
		field->set_text(dialog.get_filename());
	}
}
