#ifndef _MDQUESTIONWINDOW_H
#define _MDQUESTIONWINDOW_H

#include <gtkmm.h>

#include <vector>

#include "../molecularDynamics/molecularDynamics.h"

class MDQuestionWindow : public Gtk::Window {
	public:
		MDQuestionWindow(std::vector<struct mdQuestion> * qs);
		virtual ~MDQuestionWindow();

		// Returns true if the OK button was clicked
		const bool & ok();

		std::string getMin();
	protected:
		std::vector<struct mdQuestion> * questions;
		Gtk::Table table;

		class MinimizerColumns : public Gtk::TreeModel::ColumnRecord {
		public:
			MinimizerColumns() { add(colName); }
			Gtk::TreeModelColumn<Glib::ustring> colName;
		};

		MinimizerColumns minimizerColumns;
		Glib::RefPtr<Gtk::ListStore> minimizerModel;
		Gtk::Label minimizerLabel;
		Gtk::ComboBox minimizerCombo;

		std::vector<Gtk::Label *> questionLabels;
		std::vector<Gtk::Container *> questionFields;
		std::vector<Gtk::Button *> browseButtons;

		Gtk::HButtonBox buttonBox;
		Gtk::Button okButton, cancelButton;

		bool okClicked;

		void on_cancelButton_clicked();
		void on_okButton_clicked();
		void on_browseButton_clicked(Gtk::Entry *field);
};

#endif
