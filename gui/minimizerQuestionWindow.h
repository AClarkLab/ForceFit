#ifndef _MINIMIZERQUESTIONWINDOW_H
#define _MINIMIZERQUESTIONWINDOW_H

#include <gtkmm.h>

#include <vector>

#include "../minimizer/minimizer.h"

class MinimizerQuestionWindow : public Gtk::Window {
	public:
		MinimizerQuestionWindow(std::vector<struct minQuestion> * qs);
		virtual ~MinimizerQuestionWindow();

		// Returns true if the OK button was clicked
		const bool & ok();
	protected:
		std::vector<struct minQuestion> * questions;
		Gtk::Table table;

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
