#ifndef _ADDSET_H
#define _ADDSET_H

#include <gtkmm.h>
#include <libglademm.h>

#include <vector>

#define ScanReaderInclude
#include "../classes.h"
#undef ScanReaderInclude

class AddSet : public Gtk::Window {
public:
	AddSet(BaseObjectType *cobject, const Glib::RefPtr<Gnome::Glade::Xml>& refGlade);
	~AddSet();

	GeometrySet & getSet();
protected:
	virtual void on_addInputButton_clicked();
	virtual void on_removeInputButton_clicked();
	virtual void on_okScanButton_clicked();
	virtual void on_cancelScanButton_clicked();

	GeometrySet set;

	Glib::RefPtr<Gnome::Glade::Xml> refGlade;

	Gtk::TreeView *inputFileView;
	Glib::RefPtr<Gtk::ListStore> inputModel;

	Gtk::Button *addInputButton;
	Gtk::Button *removeInputButton;
	Gtk::Button *okScanButton;
	Gtk::Button *cancelScanButton;

	//Combo Box
	class ModelColumns : public Gtk::TreeModel::ColumnRecord
	{
	public:
		ModelColumns() { add(colName); }
		Gtk::TreeModelColumn<Glib::ustring> colName;
	};
	
	ModelColumns scanColumns;
	Glib::RefPtr<Gtk::ListStore> scanModel;
	Gtk::ComboBox *scanReaderCombo;
private:
};

#endif
