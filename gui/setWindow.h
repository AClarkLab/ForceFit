#ifndef _SETWINDOW_H
#define _SETWINDOW_H

#include <gtkmm.h>
#include <libglademm.h>

#include <vector>

#include "../geometrySet.h"

class SetWindow : public Gtk::Window {
public:
	SetWindow(BaseObjectType *cobject, const Glib::RefPtr<Gnome::Glade::Xml>& refGlade);
	~SetWindow();
protected:
	virtual void on_mdButton_clicked();
	virtual void on_variablesModel_changed(const Gtk::TreeModel::Path& path, const Gtk::TreeModel::iterator& iter);
	virtual void on_setAddButton_clicked();
	virtual void on_setRemoveButton_clicked();
	virtual void on_setView_row_activated(const Gtk::TreeModel::Path& path, Gtk::TreeViewColumn* column);
	virtual void on_openButton_clicked();
	virtual void on_saveButton_clicked();

	Glib::RefPtr<Gnome::Glade::Xml> refGlade;
	
	Gtk::Button *addSetButton;
	Gtk::Button *removeSetButton;
	Gtk::TreeView *setView;
	Gtk::ComboBox *mdCombo;
	Gtk::Button *mdButton;
	Gtk::ComboBox *gradientCombo;
	Gtk::Button *gradientButton;

	Gtk::Button *saveButton;
	Gtk::Button *openButton;

	Gtk::TreeView *variablesView;

	Glib::RefPtr<Gtk::ListStore> mdModel;

	class ModelColumns : public Gtk::TreeModel::ColumnRecord
	{
	public:
		ModelColumns() { add(colName); add(set); add(geomPtr); }
		Gtk::TreeModelColumn<Glib::ustring> colName;

		Gtk::TreeModelColumn<GeometrySet *> set;
		Gtk::TreeModelColumn<Geometry *> geomPtr;
	};
	
	ModelColumns setColumns;
	Glib::RefPtr<Gtk::TreeStore> setModel;

	class VariableColumns : public Gtk::TreeModel::ColumnRecord
	{
	public:
		VariableColumns() { add(colVar); add(colValue); }
		Gtk::TreeModelColumn<Glib::ustring> colVar;
		Gtk::TreeModelColumn<float> colValue;
	};

	VariableColumns variablesColumns;
	Glib::RefPtr<Gtk::ListStore> variablesModel;

private:
	void addGeometrySet(GeometrySet * set);
};

#endif
