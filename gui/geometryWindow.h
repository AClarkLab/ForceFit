#ifndef _GEOMETRYWINDOW_H
#define _GEOMETRYWINDOW_H

#include <gtkmm.h>
#include <libglademm.h>

#include <vector>

#include "../geometry.h"

class GeometryWindow : public Gtk::Window {
	public:
		GeometryWindow(BaseObjectType *cobject, const Glib::RefPtr<Gnome::Glade::Xml>& refGlade);

		void setup(Geometry * geom);
	protected:
		class AtomColumns : public Gtk::TreeModel::ColumnRecord {
		public:
			AtomColumns() { add(colNumber); add(colAtomicNumber); add(colSymbol); add(colMass); add(colCharge); add(colX); add(colY); add(colZ); add(colForceX); add(colForceY); add(colForceZ); }
			Gtk::TreeModelColumn<int> colNumber;
			Gtk::TreeModelColumn<int> colAtomicNumber;
			Gtk::TreeModelColumn<Glib::ustring> colSymbol;
			Gtk::TreeModelColumn<float> colMass;
			Gtk::TreeModelColumn<float> colCharge;
			Gtk::TreeModelColumn<float> colX, colY, colZ;
			Gtk::TreeModelColumn<float> colForceX, colForceY, colForceZ;
		};

		AtomColumns atomColumns;
		Glib::RefPtr<Gtk::ListStore> atomModel;

		Glib::RefPtr<Gnome::Glade::Xml> refGlade;

		Gtk::Entry *modeEntry, *stepEntry, *energyEntry;
		Gtk::TreeView *atomView;

		Gtk::Button *okButton, *cancelButton;

		void on_cancelButton_clicked();
		void on_okButton_clicked();
};

#endif
