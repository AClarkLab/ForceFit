#include "geometryWindow.h"

GeometryWindow::GeometryWindow(BaseObjectType *cobject, const Glib::RefPtr<Gnome::Glade::Xml>& refGlade) : Gtk::Window(cobject), refGlade(refGlade) {
}

void GeometryWindow::setup(Geometry * geom){
	refGlade->get_widget("modeEntry", modeEntry);
	refGlade->get_widget("stepEntry", stepEntry);
	refGlade->get_widget("energyEntry", energyEntry);
	refGlade->get_widget("atomView", atomView);
	refGlade->get_widget("geometryOkButton", okButton);
	refGlade->get_widget("geometryCancelButton", cancelButton);

	atomModel = Gtk::ListStore::create(atomColumns);

	atomView->append_column("#", atomColumns.colNumber);
	atomView->append_column("Atomic #", atomColumns.colAtomicNumber);
	atomView->append_column("Symbol", atomColumns.colSymbol);
	atomView->append_column("Mass", atomColumns.colMass);
	atomView->append_column("Charge", atomColumns.colCharge);
	atomView->append_column("X", atomColumns.colX);
	atomView->append_column("Y", atomColumns.colY);
	atomView->append_column("Z", atomColumns.colZ);
	atomView->append_column("Force X", atomColumns.colForceX);
	atomView->append_column("Force Y", atomColumns.colForceY);
	atomView->append_column("Force Z", atomColumns.colForceZ);
	
	atomView->set_model(atomModel);

	for(std::vector<Atom>::iterator it = geom->atoms.begin(); it != geom->atoms.end(); it++){
		Gtk::TreeModel::Row row;
		row = *(atomModel->append());
		
		row[atomColumns.colNumber] = it - geom->atoms.begin() + 1;
		row[atomColumns.colAtomicNumber] = it->number;
		//row[atomColumns.colSymbol] = it->symbol;
		row[atomColumns.colMass] = it->mass;
		row[atomColumns.colCharge] = it->charge;
		row[atomColumns.colX] = it->x;
		row[atomColumns.colY] = it->y;
		row[atomColumns.colZ] = it->z;
		row[atomColumns.colForceX] = it->forcex;
		row[atomColumns.colForceY] = it->forcey;
		row[atomColumns.colForceZ] = it->forcez;
	}

	cancelButton->signal_clicked().connect(sigc::mem_fun(*this, &GeometryWindow::on_cancelButton_clicked));
	okButton->signal_clicked().connect(sigc::mem_fun(*this, &GeometryWindow::on_okButton_clicked));
	show_all_children();
}

void GeometryWindow::on_cancelButton_clicked(){
	hide();
}

void GeometryWindow::on_okButton_clicked(){
	hide();
}
