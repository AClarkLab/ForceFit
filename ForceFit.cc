#include <gtkmm.h>
#include <libglademm/xml.h>
#include <cmath>

#include "gui/setWindow.h"

int main(int argc, char **argv){
	//Initialize GTK
	Gtk::Main kit(argc, argv);

	//Load the glade file
	Glib::RefPtr<Gnome::Glade::Xml> refXml = Gnome::Glade::Xml::create("interface.glade");
	
	SetWindow *setWindow;
	refXml->get_widget_derived("setWindow", setWindow);

	kit.run(*setWindow);

	delete setWindow;

	return 0;
}
