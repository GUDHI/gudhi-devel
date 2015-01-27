#include "MainWindow.h"
#include <QApplication>
#include <CGAL/Qt/resources.h>


int main(int argc, char** argv)
{
 QApplication application(argc,argv);
 application.setOrganizationDomain("inria.fr");
 application.setOrganizationName("INRIA");
 application.setApplicationName("GudhUI");
  
  MainWindow mw;
  application.setQuitOnLastWindowClosed(false);
  mw.show();
  return application.exec();
}
