/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       David Salinas
 *
 *    Copyright (C) 2014 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

#include "MainWindow.h"
#include <QApplication>
#include <CGAL/Qt/resources.h>

int main(int argc, char** argv) {
  QApplication application(argc, argv);
  application.setOrganizationDomain("inria.fr");
  application.setOrganizationName("Inria");
  application.setApplicationName("GudhUI");

  MainWindow mw;
  application.setQuitOnLastWindowClosed(false);
  mw.show();
  return application.exec();
}
