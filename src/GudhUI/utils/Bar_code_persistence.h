/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       David Salinas
 *
 *    Copyright (C) 2014 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

#include <math.h>  // isfinite

#include <QtGui/QApplication>

#include <QGraphicsView>
#include <QGraphicsScene>
#include <QPointF>
#include <QVector>
#include <QGraphicsTextItem>

#include <iostream>
#include <vector>
#include <limits>  // NaN, infinity
#include <utility>  // for pair
#include <string>

#ifndef UTILS_BAR_CODE_PERSISTENCE_H_
#define UTILS_BAR_CODE_PERSISTENCE_H_

class Bar_code_persistence {
 private:
  typedef std::vector<std::pair<double, double>> Persistence;
  Persistence persistence_vector;
  double min_birth;
  double max_death;

 public:
  Bar_code_persistence()
      : min_birth(std::numeric_limits<double>::quiet_NaN()),
      max_death(std::numeric_limits<double>::quiet_NaN()) { }

  void insert(double birth, double death) {
    persistence_vector.push_back(std::make_pair(birth, death));
    if (std::isfinite(birth)) {
      if ((birth < min_birth) || (std::isnan(min_birth)))
        min_birth = birth;
      if ((birth > max_death) || (std::isnan(max_death)))
        max_death = birth;
    }
    if (std::isfinite(death))
      if ((death > max_death) || (std::isnan(max_death)))
        max_death = death;
  }

  void show(const std::string& window_title) {
    // Create a view, put a scene in it
    QGraphicsView * view = new QGraphicsView();
    QGraphicsScene * scene = new QGraphicsScene();
    view->setScene(scene);
    double ratio = 600.0 / (max_death - min_birth);
    // std::clog << "min_birth=" << min_birth << " - max_death=" << max_death << " - ratio=" << ratio << std::endl;

    double height = 0.0, birth = 0.0, death = 0.0;
    int pers_num = 1;
    for (auto& persistence : persistence_vector) {
      height = 5.0 * pers_num;
      // std::clog << "[" << pers_num << "] birth=" << persistence.first << " - death=" << persistence.second << std::endl;
      if (std::isfinite(persistence.first))
        birth = ((persistence.first - min_birth) * ratio) + 50.0;
      else
        birth = 0.0;

      if (std::isfinite(persistence.second))
        death = ((persistence.second - min_birth) * ratio) + 50.0;
      else
        death = 700.0;

      scene->addLine(birth, height, death, height, QPen(Qt::blue, 2));
      pers_num++;
    }
    height += 10.0;
    // scale line
    scene->addLine(0, height, 700.0, height, QPen(Qt::black, 1));
    int modulo = 0;
    for (double scale = 50.0; scale < 700.0; scale += 50.0) {
      modulo++;
      // scale small dash
      scene->addLine(scale, height - 3.0, scale, height + 3.0, QPen(Qt::black, 1));
      // scale text
      QString scale_value = QString::number(((scale - 50.0) / ratio) + min_birth);
      QGraphicsTextItem* dimText = scene->addText(scale_value, QFont("Helvetica", 8));
      dimText->setPos(scale - (3.0 * scale_value.size()), height + 9.0 * (modulo % 2));
    }
    view->setWindowTitle(window_title.c_str());
    // Show the view
    view->show();
  }
};

#endif  // UTILS_BAR_CODE_PERSISTENCE_H_
