/*--------------------------------------------------------------------------*\
 |                                                                          |
 |  Copyright (C) 2023                                                      |
 |                                                                          |
 |         , __                 , __                                        |
 |        /|/  \               /|/  \                                       |
 |         | __/ _   ,_         | __/ _   ,_                                |
 |         |   \|/  /  |  |   | |   \|/  /  |  |   |                        |
 |         |(__/|__/   |_/ \_/|/|(__/|__/   |_/ \_/|/                       |
 |                           /|                   /|                        |
 |                           \|                   \|                        |
 |                                                                          |
 |      Enrico Bertolazzi                                                   |
 |      Dipartimento di Ingegneria Industriale                              |
 |      Universita` degli Studi di Trento                                   |
 |      email: enrico.bertolazzi@unitn.it                                   |
 |                                                                          |
\*--------------------------------------------------------------------------*/

///
/// file: Utils_progress_bar.cc
///

#include "Utils.hh"

namespace Utils {

  string
  progress_bar( double progress, int width ) {
    string res{"["};
    int pos = width * progress;
    for (int i = 0; i < width; ++i) {
      if      (i  < pos) res += '=';
      else if (i == pos) res += '>';
      else               res += '_';
    }
    res += fmt::format( "] {}%", int(100*progress) );
    return res;
  }

}

///
/// eof: Utils_progress_bar.c
///
