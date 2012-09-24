/*
  base/record.h
  Simply records data

  RELACS - Relaxed ELectrophysiological data Acquisition, Control, and Stimulation
  Copyright (C) 2002-2012 Jan Benda <benda@bio.lmu.de>

  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 3 of the License, or
  (at your option) any later version.
  
  RELACS is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
  
  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef _RELACS_BASE_RECORD_H_
#define _RELACS_BASE_RECORD_H_ 1

#include <relacs/repro.h>
using namespace relacs;

namespace base {


/*!
\class Record
\brief [RePro] Simply records data
\author Jan Benda

The Record-RePro simply records data without writing out any stimulus
and terminates after \c repeats times \c duration ms.

\par Options
- \c duration=1000ms: Duration (\c number)
- \c repeats=1: Repeats (\c integer)

\par Files
No output files.

\par Plots
No plot.

\par Requirements
No requirements.

\version 1.0 (Aug 13, 2012)
*/


class Record : public RePro
{
  Q_OBJECT

public:

  Record( void );
  virtual ~Record( void );
  virtual int main( void );

};


}; /* namespace base */

#endif /* ! _RELACS_BASE_RECORD_H_ */