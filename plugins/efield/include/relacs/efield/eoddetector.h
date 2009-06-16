/*
  efield/eoddetector.cc
  A detector for EOD cycles of weakly electric fish

  RELACS - Relaxed ELectrophysiological data Acquisition, Control, and Stimulation
  Copyright (C) 2002-2009 Jan Benda <j.benda@biologie.hu-berlin.de>

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

#ifndef _RELACS_EFIELD_EODDETECTOR_H_
#define _RELACS_EFIELD_EODDETECTOR_H_ 1

#include <relacs/optwidget.h>
#include <relacs/detector.h>
#include <relacs/filter.h>
using namespace relacs;

namespace efield {


/*!
\class EODDetector
\brief [Detector] A detector for EOD cycles of weakly electric fish
\author Jan Benda
\version 1.3 (Jun 16, 2009)
*/


class EODDetector : public Filter
{
  Q_OBJECT

public:

  EODDetector( const string &ident="", int mode=0 );
  ~EODDetector( void );

  virtual int init( const InData &data, EventData &outevents,
		     const EventList &other, const EventData &stimuli );
  virtual void notify( void );
  virtual int adjust( const InData &data );
    /*! Detect EODs in a single trace of the analog data \a data. */
  virtual int detect( const InData &data, EventData &outevents,
		      const EventList &other, const EventData &stimuli );

  int checkEvent( const InData::const_iterator &first, 
		  const InData::const_iterator &last,
		  InData::const_iterator &event, 
		  InDataTimeIterator &eventtime, 
		  InData::const_iterator &index,
		  InDataTimeIterator &indextime, 
		  InData::const_iterator &prevevent, 
		  InDataTimeIterator &prevtime, 
		  EventData &outevents,
		  double &threshold,
		  double &minthresh, double &maxthresh,
		  double &time, double &size, double &width );


protected:

  Detector< InData::const_iterator, InDataTimeIterator > D;

  double Threshold;
  double MinThresh;
  double MaxThresh;
  double ThreshRatio;

    /*! Maximum period of the EOD to detect in seconds. */
  double MaxEODPeriod;

  OptWidget EDW;

};


}; /* namespace efield */

#endif /* ! _RELACS_EFIELD_EODDETECTOR_H_ */