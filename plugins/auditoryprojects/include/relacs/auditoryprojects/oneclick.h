/*
  auditoryprojects/oneclick.h
  

  RELACS - Relaxed ELectrophysiological data Acquisition, Control, and Stimulation
  Copyright (C) 2002-2015 Jan Benda <jan.benda@uni-tuebingen.de>

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

#ifndef _RELACS_AUDITORYPROJECTS_ONECLICK_H_
#define _RELACS_AUDITORYPROJECTS_ONECLICK_H_ 1

#include <relacs/repro.h>
#include <relacs/rangeloop.h>
#include <relacs/multiplot.h>
#include <relacs/ephys/traces.h>
#include <relacs/acoustic/traces.h>
using namespace relacs;

namespace auditoryprojects {


/*!
\class OneClick
\brief [RePro] A single short click
\author Alexander Wolf

\par Options
- Intensities
- \c intmin=30dB SPL: Minimum Click intensity (\c number)
- \c intmax=100dB SPL: Maximum Click intensity (\c number)
- \c intstep=5dB SPL: Click intensity step (\c number)
- \c repeat=15: Number of repetitions of the whole f-I curve measurement (\c integer)
- Waveform
- \c duration=2000microsec: Duration of stimulus (\c number)
- \c latency=2ms: Latency after stimulus (\c number)
- \c pause=400ms: Pause (\c number)
- \c side=left: Speaker (\c string)

\version 0.2 (Jan 10, 2008)
*/


class OneClick : public RePro, public ephys::Traces, public acoustic::Traces
{
  Q_OBJECT

public:

    /*! Constructor. */
  OneClick( void );
    /*! Destructor. */
  ~OneClick( void );

    /*! Read options, create stimulus and output of stimuli. */
  virtual int main( void );

  void saveSpikes( const string &file );

    /*! Plot data. */
  void plot( void );
    /*! Analyze data. */
  void analyze( void );


protected:

  double MinIntensity;
  double MaxIntensity;
  double IntensityStep;
  int IntRepeat;

  double Stepsize;
  double Duration;
  double Latency;
  double PreWidth;
  double Pause;
  int Side;

  double Intensity;

  MultiPlot P;  // Plots!

  Options Header; //?????????????

};


}; /* namespace auditoryprojects */

#endif /* ! _RELACS_AUDITORYPROJECTS_ONECLICK_H_ */
