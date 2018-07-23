/*
  efish/spiketriggeredaverage.h
  Record a spike triggered average using white noise stimuli.

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

#ifndef _RELACS_EFISH_SPIKETRIGGEREDAVERAGE_H_
#define _RELACS_EFISH_SPIKETRIGGEREDAVERAGE_H_ 1

#include <relacs/repro.h>
#include <relacs/ephys/traces.h>
#include <relacs/efield/traces.h>
#include <relacs/efield/eodtools.h>
#include <relacs/plot.h>
#include <relacs/kernel.h>

using namespace relacs;

namespace efish {


/*!
\class SpikeTriggeredAverage
\brief [RePro] Record a spike triggered average using white noise stimuli.
\author Jan Grewe
\version 1.0 (Jul 17, 2018)
*/


class SpikeTriggeredAverage
  : public RePro,
  public ephys::Traces,
  public efield::Traces,
  public efield::EODTools

{
  Q_OBJECT

public:

  SpikeTriggeredAverage( void );
  virtual int main( void );

 private:
  Plot stimPlot;
  Plot spikesPlot;
  Plot staPlot;
  SampleDataD stimCopy;
  double startTime, duration, samplerate;
  double tmin, tmax;
  bool plotPsth;
  int psthIndex = -1;
  GaussKernel kernel;
  void plotStimulus( const OutData &stimulus );
  void analyze( EventList &spikes, int currentRepeat );
};


}; /* namespace efish */

#endif /* ! _RELACS_EFISH_SPIKETRIGGEREDAVERAGE_H_ */
