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
\brief [RePro] Record the spike triggered average using frozen white noise.
\author Jan Grewe
\version 1.0 (Jul 17, 2018)
\par Options
- \c Stimulation
    - \c name="": The name used during storing of the data (\c string)
    - \c duration=1seconds: Signal duration (\c number)
    - \c cutoff=150hertz: Upper cutoff frequency of the stimulus (\c number)
    - \c contrast=10percent: the intensity of the stimulus relative to EOD amplitude (\c number)
    - \c samplerate=1000Hertz: temporal resolution of the stimulus (\c number)
    - \c count=5: number of stimulus repetitions (\c number)
    - \c pause=0.2s: duration of the pause between stimulus presenstations (\c number)
- \c Analysis
    - \c tmin=-0.02s: minimum time before the spike to use for the sta (\c number)
    - \c tmax=0.02s: maximum time after the spike to use for the sta (\c number)
    - \c reconstruct=false: defines whether or not a reverse reconstruction of the stimulus is done (\c boolean)
    - \c psth=false: defines whether the firing rate is calculated and plotted (\c boolean)
    - \c kernel=0.001s: width of the Gaussian kernel used for firing rate calculation (\c number)
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
  SampleDataD stimCopy, stimReconstruct;
  SampleDataD sta;
  double startTime, duration, samplerate;
  double tmin, tmax;
  bool plotPsth, reconstruct;
  int psthIndex = -1, reconstructionIndex = -1;
  GaussKernel kernel;
  void plotStimulus( const OutData &stimulus );
  void analyze( EventList &spikes, int currentRepeat );
};


}; /* namespace efish */

#endif /* ! _RELACS_EFISH_SPIKETRIGGEREDAVERAGE_H_ */
