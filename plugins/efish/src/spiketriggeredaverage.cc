/*
  efish/spiketriggeredaverage.cc
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

#include <relacs/efish/spiketriggeredaverage.h>
using namespace relacs;

namespace efish {


SpikeTriggeredAverage::SpikeTriggeredAverage( void )
  : RePro( "SpikeTriggeredAverage", "efish", "Jan Grewe", "1.0", "Jul 17, 2018" )
{
  newSection( "Stimulus" );
  addNumber( "duration", "Stimulus duration", 1.0, 0.001, 100000.0, 0.001, "s", "ms" );
  addNumber( "cutoff", "Cutoff frequency of noise stimulus", 150.0, 1.0, 2000.0, 1.0, "Hz" );
  addNumber( "contrast", "Contrast of the stimulus, i.e. stimulus intensity", 10.0, 0.1, 100.0, .1, "%");
  addNumber( "samplerate", "temporal resolution of the stimulus", 1000.0, 100.0, 100000.0, 100, "Hz");
  addNumber( "count", "Stimulus repeats", 1, 1, 100, 1);
  addNumber( "pause", "Time between successive stimulus presentations", 0.5, 0.0, 100.0, 0.01, "s" ).setActivation( "count", ">1" );

  newSection( "Analysis" );
  addNumber( "tmin", "Time before the spike time", 0.02, 0.001, 0.15, 0.001, "s", "ms" );
  addNumber( "tmax", "Time after the spike time", 0.02, 0.001, 0.15, 0.001, "s", "ms" );
  addBoolean( "reconstruct", "Do the stimulus reconstruction", false );
  QVBoxLayout *vb = new QVBoxLayout;
  QHBoxLayout *hb = new QHBoxLayout;
  stimPlot.lock();
  stimPlot.setXLabel( "time [s]" );
  stimPlot.setYRange( -2.0, 2.0 );
  stimPlot.setXRange( 0., 1.0 );
  stimPlot.setYLabel( "voltage [mV]" );
  stimPlot.setLMarg( 6 );
  stimPlot.setRMarg( 1 );
  stimPlot.setTMarg( 3 );
  stimPlot.setBMarg( 4 );
  stimPlot.unlock();

  spikesPlot.lock();
  spikesPlot.setXLabel( "time [s]" );
  spikesPlot.setYRange( 0., 1.5 );
  spikesPlot.setXRange( 0., 1.0 );
  spikesPlot.setYLabel( "trials" );
  spikesPlot.setLMarg( 6 );
  spikesPlot.setRMarg( 1 );
  spikesPlot.setTMarg( 3 );
  spikesPlot.setBMarg( 4 );
  spikesPlot.unlock();
  vb->addWidget( &stimPlot );
  vb->addWidget( &spikesPlot );
  hb->addLayout( vb );
  hb->addWidget( &staPlot );
  setLayout( hb );
}


int SpikeTriggeredAverage::main( void )
{
  double duration = number( "duration" );
  int count = number( "count" );
  double pause = number( "pause" );
  double cutoff = number( "cutoff" );
  double contrast = number( "contrast" );
  double tmax = number( "tmax" );
  double tmin = number( "tmin" );
  double samplerate = number( "samplerate" );
  double eod = eodAmplitude( trace( LocalEODTrace[0] ),
                             currentTime() - 0.2, currentTime() );

  OutData stimulus;
  stimulus.setTrace( GlobalAMEField );
  stimulus.setStepsize( 1.0/samplerate );
  stimulus.bandNoiseWave( duration, stimulus.stepsize(), 0.0, cutoff, contrast/100.0 );
  stimulus.setIntensity(1.0);

  plotStimulus();
  for (int i = 0; i < count; ++i ) {
    write( stimulus );
    if ( !stimulus.success() ) {
      string s = "Output of stimulus failed!<br>Error is <b>";
      s += stimulus.errorText() + "</b>";
      warning( s );
      //stop();
      return Failed;
    }
   sleep( pause ); 
  }
  return Completed;
}


addRePro( SpikeTriggeredAverage, efish );

}; /* namespace efish */

#include "moc_spiketriggeredaverage.cc"
