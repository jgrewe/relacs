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
  newSection( "Stimulation" );
  addText( "name", "Name/role of the repro", "" );
  addNumber( "duration", "Stimulus duration", 1.0, 0.001, 100000.0, 0.001, "s", "ms" );
  addNumber( "cutoff", "Cutoff frequency of noise stimulus", 150.0, 1.0, 2000.0, 1.0, "Hz" );
  addNumber( "contrast", "Contrast of the stimulus, i.e. stimulus intensity", 10.0, 0.1, 100.0, .1, "%");
  addNumber( "samplerate", "temporal resolution of the stimulus", 1000.0, 100.0, 100000.0, 100, "Hz");
  addNumber( "count", "Stimulus repeats", 1, 1, 100, 1);
  addNumber( "pause", "Time between successive stimulus presentations", 0.5, 0.0, 100.0, 0.01, "s" ).setActivation( "count", ">1" );

  newSection( "Analysis" );
  addNumber( "tmin", "Time before the spike time", -0.02, -0.001, -0.15, 0.001, "s", "ms" );
  addNumber( "tmax", "Time after the spike time", 0.02, 0.001, 0.15, 0.001, "s", "ms" );
  addBoolean( "reconstruct", "Do the stimulus reconstruction", false );
  addBoolean( "psth", "Plot the firing rate.", false );
  addNumber( "kernel", "Kernel width for firing rate estimation.", 0.001, 0.00001, 0.1, 0.0005, "s" ).setActivation( "psth", "true" );

  kernel = GaussKernel(0.001);
  QVBoxLayout *vb = new QVBoxLayout;
  QHBoxLayout *hb = new QHBoxLayout;

  startTime = 0.0;
  duration = 0.0;

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
  spikesPlot.setY2Range( 0., 1.5 );
  spikesPlot.setXRange( 0., 1.0 );
  spikesPlot.setY2Label( "trials" );
  spikesPlot.setLMarg( 6 );
  spikesPlot.setRMarg( 1 );
  spikesPlot.setTMarg( 3 );
  spikesPlot.setBMarg( 4 );
  spikesPlot.unlock();

  staPlot.lock();
  staPlot.setXLabel( "time [s]" );
  staPlot.setXRange( -0.02, 0.02 );
  staPlot.unlock();

  vb->addWidget( &stimPlot );
  vb->addWidget( &spikesPlot );
  hb->addLayout( vb );
  hb->addWidget( &staPlot );
  setLayout( hb );
}

void SpikeTriggeredAverage::plotStimulus(const OutData &stimulus) {
  stimPlot.lock();
  stimPlot.clear();
  stimPlot.setTitle( "stimulus" );
  stimPlot.plot( stimulus, 1.0, Plot::Yellow, 2, Plot::Solid );
  stimPlot.setAutoScaleY();
  stimPlot.unlock();
}

void SpikeTriggeredAverage::analyze( EventList &myspikes, int currentRepeat ) {
  int spikes_trace = 0;
  for ( int k=0; k<MaxTraces; k++ ) {
    if ( SpikeEvents[k] >= 0 ) {
      spikes_trace = k;
      break;
    }
  }

  const EventData &spikeEvents = events( SpikeEvents[spikes_trace] );
  SampleDataD firingRate( int(duration * samplerate), 0.0, 1.0/samplerate );
  SampleDataD rateSd( firingRate );
  EventData sp( spikeEvents, startTime, startTime+duration, startTime );
  SampleDataD temp(sta);
  SampleDataD tempSta(sta);
  int count = 0, firstIndex = 0, lastIndex = 0, spikeIndex =0;
  SampleDataD revRec( stimCopy.size(), 0.0, stimCopy.stepsize(), 0.0 );

  myspikes[currentRepeat].clear();
  myspikes[currentRepeat].append( spikeEvents, startTime, startTime+duration, startTime );
  spikesPlot.lock();
  spikesPlot.setY2Range( 0.0, double(currentRepeat + 1.25) );
  for (int i = 0; i < myspikes[currentRepeat].size(); ++i ) {
    double x = myspikes[currentRepeat][i];
    spikesPlot.plotPoint(x, Plot::SecondX, double( currentRepeat + 1 ), Plot::SecondY, 0,
                         Plot::StrokeVertical, 0.5, Plot::SecondY, Plot::Red, Plot::Red);
    // sta estimation
    if (x - std::abs(tmin) >= 0.0 && x + tmax <= duration) {
      stimCopy.copy(x - std::abs(tmin), x + tmax, temp);
      tempSta += temp;
      count++;
    }
  }
  tempSta /= count;
  if ( currentRepeat == 0 ) {
    sta = tempSta;
  } else {
    sta += ( tempSta - sta )/( currentRepeat + 1 );
  }
  sta.setOffset(tmin);

  if ( plotPsth ) {
    if ( psthIndex >= 0 )
      spikesPlot.clearData( psthIndex );
    else {
      spikesPlot.setYLabel("firing rate [Hz]");
    }
    spikesPlot.setAutoScaleY();
    myspikes.rate( firingRate, rateSd, kernel );
    Plot::LineStyle l( Plot::Yellow );
    psthIndex = spikesPlot.plot( firingRate, 1.0, Plot::Green, 1 );
  }
  spikesPlot.draw();
  spikesPlot.unlock();

  staPlot.lock();
  staPlot.clear();
  staPlot.plot( sta, 1.0, Plot::Green, 2. );
  staPlot.setXRange( tmin, tmax );
  staPlot.setAutoScaleY();
  staPlot.draw();
  staPlot.unlock();

  if ( reconstruct ) {
    for (int i = 0; i < myspikes[currentRepeat].size(); ++i ) {
      double x = myspikes[currentRepeat][i];
      if ( reconstruct ) {
        spikeIndex = stimCopy.index( x );
        firstIndex = stimCopy.index( x - std::abs( tmin ) );
        lastIndex = stimCopy.index( x + tmax );
        if ( firstIndex >= 0  && lastIndex < stimCopy.size() ){
          for ( int index = 0; index < sta.size(); ++index ){
            revRec[firstIndex + index] += sta[index];
          }
        }
      }
    }
    if ( currentRepeat == 0 )
      stimReconstruct = revRec;
    else
      stimReconstruct += ( revRec - stimReconstruct )/( currentRepeat + 1);
    stimPlot.lock();
    if ( reconstructionIndex >= 0 )
      stimPlot.clearData( reconstructionIndex );
    reconstructionIndex = stimPlot.plot( stimReconstruct, 1.0, Plot::Red, 1 );
    stimPlot.draw();
    stimPlot.unlock();
  }
}


int SpikeTriggeredAverage::main( void )
{
  duration = number( "duration" );
  int count = number( "count" );
  double pause = number( "pause" );
  double cutoff = number( "cutoff" );
  double contrast = number( "contrast" );
  tmax = number( "tmax" );
  tmin = number( "tmin" );
  samplerate = number( "samplerate" );
  double eod = eodAmplitude( trace( LocalEODTrace[0] ),
                             currentTime() - 0.2, currentTime() );
  plotPsth = boolean( "psth" );
  reconstruct = boolean( "reconstruct" );
  psthIndex = -1;
  reconstructionIndex = -1;
  if ( plotPsth )
    kernel = GaussKernel( number( "kernel" ) );

  EventList spikes( count );
  OutData stimulus;
  stimulus.setTrace( GlobalAMEField );
  stimulus.setStepsize( 1.0/samplerate );
  stimulus.bandNoiseWave( duration, stimulus.stepsize(), 0.0, cutoff, contrast/100.0 );
  stimulus.setIntensity(1.0);  // FIXME amplitude is independent of Fish, so far
  stimCopy.resize( stimulus.size(), 0.0 );
  for ( int i = 0; i < stimulus.size(); ++i ) {
    stimCopy[i] = stimulus[i];
  }
  stimCopy.setStepsize( stimulus.stepsize() );
  stimReconstruct.resize( stimulus.size(), 0.0, stimulus.stepsize() );
  plotStimulus( stimulus );

  spikesPlot.clear();
  staPlot.clear();
  SampleDataD sta( int((tmax - tmin) * samplerate), tmin, 1.0/samplerate, 0.0);

  for (int i = 0; i < count; ++i ) {
    startTime = currentTime();
    write( stimulus );
    if ( !stimulus.success() ) {
      string s = "Output of stimulus failed!<br>Error is <b>";
      s += stimulus.errorText() + "</b>";
      warning( s );
      //stop();
      return Failed;
    }
    analyze( spikes, i );
    sleep( pause );
  }
  return Completed;
}


addRePro( SpikeTriggeredAverage, efish );

}; /* namespace efish */

#include "moc_spiketriggeredaverage.cc"
