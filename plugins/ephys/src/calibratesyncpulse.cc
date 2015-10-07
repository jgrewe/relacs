/*
  ephys/calibratesyncpulse.cc
  Sets the width of the pulse for synchronizing an amplifier.

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

#include <relacs/rangeloop.h>
#include <relacs/map.h>
#include <relacs/ephys/amplifiercontrol.h>
#include <relacs/ephys/calibratesyncpulse.h>
using namespace relacs;

namespace ephys {


CalibrateSyncPulse::CalibrateSyncPulse( void )
  : RePro( "CalibrateSyncPulse", "ephys", "Jan Benda", "1.0", "Oct 07, 2015" ),
    IUnit( "nA" )
{
  // add some options:
  addNumber( "imin", "Minimum injected current", -1.0, -1000.0, 1000.0, 0.001 );
  addNumber( "imax", "Maximum injected current", 1.0, -1000.0, 1000.0, 0.001 );
  addNumber( "istep", "Step-size of current increments", 0.001, 0.001, 1000.0, 0.001 );
  addNumber( "skipwin", "Initial portion of the response not used for analysis", 1.0, 0.001, 100000.0, 0.001, "s", "ms" );
  addNumber( "duration", "Stimulus duration used for analysis", 1.0, 0.001, 100000.0, 0.001, "s", "ms" );
}


void CalibrateSyncPulse::preConfig( void )
{
  if ( CurrentOutput[0] >= 0 ) {
    IUnit = outTrace( CurrentOutput[0] ).unit();
    setUnit( "imin", IUnit );
    setUnit( "imax", IUnit );
    setUnit( "istep", IUnit );
  }
}


int CalibrateSyncPulse::main( void )
{
  // get options:
  double imin = number( "duration" );
  double imax = number( "duration" );
  double istep = number( "duration" );
  double skipwin = number( "skipwin" );
  double duration = number( "duration" );

  // in- and outtrace:
  const InData &intrace = trace( SpikeTrace[0] >= 0 ? SpikeTrace[0] : 0 );
  int outtrace = CurrentOutput[0] >= 0 ? CurrentOutput[0] : 0;

  lockStimulusData();
  double dccurrent = stimulusData().number( outTraceName( outtrace ) );
  unlockStimulusData();

  // get amplifier:
  AmplifierControl *ampl = dynamic_cast< AmplifierControl* >( control( "AmplifierControl" ) );
  if ( ampl == 0 ) {
    warning( "No amplifier found." );
    return Failed;
  }

  // dc signal:
  OutData dcsignal;
  dcsignal.setTrace( outtrace );
  dcsignal.constWave( dccurrent );
  dcsignal.setIdent( "DC=" + Str( dccurrent ) + IUnit );

  MapD volt;

  // write stimulus:
  RangeLoop range( imin, imax, istep );
  for ( range.reset(); ! range && softStop() == 0; ) {

    ampl->activateCurrentClampMode();

    double amplitude = *range;
    if ( fabs( amplitude ) < 1.0e-8 )
      amplitude = 0.0;

    Str s = "Current <b>" + Str( amplitude ) + " " + IUnit +"</b>";
    message( s );

    // signal:
    OutData signal;
    signal.setTrace( outtrace );
    signal.constWave( amplitude );
    signal.setIdent( "DC=" + Str( amplitude ) + IUnit );
    write( signal );
    if ( signal.failed() ) {
      if ( signal.overflow() ) {
	printlog( "Requested amplitude I=" + Str( amplitude ) + IUnit + "too high!" );
	break;
      }
      else if ( signal.underflow() ) {
	printlog( "Requested amplitude I=" + Str( amplitude ) + IUnit + "too small!" );
	continue;
      }
      else {
	warning( signal.errorText() );
	directWrite( dcsignal );
	return Failed;
      }
    }
    sleep( skipwin + duration );
    double ccvolt = intrace.mean( currentTime()-duration, currentTime() );

    ampl->activateDynamicClampMode();
    sleep( skipwin + duration );
    double dcvolt = intrace.mean( currentTime()-duration, currentTime() );

    volt.push( ccvolt, dcvolt );

    // plot

  }

  directWrite( dcsignal );
  sleep( duration );

  // analyze and save

  return Completed;
}


addRePro( CalibrateSyncPulse, ephys );

}; /* namespace ephys */

#include "moc_calibratesyncpulse.cc"
