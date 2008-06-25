/*
  tracespec.cc
  Specification of an output signal.

  RELACS - RealTime ELectrophysiological data Acquisition, Control, and Stimulation
  Copyright (C) 2002-2008 Jan Benda <j.benda@biologie.hu-berlin.de>

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

#include <relacs/outdata.h>
#include <relacs/tracespec.h>
using namespace std;

namespace relacs {


TraceSpec::TraceSpec( void )
  : Trace( -1 ),
    TraceName( "" ),
    Device( 0 ),
    Channel( 0 ),
    Scale( 1.0 ),
    Unit( "V" ),
    Reglitch( false ),
    MaxRate( 1000.0 ),
    SignalDelay( 0.0 )
{
}


TraceSpec::TraceSpec( int index, const string &name,
		      int device, int channel, 
		      double scale, const string &unit,
		      bool reglitch, double maxrate, double signaldelay )
  : Trace( index ),
    TraceName( name ),
    Device( device ),
    Channel( channel ),
    Scale( scale ),
    Unit( unit ),
    Reglitch( reglitch ),
    MaxRate( maxrate ),
    SignalDelay( signaldelay )
{
}


TraceSpec::TraceSpec( const TraceSpec &trace )
  : Trace( trace.Trace ),
    TraceName( trace.TraceName ),
    Device( trace.Device ),
    Channel( trace.Channel ),
    Scale( trace.Scale ),
    Unit( trace.Unit ),
    Reglitch( trace.Reglitch ),
    MaxRate( trace.MaxRate ),
    SignalDelay( trace.SignalDelay )
{
}


int TraceSpec::device( void ) const
{
  return Device;
}


void TraceSpec::setDevice( int device )
{
  Device = device;
}


int TraceSpec::channel( void ) const
{
  return Channel;
}


void TraceSpec::setChannel( int channel )
{
  Channel = channel;
}


void TraceSpec::setChannel( int channel, int device )
{
  Channel = channel;
  Device = device;
}


int TraceSpec::trace( void ) const
{
  return Trace;
}


void TraceSpec::setTrace( int index )
{
  Trace = index;
}


string TraceSpec::traceName( void ) const
{
  return TraceName;
}


void TraceSpec::setTraceName( const string &name )
{
  TraceName = name;
}


double TraceSpec::scale( void ) const
{
  return Scale;
}


void TraceSpec::setScale( double scale )
{
  Scale = scale;
}


string TraceSpec::unit( void ) const
{
  return Unit;
}


void TraceSpec::setUnit( const string &unit )
{
  Unit = unit;
}


void TraceSpec::setUnit( double scale, const string &unit )
{
  Scale = scale;
  Unit = unit;
}


bool TraceSpec::reglitch( void ) const
{
  return Reglitch;
}


void TraceSpec::setReglitch( bool reglitch )
{
  Reglitch = reglitch;
}


double TraceSpec::maxSampleRate( void )
{
  return MaxRate;
}


void TraceSpec::setMaxSampleRate( double maxrate )
{
  MaxRate = maxrate;
}


double TraceSpec::signalDelay( void ) const
{
  return SignalDelay;
}


void TraceSpec::setSignalDelay( double sigdelay )
{
  SignalDelay = sigdelay;
}


int TraceSpec::apply( OutData &signal ) const
{
  if ( ( ! signal.traceName().empty() && signal.traceName() == TraceName ) ||
       ( signal.trace() >= 0 && signal.trace() == Trace ) ) {
    signal.setDevice( Device );
    signal.setChannel( Channel );
    signal.setScale( Scale );
    signal.setUnit( Unit );
    signal.setReglitch( Reglitch );
    signal.setMaxSampleRate( MaxRate );
    signal.setSignalDelay( SignalDelay );
    return 0;
  }
  else {
    signal.addError( DaqError::InvalidTrace );
    return -1;
  }
}


bool operator==( const TraceSpec &trace, const OutData &signal )
{
  return ( trace.device() == signal.device() &&
	   trace.channel() == signal.channel() );
}


bool operator==( const OutData &signal, const TraceSpec &trace )
{
  return ( trace.device() == signal.device() &&
	   trace.channel() == signal.channel() );
}


}; /* namespace relacs */

