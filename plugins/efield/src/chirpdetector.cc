/*
  efield/chirpdetector.cc
  Detects chirps of wave-type weakly electric fish

  RELACS - Relaxed ELectrophysiological data Acquisition, Control, and Stimulation
  Copyright (C) 2002-2010 Jan Benda <benda@bio.lmu.de>

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

#include <relacs/efield/chirpdetector.h>
using namespace relacs;

namespace efield {


ChirpDetector::ChirpDetector( const string &ident, int mode )
  : Filter( ident, mode, SingleEventDetector, 1,
	    "ChirpDetector", "EField", 
	    "Jan Benda", "1.2", "Jun 17, 2009" )
{
  // parameter:
  Threshold = 8.0;
  MinThresh = 0.0;
  ChirpMinWidth = 0.003;
  ChirpMaxWidth = 0.05;
  ChirpCycles = 20;
  AverageCycles = 40;
  Other = 0;

  // options:
  addNumber( "threshold", "Threshold", Threshold, 0.0, 10000.0, 2.0, "Hz", "Hz", "%.0f", 2+8+32 );
  addNumber( "minthresh", "Minimum threshold", MinThresh, 0.0, 10000.0, 2.0, "Hz", "Hz", "%.0f", 2+8+32 );
  addNumber( "minwidth", "Minimum width", ChirpMinWidth, 0.0, 100.0, 0.002, "s", "ms", "%.0f", 2+8+32 );
  addNumber( "maxwidth", "Maximum width", ChirpMaxWidth, 0.002, 1000.0, 0.002, "s", "ms", "%.0f", 2+8+32 );
  addNumber( "rate", "Rate", 0.0, 0.0, 10000.0, 1.0, "Hz", "Hz", "%.0f", 2+4 );
  addNumber( "size", "Size", 0.0, 0.0, 10000.0, 1.0, "Hz", "Hz", "%.0f", 2+4 );
  addNumber( "width", "Width", 0.0, 0.0, 100000.0, 0.1, "ms", "ms", "%.0f", 2+4 );
  addStyle( OptWidget::ValueLarge + OptWidget::ValueBold + OptWidget::ValueGreen + OptWidget::ValueBackBlack, 4 );

  CDW.assign( ((Options*)this), 2, 4, true, 0, mutex() );
  CDW.setVerticalSpacing( 4 );
  CDW.setMargins( 4 );
  setWidget( &CDW );
  QObject::connect( (QWidget*)this, SIGNAL( dialogAccepted( void ) ),
		    &CDW, SLOT( updateValues( void ) ) );

  setDialogSelectMask( 8 );
  setDialogReadOnlyMask( 16 );
  setConfigSelectMask( -32 );
}


ChirpDetector::~ChirpDetector( void )
{
}


int ChirpDetector::init( const EventData &inevents, EventData &outevents,
			 const EventList &other, const EventData &stimuli )
{
  outevents.setSizeScale( 1.0 );
  outevents.setSizeUnit( "Hz" );
  outevents.setSizeFormat( "%6.1f" );
  outevents.setWidthScale( 1000.0 );
  outevents.setWidthUnit( "ms" );
  outevents.setWidthFormat( "%6.1f" );

  D.init( EventFrequencyIterator( inevents.begin() + 1 ),
	  EventFrequencyIterator( inevents.end() ),
	  EventIterator( inevents.begin() + 1 ) );

  Other = &other;

  return 0;
}


void ChirpDetector::notify( void )
{
  Threshold = number( "threshold" );
  MinThresh = number( "minthresh" );
  ChirpMinWidth = number( "minwidth" );
  ChirpMaxWidth = number( "maxwidth" );

  if ( Threshold < MinThresh ) {
    Threshold = MinThresh;
    setNumber( "threshold", Threshold );
  }
  CDW.updateValues( OptWidget::changedFlag() );
}


int ChirpDetector::detect( const EventData &inevents, EventData &outevents,
			    const EventList &other, const EventData &stimuli )
{
  long lastsize = outevents.size();

  D.peak( EventFrequencyIterator( inevents.begin() + 1 ), 
	  EventFrequencyIterator( inevents.end() ),
	  outevents, Threshold, MinThresh, 2.0*Threshold, *this );

  if ( outevents.size() - lastsize <= 0 )
    outevents.updateMean();

  unsetNotify();
  setNumber( "rate", outevents.meanRate() );
  setNumber( "size", outevents.meanSize() );
  setNumber( "width", outevents.meanWidth() );
  setNotify();
  CDW.updateValues( OptWidget::changedFlag() );
  return 0;
}


int ChirpDetector::checkEvent( const EventFrequencyIterator &first, 
			       const EventFrequencyIterator &last,
			       EventFrequencyIterator &event, 
			       EventIterator &eventtime, 
			       EventFrequencyIterator &index,
			       EventIterator &indextime,
			       EventFrequencyIterator &prevevent, 
			       EventIterator &prevtime, 
			       EventData &outevents, 
			       double &threshold,
			       double &minthresh, double &maxthresh,
			       double &time, double &size, double &width )
{
  // accept everything as a chirp:
  time = *eventtime;
  size = *event;
  width = *indextime - *eventtime;
  return 1;

  // XXX The following algorithm is quite broken: XXX

  // chirp too long:
  if ( *indextime - *eventtime > ChirpMaxWidth )
    return 0;

  // store event:
  EventFrequencyIterator orgevent = event;
  EventIterator orgtime = eventtime;

  // meanrate after chirp:
  double ameanrate = 0.0;
  EventFrequencyIterator cindex = event + ChirpCycles;
  for ( int j=0; j<AverageCycles; j++ ) {
    if ( cindex >= last )
      return -1;
    ameanrate += ( *cindex - ameanrate )/(j+1);
    ++cindex;
  }

  // meanrate before chirp:
  double meanrate = 0.0;
  cindex = event - ChirpCycles;
  for ( int j=0; j<AverageCycles && !cindex; j++ ) {
    meanrate += ( *cindex - meanrate )/(j+1);
    --cindex;
  }

  // eod detector was down:
  if ( ameanrate < 0.3*meanrate || meanrate < 0.3*ameanrate )
    return 0;
  for ( cindex = event - ChirpCycles; cindex != event + ChirpCycles && !cindex; ++cindex ) {
    if ( *cindex < 0.3*meanrate )
      return 0;
  }


  // find end of chirp (ChirpCycles data points within +- 1/2 threshold in a row):
  EventFrequencyIterator findex;
  EventIterator ftime;
  int n = 0;
  double rate = *index;
  ftime = indextime+1;
  for ( findex = index+1; !findex; ++findex, ++ftime ) {
    if ( *event < *findex )
      event = findex;
    if ( fabs( *findex - rate ) < 0.5*threshold ) 
      n++;
    else {
      n = 0;
      rate = *findex;
    }
    if ( n > ChirpCycles )
      break;
  }

  // end of chirp not contained in data:
  if ( findex >= last || *findex > rate + 0.5*threshold )
    return -1;

  // size of chirp:
  size = *event - meanrate;
  if ( size < MinThresh )
    return 0;

  // chirp:
  time = *eventtime;
  double minrate = meanrate + 0.1 * size;
  index = findex;
  indextime = ftime;

  // find begin of chirp:
  EventFrequencyIterator lindex;
  EventIterator ltime;
  if ( *orgevent > minrate ) {
    ltime = orgtime-1;
    for ( lindex = orgevent-1; !lindex; --lindex, --ltime )
      if ( *lindex <= minrate )
	break;
  }
  else {
    ltime = orgtime+1;
    for ( lindex = orgevent+1; !lindex; ++lindex, --ltime )
      if ( *lindex >= minrate )
	break;
  }

  // lindex not valid:
  if ( !(!lindex) )
    return -1;

  // find end of chirp:
  for ( --findex, --ftime; !findex && findex > event; --findex, --ftime )
    if ( *findex >= minrate )
      break;

  width = *ftime - *ltime;

  // chirp too short:
  if ( width < ChirpMinWidth )
    return 0;

  // chirp too long:
  if ( width > ChirpMaxWidth )
    return 0;

  // check for chirps in the other event traces:
  for ( int k=0; k<Other->size(); k++ )
    if ( (*Other)[k].within( time, width ) )
      return 0;

  // this is a chirp which occurs only in this trace:
  return 1;
}


addDetector( ChirpDetector );

}; /* namespace efield */

#include "moc_chirpdetector.cc"
