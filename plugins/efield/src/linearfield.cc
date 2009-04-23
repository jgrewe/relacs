/*
  efield/linearfield.cc
  Measure the electric field manually with a single electrode in one direction

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

#include <qvbox.h>
#include <qlabel.h>
#include <qpushbutton.h>
#include <relacs/efield/linearfield.h>
using namespace relacs;

namespace efield {


LinearField::LinearField( void )
  : RePro( "LinearField", "LinearField", "Efield",
	   "Jan Benda", "1.0", "Apr 23, 2009" ),
    P( (QWidget*)this )
{
  // add some options:

  QVBox *bb = new QVBox( this );
  bb->setSpacing( 4 );

  // distance input:
  QLabel *l = new QLabel( "Distance", bb );
  DSB = new DoubleSpinBox( 0.0, -10000.0, 10000.0, 1.0, "%g", bb );

  // Measure button:
  QPushButton *measurebutton = new QPushButton( "&Measure", bb, "MeasureButton" );
  connect( measurebutton, SIGNAL( clicked() ),
	   this, SLOT( measure() ) );

  // Quit button:
  QPushButton *quitbutton = new QPushButton( "&Quit", bb, "QuitButton" );
  connect( quitbutton, SIGNAL( clicked() ),
	   this, SLOT( quit() ) );
}


void LinearField::measure( void )
{
  Measure = true;
  wake();
}


void LinearField::quit( void )
{
  Measure = false;
  wake();
}


int LinearField::main( void )
{
  // get options:

  noMessage();

  postCustomEvent( 1 ); // DSB->setFocus();
  do {
    // wait for input:
    Measure = false;
    unlockAll();
    sleepWait();
    lockAll();
    if ( Measure ) {
      message( "measure" );
      // analyse:
      /*
      OutList sigs;
      for ( int k=0; k<OutOpts.size(); k++ ) {
	if ( OutOpts[k].changed() ) {
	  double value = OutOpts[k].number();
	  OutData sig( value );
	  sig.setTraceName( OutOpts[k].ident() );
	  sig.setIdent( "value=" + Str( value ) );
	  sigs.push( sig );
	}
      }
      if ( sigs.size() > 0 ) {
	directWrite( sigs );
	if ( sigs.failed() ) {
	  warning( sigs.errorText() );
	  return Failed;
	}
      }
      */
    }
  } while ( Measure );
  postCustomEvent( 2 ); // DSB->clearFocus();
  return Completed;
}


void LinearField::customEvent( QCustomEvent *qce )
{
  if ( qce->type() == QEvent::User+1 ) {
    DSB->setFocus();
  }
  else if ( qce->type() == QEvent::User+2 ) {
    DSB->clearFocus();
  }
}


addRePro( LinearField );

}; /* namespace efield */

#include "moc_linearfield.cc"

