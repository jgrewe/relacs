/*
  plottrace.cc
  Plot traces and events.

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

#include <iostream>
#include <fstream>
#include <QPixmap>
#include <QBitmap>
#include <QPainter>
#include <QPolygon>
#include <QDateTime>
#include <QLabel>
#include <QToolTip>
#include <QApplication>
#include <relacs/str.h>
#include <relacs/datafile.h>
#include <relacs/filterdetectors.h>
#include <relacs/relacswidget.h>
#include <relacs/plottrace.h>

namespace relacs {


PlotTrace::PlotTrace( RELACSWidget *rw, QWidget* parent )
  : RELACSPlugin( "PlotTrace: PlotTrace", RELACSPlugin::Plugins,
		  "PlotTrace", "base" ),
    Menu( 0 ),
    PlotTimer( this ),
    P( 1, Plot::Pointer, parent ),
    PlotElements( 1, -1 )
{
  setRELACSWidget( rw );

  ViewMode = SignalView;
  ContinuousView = WrapView;
  setView( WrapView );
  Manual = false;
  Trigger = true;
  TriggerSource = -1;
  Plotting = true;

  connect( &PlotTimer, SIGNAL( timeout() ),
	   this, SLOT( plot() ) );
  connect( &P, SIGNAL( changedRanges( int ) ),
	   this, SLOT( updateRanges( int ) ) );
  connect( &P, SIGNAL( resizePlots( QResizeEvent* ) ),
	   this, SLOT( resizePlots( QResizeEvent* ) ) );

  TimeWindow = 1.0;
  TimeOffs = 0.0;
  Offset = 0.0;

  AutoOn = true;
  AutoFixed = false;
  AutoTime = 0.1;
  AutoOffs = 0.0;

  ButtonBox = 0;
  ViewButton = 0;
  ManualButton = 0;
  OnOffButton = 0;

  FilePlot = false;
  FilePath = "";
  FileTraces.clear();
  FileTracesNames.clear();
  FileSizes.clear();
  FileEvents.clear();
  FileEventsNames.clear();

  PlotTraces.clear();
  PlotEvents.clear();

  setWidget( &P );

  VP.clear();

  ButtonBox = 0;
  /*
  ButtonBox = new QWidget( this );
  ButtonBoxLayout = new QHBoxLayout;
  ButtonBoxLayout->setContentsMargins( 0, 0, 0, 0 );
  ButtonBoxLayout->setSpacing( 0 );
  ButtonBox->setLayout( ButtonBoxLayout );
  */

  int s = P.fontInfo().pixelSize();

  SignalViewIcon = QPixmap( s, s );
  QPainter p;
  p.begin( &SignalViewIcon );
  p.eraseRect( SignalViewIcon.rect() );
  p.setPen( QPen() );
  p.setBrush( Qt::black );
  QPolygon pa( 3 );
  pa.setPoint( 0, s-3, 2 );
  pa.setPoint( 1, s-3, s-2 );
  pa.setPoint( 2, 4, s/2 );
  p.drawPolygon( pa );
  p.setPen( QPen( Qt::black, 2 ) );
  p.drawLine( 3, 2, 3, s-1 );
  p.end();
  SignalViewIcon.setMask( SignalViewIcon.createHeuristicMask() );

  EndViewIcon = QPixmap( s, s );
  p.begin( &EndViewIcon );
  p.eraseRect( EndViewIcon.rect() );
  p.setPen( QPen() );
  p.setBrush( Qt::black );
  pa.setPoint( 0, 3, 2 );
  pa.setPoint( 1, 3, s-2 );
  pa.setPoint( 2, s-4, s/2 );
  p.drawPolygon( pa );
  p.setPen( QPen( Qt::black, 2 ) );
  p.drawLine( s-2, 2, s-2, s-1 );
  p.end();
  EndViewIcon.setMask( EndViewIcon.createHeuristicMask() );

  /*
  ViewButton = new QPushButton;
  ButtonBoxLayout->addWidget( ViewButton );
  ViewButton->setIcon( SignalViewIcon );
  ViewButton->setToolTip( "F: fixed (Pos1), C: continous (End)" );
  connect( ViewButton, SIGNAL( clicked() ), this, SLOT( viewToggle() ) );
  */

  QBitmap manualmask( s, s );
  p.begin( &manualmask );
  p.setPen( QPen( Qt::color1 ) );
  p.setBrush( QBrush() );
  p.setFont( QFont( "Helvetica", s, QFont::Bold ) );
  p.drawText( manualmask.rect(), Qt::AlignCenter, "M" );
  p.end();
  QPixmap manualicon( s, s );
  p.begin( &manualicon );
  p.eraseRect( manualicon.rect() );
  p.setPen( QPen( Qt::black ) );
  p.setBrush( QBrush() );
  p.setFont( QFont( "Helvetica", s, QFont::Bold ) );
  p.drawText( manualicon.rect(), Qt::AlignCenter, "M" );
  p.end();
  manualicon.setMask( manualmask );

  /*
  ManualButton = new QPushButton;
  ButtonBoxLayout->addWidget( ManualButton );
  ManualButton->setCheckable( true );
  ManualButton->setIcon( manualicon );
  ManualButton->setDown( Manual );
  ManualButton->setToolTip( "Manual or Auto" );
  connect( ManualButton, SIGNAL( clicked() ), this, SLOT( toggleManual() ) );
  */

  /*
  QPixmap onofficon( s, s );
  p.begin( &onofficon, this );
  p.eraseRect( onofficon.rect() );
  p.setPen( QPen( black, 2 ) );
  p.setBrush( QBrush() );
  p.drawArc( 2, 2, s-4, s-4, 1920, 4800 );
  p.drawLine( s/2, 0, s/2, s/2 );
  p.end();
  onofficon.setMask( onofficon.createHeuristicMask() );

  OnOffButton = new QPushButton;
  ButonBoxLayout->addWidget( OnOffButton );
  OnOffButton->setToggleButton( true );
  OnOffButton->setPixmap( onofficon );
  OnOffButton->setOn( ! Plotting );
  QToolTip::add( OnOffButton, "Switch the plot on or off" );
  connect( OnOffButton, SIGNAL( clicked() ), this, SLOT( plotOnOff() ) );
  */
}


PlotTrace::~PlotTrace( void )
{
}


void PlotTrace::resize( void )
{
  P.lock();

  // count active plots:
  VP.clear();
  for ( unsigned int c=0; c<PlotProps.size(); c++ ) {
    if ( PlotProps[c].Visible )
      VP.push_back( c );
  }

  // setup plots:
  if ( PlotTraces.size() != P.size() )
    P.resize( PlotTraces.size(), Plot::Pointer );
  P.setDataMutex( mutex() );
  P.setCommonXRange();
  PlotElements.clear();
  PlotElements.resize( VP.size(), -1 );

  if ( PlotTraces.size() > 0 ) {

    for ( int c=0; c<PlotTraces.size(); c++ ) {
      P[c].clear();
      P[c].setSkip( ! PlotProps[c].Visible );
    }
    
    double lmarg = 11.0;
    if ( VP.size() == 1 )
      lmarg = 8.0;
    int o = 0;
    if ( VP.size() > 1 && VP.size() % 2 == 1 )
      o = 1;
    
    for ( unsigned int c=0; c<VP.size(); c++ ) {
      P[VP[c]].setLMarg( lmarg );
      P[VP[c]].setRMarg( 2.0 );
      P[VP[c]].setTMarg( 0.2 );
      P[VP[c]].setBMarg( 0.2 );
      P[VP[c]].noXTics();
      P[VP[c]].setXLabel( "" );
      P[VP[c]].setYTics();
      P[VP[c]].setYLabelPos( 2.0 + ((c+o)%2)*3.0, Plot::FirstMargin,
			     0.5, Plot::Graph, Plot::Center, -90.0 );
    }
    
    if ( VP.size() > 0 ) {
      P[VP.front()].setTMarg( 1.0 );
      P[VP.back()].setXTics();
      P[VP.back()].setXLabel( "msec" );
      P[VP.back()].setXLabelPos( 1.0, Plot::FirstMargin, 0.0, Plot::FirstAxis, 
				 Plot::Left, 0.0 );
      P[VP.back()].setBMarg( 2.5 );
    }
    
    if ( VP.size() > 6 ) {
      P[VP[(VP.size()+1)/2]].setTMarg( 1.0 );
      P[VP[(VP.size()-1)/2]].setXTics();
      P[VP[(VP.size()-1)/2]].setXLabel( "msec" );
      P[VP[(VP.size()-1)/2]].setXLabelPos( 1.0, Plot::FirstMargin,
					   0.0, Plot::FirstAxis, Plot::Left, 0.0 );
      P[VP[(VP.size()-1)/2]].setBMarg( 2.5 );
    }
    
    resizeLayout();
  }
  
  P.unlock();
  
  PlotChanged = true;
}


void PlotTrace::toggle( QAction *mtrace )
{
  // check for valid trace:
  bool nodata = true;
  unsigned int i=0;
  for ( i=0; i<PlotProps.size(); i++ ) {
    if ( PlotProps[i].Action == mtrace ) {
      nodata = false;
      break;
    }
  }
  if ( nodata )
    return;

  if ( PlotProps[i].Visible ) {
    for ( unsigned int k=0; k<PlotProps.size(); k++ ) {
      if ( k != i && PlotProps[k].Visible ) {
	PlotProps[i].Visible = false;
	break;
      }
    }
  }
  else
    PlotProps[i].Visible = true;
  PlotProps[i].Action->setChecked( PlotProps[i].Visible );
  resize();
  plot();
}


void PlotTrace::init( void )
{
  P.lock();

  int origin = ViewMode == FixedView ? 3 : 2;
  double tfac = 1000.0;
  string tunit = "ms";

  for ( unsigned int c=0; c<VP.size(); c++ ) {

    int vp = VP[c];

    // clear plot:
    P[vp].clear();

    // y-label:
    string s = PlotTraces[vp].ident() + " [" + PlotTraces[vp].unit() + "]";
    P[vp].setYLabel( s );
	
    // plot stimulus events:
    for ( int s=0; s<PlotEvents.size(); s++ ) {
      if ( (PlotEvents[s].mode() & PlotTraceMode) &&
	   (PlotEvents[s].mode() & StimulusEventMode) ) {
	P[vp].plot( PlotEvents[s], origin, Offset, tfac,
		    0.0, Plot::Graph, 2,
		    Plot::StrokeUp, 1.0, Plot::GraphY,
		    Plot::White );
	break;
      }
    }
    // plot restart events:
    for ( int s=0; s<PlotEvents.size(); s++ ) {
      if ( (PlotEvents[s].mode() & PlotTraceMode) &&
	   (PlotEvents[s].mode() & RestartEventMode) ) {
	P[vp].plot( PlotEvents[s], origin, Offset, tfac,
		    1.0, Plot::Graph, 1,
		    Plot::TriangleNorth, 0.07, Plot::GraphY,
		    Plot::Orange, Plot::Orange );
	break;
      }
    }
    // plot recording events:
    for ( int s=0; s<PlotEvents.size(); s++ ) {
      if ( (PlotEvents[s].mode() & PlotTraceMode) &&
	   (PlotEvents[s].mode() & RecordingEventMode) ) {
	P[vp].plot( PlotEvents[s], origin, Offset, tfac,
		    0.0, Plot::Graph, 4,
		    Plot::StrokeUp, 1.0, Plot::GraphY,
		    Plot::Red );
	break;
      }
    }
    // plot events:
    int sn = 0;
    for ( int s=0; s<PlotEvents.size(); s++ ) {
      if ( (PlotEvents[s].mode() & PlotTraceMode) &&
	   !(PlotEvents[s].mode() & StimulusEventMode) &&
	   !(PlotEvents[s].mode() & RestartEventMode) &&
	   !(PlotEvents[s].mode() & RecordingEventMode) ) {
	
	if ( RW->FD->eventInputTrace( s ) == int( vp ) ) {

	  if ( sn == 0 ) {
	    /*
	      P[vp].plot( PlotEvents[s], origin, Offset, tfac,
	      0.05, Plot::Graph, 1,
	      Plot::StrokeUp, 20, Plot::Pixel,
	      Plot::Red );
	    */
	    P[vp].plot( PlotEvents[s], PlotTraces[vp],
			origin, Offset, tfac,
			1, Plot::Circle, 6, Plot::Pixel,
			Plot::Gold, Plot::Gold );
	  }
	  else if ( sn == 1 )
	    P[vp].plot( PlotEvents[s], origin, Offset, tfac,
			0.1, Plot::Graph, 1,
			Plot::Circle, 6, Plot::Pixel,
			Plot::Yellow, Plot::Yellow );
	  else if ( sn == 2 )
	    P[vp].plot( PlotEvents[s], origin, Offset, tfac, 0.2,
			Plot::Graph, 1,
			Plot::Diamond, 6, Plot::Pixel,
			Plot::Blue, Plot::Blue );
	  
	  else
	    P[vp].plot( PlotEvents[s], origin, Offset, tfac,
			0.3, Plot::Graph, 1,
			Plot::TriangleUp, 6, Plot::Pixel,
			Plot::Red, Plot::Red );
	  
	  sn++;
	}
      }
    }
    // plot voltage trace:
    PlotElements[c] = P[vp].plot( PlotTraces[vp], origin, Offset, tfac,
				  Plot::Green, 1, Plot::Solid,
				  Plot::Circle, 0, Plot::Green, Plot::Green );
  }
  updateStyle();

  // set xlabel:
  if ( VP.size() > 0 )
    P[VP.back()].setXLabel( tunit );

  // trigger source:
  TriggerSource = -1;
  for ( int s=0; s<PlotEvents.size(); s++ ) {
    if ( (PlotEvents[s].mode() & PlotTriggerMode) ) {
      TriggerSource = s;
      break;
    }
  }


  P.unlock();
	
}


void PlotTrace::plot( void )
{
  P.lock();
  bool plotting = Plotting;
  P.unlock();
  if ( !plotting ) 
    return;

  // get data:
  getData();

  if ( PlotChanged ) {
    init();
    PlotChanged = false;
  }

  P.lock();

  // set left- and rightmargin:
  double tfac = 1000.0;
  double leftwin = 0.0;
  double rightwin = tfac * TimeWindow;
  double sigtime = signalTime();
  if ( sigtime < 0.0 )
    sigtime = 0.0;
  if ( ViewMode == SignalView ) {
    // offset fixed at signalTime():
    leftwin = -tfac * TimeOffs;
    rightwin = leftwin + tfac * TimeWindow;
    LeftTime = leftwin/tfac + sigtime;
    Offset = sigtime;
  }
  else if ( ViewMode == EndView || ( ViewMode == WrapView && TimeWindow < 0.5 ) ) {
    // offset continuous at currentTime():
    rightwin = tfac * ( currentTime() - sigtime );
    leftwin = rightwin - tfac * TimeWindow;
    LeftTime = leftwin/tfac + sigtime;
    Offset = sigtime;
  }
  else if ( ViewMode == WrapView ) {
    // offset wrapped at signalTime():
    leftwin = tfac * ::floor( ( currentTime() - sigtime ) / TimeWindow ) * TimeWindow;
    LeftTime = leftwin/tfac + sigtime;
    rightwin += leftwin;
    Offset = sigtime;
  }
  else {
    // offset fixed at LeftTime:
    leftwin = tfac * ( LeftTime - Offset );
    rightwin = leftwin + tfac * TimeWindow;
  }

  // align to trigger:
  if ( ( ViewMode == EndView || ViewMode == WrapView ) &&
       Trigger && TriggerSource >= 0 ) {
    int ninx = PlotEvents[TriggerSource].next( LeftTime );
    int pinx = PlotEvents[TriggerSource].previous( LeftTime );      
    if ( ninx >= PlotEvents[TriggerSource].size() ) {
      ninx = pinx;
      pinx--;
    }
    if ( ninx > 0 && pinx >= 0 ) {
      // get average period:
      pinx = ninx - 10;
      if ( pinx < 0 ) 
	pinx = 0;
      double dt = (PlotEvents[TriggerSource][ninx] - PlotEvents[TriggerSource][pinx])/(ninx-pinx);
      if ( dt/TimeWindow > 0.02 )
	ninx--;
      double nt = PlotEvents[TriggerSource][ninx];
      if ( TimeWindow < dt || fabs( nt - LeftTime ) <= 2.0*dt ) {
	LeftTime = nt;
	leftwin = (LeftTime - sigtime)*tfac;
	rightwin = leftwin + tfac * TimeWindow;
      }
    }
  }

  // set xrange:
  for ( unsigned int c=0; c<VP.size(); c++ ) {
    // setting axis:
    P[VP[c]].setXRange( leftwin, rightwin );
    if ( ! P[VP[c]].zoomedYRange() )
      P[VP[c]].setYRange( PlotTraces[VP[c]].minValue(), PlotTraces[VP[c]].maxValue() );
  }

  updateStyle();

  // plot:
  P.draw();

  P.unlock();
}


void PlotTrace::updateStyle( void )
{
  // line and pointstyle:
  for ( unsigned int c=0; c<VP.size(); c++ ) {
    if ( PlotElements[c] >= 0 ) {
      double width = P[VP[c]].pixelPlotWidth();
      if ( width < 10.0 )
	width = P[VP[c]].pixelScreenWidth();
      if ( width < 10.0 )
	width = 10.0;
      if ( PlotTraces[VP[c]].indices( TimeWindow )/width > 0.2 )
	P[VP[c]][PlotElements[c]].setPoint( Plot::Circle, 0,
					    Plot::Green, Plot::Green );
      else if ( PlotTraces[VP[c]].indices( TimeWindow )/width > 0.05 )
	P[VP[c]][PlotElements[c]].setPoint( Plot::Circle, 4,
					    Plot::Green, Plot::Green );
      else
	P[VP[c]][PlotElements[c]].setPoint( Plot::Circle, 8,
					    Plot::Green, Plot::Green );
      if ( PlotTraces[VP[c]].indices( TimeWindow )/width > 2.0 )
	P[VP[c]][PlotElements[c]].setLine( Plot::Green, 1, Plot::Solid );
      else
	P[VP[c]][PlotElements[c]].setLine( Plot::Green, 2, Plot::Solid );
    }
  }
}


void PlotTrace::updateRanges( int id )
{
  // P is already locked!
  double tfac = 0.001;
  TimeWindow = tfac * ( P[id].xmaxRange() - P[id].xminRange() );
  TimeOffs = -tfac * P[id].xminRange();
  LeftTime = Offset - TimeOffs;
}


void PlotTrace::addMenu( QMenu *menu )
{
  Menu = menu;

  Menu->addAction( "Zoom &in", this, SLOT( zoomIn() ), Qt::Key_Plus );
  Menu->addAction( "Zoom &out", this, SLOT( zoomOut() ), Qt::Key_Minus );
  Menu->addAction( "Move &left", this, SLOT( moveLeft() ), Qt::Key_PageUp );
  Menu->addAction( "Move &right", this, SLOT( moveRight() ), Qt::Key_PageDown );
  Menu->addAction( "&Begin", this, SLOT( moveStart() ), Qt::CTRL + Qt::Key_PageUp );
  Menu->addAction( "&End", this, SLOT( moveEnd() ), Qt::CTRL + Qt::Key_PageDown );
  Menu->addAction( "Move to signal", this, SLOT( moveToSignal() ), Qt::CTRL + Qt::Key_Home );
  Menu->addAction( "&Signal view", this, SLOT( viewSignal() ), Qt::Key_Home );
  Menu->addAction( "Move offset left", this, SLOT( moveSignalOffsLeft() ), Qt::SHIFT + Qt::Key_PageUp );
  Menu->addAction( "Move offset right", this, SLOT( moveSignalOffsRight() ), Qt::SHIFT + Qt::Key_PageDown );
  Menu->addAction( "&End view", this, SLOT( viewEnd() ), Qt::Key_End );
  Menu->addAction( "&Wrapped view", this, SLOT( viewWrapped() ), Qt::Key_Insert );
  Menu->addAction( "&Trigger", this, SLOT( toggleTrigger() ), Qt::CTRL + Qt::Key_T );
  Menu->addAction( "&Manual", this, SLOT( manualRange() ), Qt::CTRL + Qt::Key_M );
  Menu->addAction( "&Auto", this, SLOT( autoRange() ), Qt::CTRL + Qt::Key_A );
  Menu->addAction( "Center &vertically", this, SLOT( centerVertically() ) );
  Menu->addAction( "&Zoom vertically", this, SLOT( centerZoomVertically() ) );
  Menu->addAction( "&Toggle Plot", this, SLOT( plotOnOff() ) );

  Menu->addSeparator();

  PlotProps.clear();
  for ( int k=0; k<traces().size(); k++ ) {
    string s = "&" + Str( k+1 );
    s += " ";
    s += trace(k).ident();
    PlotProps.push_back( PlotProperties() );
    PlotProps.back().Action = Menu->addAction( s.c_str() );
    PlotProps.back().Action->setCheckable( true );
    PlotProps.back().Visible = true;
    PlotProps.back().Action->setChecked( PlotProps.back().Visible );
  }
  connect( Menu, SIGNAL( triggered( QAction* ) ),
	   this, SLOT( toggle( QAction* ) ) );
}


void PlotTrace::updateMenu( void )
{
  if ( Menu != 0 ) {

    P.lock();
  
    // remove old traces:
    for ( unsigned int k=0; k<PlotProps.size(); k++ )
      Menu->removeAction( PlotProps[k].Action );

    // get traces and events:
    PlotTraces.clear();
    PlotEvents.clear();
    if ( FilePlot ) {
      PlotTraces.add( FileTraces );
      PlotEvents.add( FileEvents );
    }
    else {
      PlotTraces.add( traces() );
      PlotEvents.add( events() );
    }

    // add new traces:
    PlotProps.resize( PlotTraces.size() );
    for ( int k=0; k<PlotTraces.size(); k++ ) {
      string s = "&" + Str( k+1 );
      s += " ";
      s += PlotTraces[k].ident();
      PlotProps[k].Action = Menu->addAction( s.c_str() );
      PlotProps[k].Action->setCheckable( true );
      PlotProps[k].Action->setChecked( PlotProps.back().Visible );
    }

    P.unlock();
  }
}


void PlotTrace::start( double time )
{
  PlotTimer.start( (int)::rint( 1000.0*time ) );
}


void PlotTrace::stop( void )
{
  PlotTimer.stop();
}


double PlotTrace::signalTime( void ) const
{
  if ( FilePlot && FileTraces.size() > 0 )
    return FileTraces[0].signalTime();
  else
    return RELACSPlugin::signalTime();
}


double PlotTrace::currentTime( void ) const
{
  if ( FilePlot && FileTraces.size() > 0 )
    return FileTraces[0].currentTime();
  else
    return RELACSPlugin::currentTime();
}


void PlotTrace::setPlotOn( bool on )
{
  P.lock();
  AutoOn = on;
  bool man = Manual;
  P.unlock();

  if ( man )
    return;

  P.lock();

  // toggle plot:
  Plotting = on;
  postCustomEvent( 11 );  // animate on/off button
  P.unlock();
}


void PlotTrace::setPlotOff( void )
{
  setPlotOn( false );
}


void PlotTrace::setPlotSignal( double length, double offs )
{
  P.lock();
  AutoOn = true;
  AutoFixed = true;
  AutoTime = length;
  AutoOffs = offs;
  bool man = Manual;
  P.unlock();

  if ( man )
    return;

  P.lock();

  // toggle plot:
  Plotting = true;
  postCustomEvent( 11 );  // animate on/off button

  // toggle fixed offset:
  setView( SignalView );

  // length of total time window:
  TimeWindow = length;

  // left offset to signal:
  TimeOffs = offs;

  updateStyle();

  P.unlock();
}


void PlotTrace::setPlotSignal( void )
{
  P.lock();
  AutoOn = true;
  AutoFixed = true;
  bool man = Manual;
  P.unlock();

  if ( man )
    return;

  P.lock();

  // toggle plot:
  Plotting = true;
  postCustomEvent( 11 );  // animate on/off button

  // toggle fixed offset:
  setView( SignalView );

  P.unlock();
}


void PlotTrace::setPlotContinuous( double length )
{
  P.lock();
  AutoOn = true;
  AutoFixed = false;
  AutoTime = length;
  bool man = Manual;
  P.unlock();

  if ( man )
    return;

  P.lock();

  // toggle plot:
  Plotting = true;
  postCustomEvent( 11 );  // animate on/off button

  // toggle fixed offset:
  setView( ContView );

  // length of total time window:
  TimeWindow = length;

  updateStyle();

  P.unlock();
}


void PlotTrace::setPlotContinuous( void )
{
  P.lock();
  AutoOn = true;
  AutoFixed = false;
  bool man = Manual;
  P.unlock();

  if ( man )
    return;

  P.lock();

  // toggle plot:
  Plotting = true;
  postCustomEvent( 11 );  // animate on/off button

  // toggle fixed offset:
  setView( ContView );

  P.unlock();
}


void PlotTrace::zoomOut( void )
{
  P.lock();
  TimeWindow *= 2;
  TimeOffs *= 2;
  P.unlock();
  if ( RW->idle() )
    plot();
  else {
    P.lock();
    updateStyle();
    P.unlock();
  }
}


void PlotTrace::zoomIn( void )
{
  P.lock();
  TimeWindow /= 2.0;
  TimeOffs /= 2.0;
  P.unlock();
  if ( RW->idle() )
    plot();
  else {
    P.lock();
    updateStyle();
    P.unlock();
  }
}


void PlotTrace::moveLeft( void )
{
  P.lock();
  if ( ViewMode != FixedView )
    setView( FixedView );
  else
    LeftTime -= 0.5 * TimeWindow;
  P.unlock();
  if ( RW->idle() )
    plot();
}


void PlotTrace::moveRight( void )
{
  P.lock();
  if ( ViewMode != FixedView )
    setView( FixedView );
  else
    LeftTime += 0.5 * TimeWindow;
  P.unlock();
  if ( RW->idle() )
    plot();
}


void PlotTrace::moveStart( void )
{
  P.lock();
  if ( ViewMode != FixedView )
    setView( FixedView );
  LeftTime = 0.0;
  P.unlock();
  if ( RW->idle() )
    plot();
}


void PlotTrace::moveEnd( void )
{
  P.lock();
  if ( ViewMode != FixedView )
    setView( FixedView );
  LeftTime = currentTime() - TimeWindow;
  P.unlock();
  if ( RW->idle() )
    plot();
}


void PlotTrace::moveToSignal( void )
{
  P.lock();
  if ( ViewMode == SignalView )
    TimeOffs = 0.0;
  else {
    if ( ViewMode != FixedView )
      setView( FixedView );
    double sigtime = signalTime();
    if ( sigtime < 0.0 )
      sigtime = 0.0;
    LeftTime = sigtime;
  }
  P.unlock();
  if ( RW->idle() )
    plot();
}


void PlotTrace::viewSignal( void )
{
  P.lock();
  setView( SignalView );
  P.unlock();
  if ( RW->idle() )
    plot();
}


void PlotTrace::moveSignalOffsLeft( void )
{
  P.lock();
  if ( ViewMode == SignalView || ViewMode == FixedView ) {
    TimeOffs += 0.25 * TimeWindow;
    if ( ViewMode == FixedView )
      LeftTime = Offset - TimeOffs;
  }
  else
    setView( SignalView );
  P.unlock();
  if ( RW->idle() )
    plot();
}


void PlotTrace::moveSignalOffsRight( void )
{
  P.lock();
  if ( ViewMode == SignalView || ViewMode == FixedView ) {
    TimeOffs -= 0.25 * TimeWindow;
    if ( ViewMode == FixedView )
      LeftTime = Offset - TimeOffs;
  }
  else
    setView( SignalView );
  P.unlock();
  if ( RW->idle() )
    plot();
}


void PlotTrace::viewEnd( void )
{
  P.lock();
  setView( EndView );
  P.unlock();
  if ( RW->idle() )
    plot();
}


void PlotTrace::viewWrapped( void )
{
  P.lock();
  setView( WrapView );
  P.unlock();
  if ( RW->idle() )
    plot();
}


void PlotTrace::toggleTrigger( void )
{
  P.lock();
  Trigger = ! Trigger;
  P.unlock();
  if ( RW->idle() )
    plot();
}


void PlotTrace::plotOnOff( )
{
  P.lock();
  Plotting = ! Plotting;
  int p = Plotting;
  P.unlock();
  if ( OnOffButton != 0 )
    OnOffButton->setDown( ! p );
}


void PlotTrace::toggleManual( void )
{
  P.lock();
  bool man = Manual;
  P.unlock();
  if ( man ) 
    autoRange();
  else
    manualRange();
}


void PlotTrace::manualRange( void )
{
  P.lock();
  Manual = true;
  if ( ManualButton != 0 )
    ManualButton->setDown( Manual );
  P.unlock();
}


void PlotTrace::autoRange( void )
{
  P.lock();
  Manual = false;
  if ( ManualButton != 0 )
    ManualButton->setDown( Manual );

  // toggle plot:
  Plotting = AutoOn;
  if ( OnOffButton != 0 ) {
    if ( Plotting )
      OnOffButton->setDown( false );
    else
      OnOffButton->setDown( true );
  }

  // toggle fixed offset:
  setView( AutoFixed ? SignalView : WrapView );

  // length of total time window:
  TimeWindow = AutoTime;

  // left offset to signal:
  TimeOffs = AutoOffs;

  updateStyle();

  P.unlock();
}


void PlotTrace::centerVertically( void )
{
  // select plots to be centered:
  vector<int> cp;
  cp.reserve( VP.size() );
  cp.clear();
  P.lock();
  if ( VP.size() == 1 )
    cp.push_back( VP[0] );
  else {
    for ( unsigned int c=0; c<VP.size(); c++ ) {
      if ( PlotTraces[c].mode() & PlotTraceCenterVertically )
	cp.push_back( VP[c] );
    }
  }
  P.unlock();
  if ( cp.empty() )
    return;

  // center plots:
  P.lock();
  double tfac = 1000.0;
  for ( unsigned int c=0; c<cp.size(); c++ ) {
    double xmin = P[cp[c]].xminRange()/tfac + Offset;
    double xmax = P[cp[c]].xmaxRange()/tfac + Offset;
    // make sure, there are enough data shown in the window:
    if ( PlotTraces[cp[c]].currentTime() < 0.5*(xmin+xmax) ) {
      xmin = PlotTraces[cp[c]].currentTime() - (xmax-xmin);
      xmax = PlotTraces[cp[c]].currentTime();
    }
    float min = 0.0;
    float max = 0.0;
    PlotTraces[cp[c]].minMax( min, max, xmin, xmax );
    if ( P[cp[c]].ranges() == 0 )
      P[cp[c]].pushRanges();
    double center = 0.5*(min+max);
    double range = 0.5*(P[cp[c]].ymaxRange() - P[cp[c]].yminRange());
    double nmin = center-range;
    double nmax = center+range;
    if ( nmin < PlotTraces[cp[c]].minValue() ) {
      nmin = PlotTraces[cp[c]].minValue();
      nmax = nmin + 2.0*range;
    }
    if ( nmax > PlotTraces[cp[c]].maxValue() ) {
      nmax = PlotTraces[cp[c]].maxValue();
      nmin = nmax - 2.0*range;
    }
    P[cp[c]].setYRange( nmin, nmax );
  }
  P.unlock();
}


void PlotTrace::centerZoomVertically( void )
{
  // select plots to be centered:
  vector<int> cp;
  cp.reserve( VP.size() );
  cp.clear();
  P.lock();
  if ( VP.size() == 1 )
    cp.push_back( VP[0] );
  else {
    for ( unsigned int c=0; c<VP.size(); c++ ) {
      if ( PlotTraces[c].mode() & PlotTraceCenterVertically )
	cp.push_back( VP[c] );
    }
  }
  P.unlock();
  if ( cp.empty() )
    return;

  // center plots:
  P.lock();
  double tfac = 1000.0;
  for ( unsigned int c=0; c<cp.size(); c++ ) {
    double xmin = P[cp[c]].xminRange()/tfac + Offset;
    double xmax = P[cp[c]].xmaxRange()/tfac + Offset;
    // make sure, there are enough data shown in the window:
    if ( PlotTraces[cp[c]].currentTime() < 0.5*(xmin+xmax) ) {
      xmin = PlotTraces[cp[c]].currentTime() - (xmax-xmin);
      xmax = PlotTraces[cp[c]].currentTime();
    }
    float min = 0.0;
    float max = 0.0;
    PlotTraces[cp[c]].minMax( min, max, xmin, xmax );
    if ( P[cp[c]].ranges() == 0 )
      P[cp[c]].pushRanges();
    double center = 0.5*(min+max);
    double range = 0.6*(max-min);
    double nmin = center-range;
    double nmax = center+range;
    if ( nmin < PlotTraces[cp[c]].minValue() ) {
      nmin = PlotTraces[cp[c]].minValue();
      nmax = nmin + 2.0*range;
    }
    if ( nmax > PlotTraces[cp[c]].maxValue() ) {
      nmax = PlotTraces[cp[c]].maxValue();
      nmin = nmax - 2.0*range;
    }
    P[cp[c]].setYRange( nmin, nmax );
  }
  P.unlock();
}


void PlotTrace::viewToggle( void )
{
  P.lock();
  if ( ViewMode == EndView )
    setView( SignalView );
  else if ( ViewMode == SignalView )
    setView( WrapView );
  else if ( ViewMode == WrapView )
    setView( EndView );
  P.unlock();
}


void PlotTrace::setView( Views mode )
{
  if ( ViewMode != mode ) {
    if ( mode == EndView || mode == WrapView )
      ContinuousView = mode;
    if ( mode == ContView )
      mode = ContinuousView;
    ViewMode = mode;
    PlotChanged = true;
    postCustomEvent( 12 );
  }
}


void PlotTrace::displayIndex( const string &fpath, const deque<int> &traceindex,
			      const deque<int> &eventsindex, double time )
{
  bool success = true;
  bool newfile = false;
  P.lock();
  // XXX this comparison needs the absolute path of both sides:
  if ( Str( fpath ).dir() == defaultPath() ) {
    if ( FilePlot ) {
      newfile = true;
      FilePlot = false;
      FileHeader.clear();
      FileTraces.clear();
      FileTracesNames.clear();
      FileSizes.clear();
      FileEvents.clear();
      FileEventsNames.clear();
    }
  }
  else {
    if ( FilePath != fpath ) {
      newfile = true;
      FileHeader.clear();
      FileTraces.clear();
      FileTracesNames.clear();
      FileSizes.clear();
      FileEvents.clear();
      FileEventsNames.clear();

      // read in meta data:
      DataFile sf( fpath );
      sf.readMetaData();
      Options header = sf.metaDataOptions( sf.levels()-1 );
      if ( header.empty() )
	success = false;
      else
	FileHeader = header;

      // create traces:
      for ( int k=1; ; k++ ) {
	string tracename = FileHeader.text( "identifier" + Str( k ), "" );
	if ( tracename.empty() )
	  break;
	string tracefile = FileHeader.text( "data file" + Str( k ) );
	double tracestep = FileHeader.number( "sample interval" + Str( k ), "ms" );
	// we need to guess the exact rate:
	double rate = ::round( 1.0/tracestep );
	tracestep = 0.001/rate;
	string traceunit = FileHeader.text( "unit" + Str( k ) );

	InData id;
	id.setIdent( tracename );
	id.setTrace( k-1 );
	id.setStepsize( tracestep );
	id.setStartSource( 0 );
	id.setUnipolar( false );
	id.setUnit( 1.0, traceunit );
	id.setChannel( k-1 );
	id.setMinValue( -150.0 );  // XXX get this somehow from the configuration file?
	id.setMaxValue( +150.0 );
	id.setDevice( 0 );
	id.setContinuous();
	int m = PlotTraceMode;
	/*
	  if ( boolean( "inputtracecenter", k, false ) )
	  m |= PlotTraceCenterVertically;
	*/
	id.setMode( m );
	// id.setReference( text( "inputtracereference", k, InData::referenceStr( InData::RefGround ) ) );
	id.setGainIndex( 0 );
	FileTraces.push( id );
	// FileTraces.back().reserve( id.indices( number( "inputtracecapacity", 0, 1000.0 ) ) );
	FileTraces.back().reserve( 60.0 );
	FileTracesNames.push_back( Str( fpath ).dir() + tracefile );
	ifstream is( FileTracesNames.back().c_str() );
	if ( k-1 >= (int)FileSizes.size() ) {
	  is.seekg( 0, is.end );
	  FileSizes.push_back( is.tellg()/sizeof(float) );
	}
      }
      for ( int k=1; ; k++ ) {
	string eventfile = FileHeader.text( "event file" + Str( k ), "" );
	if ( eventfile.empty() )
	  break;
	FileEventsNames.push_back( Str( fpath ).dir() + eventfile );
	DataFile sf( FileEventsNames.back() );
	sf.readMetaData();
	Options header = sf.metaDataOptions( sf.levels()-1 );
	bool eventsizes = ( sf.key().columns() > 1 );
	bool eventwidths = ( sf.key().columns() > 2 );
	EventData ed( 1000, eventsizes, eventwidths );
	ed.setCyclic();
	string eventname = header.text( "events" );
	ed.setIdent( eventname );
	if ( eventname == "Stimulus" )
	  ed.setMode( StimulusEventMode );
	else if ( eventname == "Restart" )
	  ed.setMode( RestartEventMode );
	else if ( eventname == "Recording" )
	  ed.setMode( RecordingEventMode );
	ed.setMode( ed.mode() | PlotTraceMode );
	if ( eventsizes ) {
	  ed.setSizeName( sf.key().name( 1 ) );
	  ed.setSizeUnit( sf.key().unit( 1 ) );
	}
	if ( eventwidths ) {
	  ed.setWidthName( sf.key().name( 2 ) );
	  ed.setWidthUnit( sf.key().unit( 2 ) );
	}
	FileEvents.push( ed );
      }
    }
    if ( FilePlot && FileTraces.size() > 0 && FileTraces[0].size() > 0 &&
	 FileTraces[0].minIndex() <= traceindex[0] && 
	 ( FileTraces[0].size() == FileSizes[0] ||
	   FileTraces[0].size() >= traceindex[0] + FileTraces[0].indices( 10.0 ) ) ) {
      time = FileTraces[0].pos( traceindex[0] );
    }
    else {
      // load data:
      for ( int k=0; k<FileTraces.size(); k++ ) {
	int index = traceindex[k];
	index -= FileTraces.back().indices( 10.0 );  // we want to have some time before as well, at least TimeOffs!
	if ( index + FileTraces.back().capacity() > FileSizes[k] )
	  index = FileSizes[k] - FileTraces.back().capacity();
	if ( index < 0 )
	  index = 0;
	ifstream is( FileTracesNames[k].c_str() );
	FileTraces[k].loadBinary( is, index );
	FileTraces[k].setSignalIndex( traceindex[k] );
	time = FileTraces[k].pos( traceindex[k] );
      }
      for ( int k=0; k<FileEvents.size(); k++ ) {
	DataFile sf( FileEventsNames[k] );
	sf.read( 10 );
	int index = eventsindex[k];
	if ( index >=0 && index < sf.data().rows() )
	  for ( ; index>=0 && sf.data( 0, index ) > time-10.0; index-- );
	if ( index + FileEvents[k].capacity() > sf.data().rows() )
	  index = sf.data().rows() - FileEvents[k].capacity();
	if ( index < 0 )
	  index = 0;
	FileEvents[k].set( index, sf.col( 0 ), sf.col( 1 ), sf.col( 2 ) );  
      }
    }
    FilePlot = true;
  }
  if ( success ) {
    if ( ViewMode != FixedView ) {
      TimeOffs = 0.0;
      setView( FixedView );
    }
    Offset = time;
    LeftTime = Offset - TimeOffs;
    PlotChanged = true;
  }
  P.unlock();
  if ( success && newfile ) {
    updateMenu();
    resize();
  }
  FilePath = fpath;
  if ( RW->idle() )
    plot();
}


void PlotTrace::displayData( void )
{
  P.lock();
  if ( FilePlot ) {
    FilePlot = false;
    FilePath = path();
    FileHeader.clear();
    FileTraces.clear();
    FileTracesNames.clear();
    FileSizes.clear();
    FileEvents.clear();
    FileEventsNames.clear();
    setView( WrapView );
    Offset = trace( 0 ).currentTime() - TimeWindow;
    LeftTime = Offset - TimeOffs;
    PlotChanged = true;
    P.unlock();
    updateMenu();
    resize();
    if ( RW->idle() )
      plot();
  }
}


void PlotTrace::keyPressEvent( QKeyEvent *event )
{
  event->accept();
  switch ( event->key() ) {

  case Qt::Key_1: case Qt::Key_2: case Qt::Key_3: case Qt::Key_4: case Qt::Key_5: 
  case Qt::Key_6: case Qt::Key_7: case Qt::Key_8: case Qt::Key_9: {
    int n = event->key() - Qt::Key_1;
    if ( n >=0 && n < (int)PlotProps.size() )
      toggle( PlotProps[n].Action );
    break;
  }

  case Qt::Key_Equal:
    zoomIn();
    break;

  case Qt::Key_V:
    if ( event->modifiers() & Qt::ShiftModifier )
      centerZoomVertically();
    else
      centerVertically();
    break;

  default: {
    RELACSPlugin::keyPressEvent( event );
  }

  }
}


void PlotTrace::resizeLayout( void )
{
  if ( VP.empty() )
    return;

  if ( VP.size() == 1 ) {
    P[VP.front()].setOrigin( 0.0, 0.0 );
    P[VP.front()].setSize( 1.0, 1.0 );
    return;
  }

  int columns = 1;
  if ( VP.size() > 6 )
    columns = 2;
  int rows = VP.size()/columns;
  if ( VP.size()%columns > 0 )
    rows++;
  double xsize = 1.0/columns;
  double yboffs = double( P[VP.front()].fontPixel( 2.3 ) ) / double( P.height() );
  double ytoffs = double( P[VP.front()].fontPixel( 0.8 ) ) / double( P.height() );
  double yheight = (1.0-yboffs-ytoffs)/rows;
    
  int c = 0;
  int r = 0;
  for ( unsigned int k=0; k<VP.size(); k++ ) {
    P[VP[k]].setOrigin( c*xsize, yboffs+(rows-r-1)*yheight );
    P[VP[k]].setSize( xsize, yheight );
    r++;
    if ( r >= rows ) {
      c++;
      r = 0;
    }
  }
  P[VP.front()].setSize( xsize, yheight+ytoffs );
  P[VP.back()].setOrigin( (columns-1)*xsize, 0.0 );
  P[VP.back()].setSize( xsize, yheight+yboffs );
  if ( columns > 1 ) {
    P[VP[(VP.size()+1)/2]].setSize( xsize, yheight+ytoffs );
    P[VP[(VP.size()-1)/2]].setOrigin( 0.0, 0.0 );
    P[VP[(VP.size()-1)/2]].setSize( xsize, yheight+yboffs );
  }
}


void PlotTrace::resizeEvent( QResizeEvent *qre )
{
  /*
  P.lock();
  resizeLayout();
  P.unlock();

  P.resizeEvent( qre );
  */
  if ( ButtonBox == 0 )
    return;

  QSize bs = ButtonBox->sizeHint();
  /*
  ButtonBox->setGeometry( width() - bs.width(), 0,
			  bs.width(), bs.height() );
  */
  ButtonBox->setGeometry( 0, 0, bs.width(), bs.height() );
}


void PlotTrace::resizePlots( QResizeEvent *qre )
{
  P.lock();
  resizeLayout();
  P.unlock();
}


void PlotTrace::customEvent( QEvent *qce )
{
  switch ( qce->type() - QEvent::User ) {

  case 11 :
    if ( OnOffButton != 0 ) {
      P.lock();
      bool plotting = Plotting;
      P.unlock();
      if ( plotting )
	OnOffButton->setDown( false );
      else
	OnOffButton->setDown( true );
    }
    break;

  case 12:
    if ( ViewButton != 0 ) {
      if ( ViewMode == SignalView )
	ViewButton->setIcon( SignalViewIcon );
      else if ( ViewMode == WrapView )
	printlog( "warnig: PlotTrace WrapIcon not implemented yet!" );
      else
	ViewButton->setIcon( EndViewIcon );
    }
    break;

  }
}



}; /* namespace relacs */

#include "moc_plottrace.cc"

