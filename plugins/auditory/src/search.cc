/*
  auditory/search.cc
  Periodically emits a search stimulus.

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

#include <QWidget>
#include <QWidget>
#include <QLabel>
#include <QGroupBox>
#include <QButtonGroup>
#include <relacs/tablekey.h>
#include <relacs/auditory/search.h>
using namespace relacs;

namespace auditory {


const double Search::ShortIntensityStep = 1.0;
const double Search::LongIntensityStep = 5.0;
const double Search::MaxIntensity = 100.0;
const double Search::MinIntensity = 0.0;  

const double Search::ShortDurationStep = 0.005;
const double Search::LongDurationStep = 0.05;
const double Search::MaxDuration = 10.0;
const double Search::MinDuration = 0.005;

const double Search::ShortPauseStep = 0.005;
const double Search::LongPauseStep = 0.05;
const double Search::MaxPause = 10.0;
const double Search::MinPause = 0.005;

const double Search::ShortFrequencyStep = 500.0;
const double Search::LongFrequencyStep = 5000.0;
const double Search::MaxFrequency = 40000.0;
const double Search::MinFrequency = 2000.0;


Search::Search( void )
  : RePro( "Search", "Auditory", "Jan Benda and Christian Machens", "2.2", "Jan 10, 2008" )
{
  // parameter:
  Intensity = 80.0;
  Duration = 0.05;
  Pause = 0.5;
  PrePause = 0.05;
  Frequency = 5000.0;
  Waveform = 0;
  SearchLeft = false;
  SetBestSide = 1;

  // options:
  addNumber( "intensity", "Intensity",  Intensity, MinIntensity, MaxIntensity, ShortIntensityStep, "dB", "dB", "%.1f" ).setActivation( "mute", "false" );
  addBoolean( "mute", "No stimulus", false );
  addNumber( "duration", "Duration of stimulus", Duration, MinDuration, MaxDuration, ShortDurationStep, "sec", "ms" );
  addNumber( "pause", "Duration of pause", Pause, MinPause, MaxPause, ShortPauseStep, "sec", "ms" );
  addNumber( "prepause", "Part of pause before stimulus", PrePause, 0.0, MaxPause, ShortPauseStep, "sec", "ms" );
  addNumber( "frequency", "Frequency of stimulus", Frequency, MinFrequency, MaxFrequency, ShortFrequencyStep, "Hz", "kHz" );
  addSelection( "waveform", "Waveform of stimulus", "sine|noise" );
  addNumber( "ramp", "Ramp", 0.002, 0.0, 10.0, 0.001, "sec", "ms" );
  addSelection( "side", "Speaker", "left|right|best" );
  addInteger( "repeats", "Number of repetitions", 0, 0, 10000, 2 );
  addBoolean( "adjust", "Adjust input gains", true );
  addBoolean( "saving", "Save raw data", true );
  addSelection( "setbestside", "Set the sessions's best side", "never|no session|always" );
  addBoolean( "keep", "Keep changes", true );

  // variables:
  NewSignal = true;
  Mute = false;

  // layout:
  QGridLayout *grid = new QGridLayout;
  setLayout( grid );

  // Intensity Settings:
  ILCD = new LCDRange( "Intensity (dB)", 3 );
  ILCD->setRange( int(MinIntensity), int(MaxIntensity) );
  ILCD->setValue( int(Intensity) );
  ILCD->setSteps( int(ShortIntensityStep), int(LongIntensityStep) );
  grid->addWidget( ILCD, 0, 0 );
  QObject::connect( ILCD, SIGNAL( valueChanged( int ) ), 
		    (QWidget*)this, SLOT( setIntensity( int ) ) );

  QGridLayout *sgrid = new QGridLayout;
  grid->addLayout( sgrid, 0, 1 );

  // Duration Settings:
  DLCD = new LCDRange( "Stimulus (msec)", 4 );
  DLCD->setRange( int(1000.0*MinDuration), int(1000.0*MaxDuration) );
  DLCD->setValue( int(1000.0*Duration) );
  DLCD->setSteps( int(1000.0*ShortDurationStep), int(1000.0*LongDurationStep) );
  sgrid->addWidget( DLCD, 0, 0 );
  QObject::connect( DLCD, SIGNAL( valueChanged( int ) ), 
		    (QWidget*)this, SLOT( setDuration( int ) ) );

  // Pause Settings:
  PLCD = new LCDRange( "Pause (msec)", 4 );
  PLCD->setRange( int(1000.0*MinPause), int(1000.0*MaxPause) );
  PLCD->setValue( int(1000.0*Pause) );
  PLCD->setSteps( int(1000.0*ShortPauseStep), int(1000.0*LongPauseStep) );
  sgrid->addWidget( PLCD, 0, 1 );
  QObject::connect( PLCD, SIGNAL( valueChanged( int ) ), 
		    (QWidget*)this, SLOT( setPause( int ) ) );

  // Waveform:
  QGroupBox *gb = new QGroupBox( "Waveform" );
  sgrid->addWidget( gb, 1, 0 );
  SineButton = new QRadioButton( "Sine" );
  NoiseButton = new QRadioButton( "Noise" );
  SineButton->setChecked( true );
  QVBoxLayout *gbl = new QVBoxLayout;
  gbl->addWidget( SineButton );
  gbl->addWidget( NoiseButton );
  gb->setLayout( gbl );
  QButtonGroup *WaveformButtons = new QButtonGroup;
  WaveformButtons->addButton( SineButton, 0 );
  WaveformButtons->addButton( NoiseButton, 1 );
  QObject::connect( WaveformButtons, SIGNAL( buttonClicked( int ) ), 
		    (QWidget*)this, SLOT( setWaveform( int ) ) );

  // Frequency Settings:
  FLCD = new LCDRange( "Frequency (Hz)", 5 );
  FLCD->setRange( int(MinFrequency), int(MaxFrequency) );
  FLCD->setValue( int(Frequency) );
  FLCD->setSteps( int(ShortFrequencyStep), int(LongFrequencyStep) );
  sgrid->addWidget( FLCD, 1, 1 );
  QObject::connect( FLCD, SIGNAL( valueChanged( int ) ), 
		    (QWidget*)this, SLOT( setFrequency( int ) ) );

  // mute button:
  MuteButton = new QPushButton;
  MuteButton->setCheckable( true );
  MuteButton->setText( "Mute" );
  grid->addWidget( MuteButton, 1, 0 );
  QObject::connect( MuteButton, SIGNAL( clicked() ), 
		    (QWidget*)this, SLOT( toggleMute() ) ); 

  // SearchSide Settings:
  QHBoxLayout *hbox = new QHBoxLayout;
  grid->addLayout( hbox, 1, 1 );
  QLabel *label = new QLabel( "Speaker:" );
  hbox->addWidget( label );
  LeftButton = new QRadioButton( "left" );
  hbox->addWidget( LeftButton );
  RightButton = new QRadioButton( "right" );
  hbox->addWidget( RightButton );
  QObject::connect( LeftButton, SIGNAL( clicked() ), 
		    (QWidget*)this, SLOT( setSpeakerLeft() ) ); 
  QObject::connect( RightButton, SIGNAL( clicked() ), 
		    (QWidget*)this, SLOT( setSpeakerRight() ) ); 
}


Search::~Search( void )
{
}


int Search::main( void )
{
  // get options:
  if ( ! boolean( "saving" ) )
    noSaving();
  Intensity = int( number( "intensity" ) );
  Mute = boolean( "mute" );
  Duration = number( "duration" );
  Pause = number( "pause" );
  PrePause = number( "prepause" );
  Frequency = int( number( "frequency" ) );
  Waveform = index( "waveform" );
  double ramp = number( "ramp" );
  int side = index( "side" );
  int repeats = integer( "repeats" );
  bool adjustgain = boolean( "adjust" );
  SetBestSide = index( "setbestside" );
  bool keepchanges = boolean( "keep" );

  if ( side > 1 )
    side = metaData( "Cell" ).index( "best side" );
  SearchLeft = ( side == 0 );

  // update widgets:
  postCustomEvent( 11 );

  // don't print repro message:
  if ( repeats <= 0 )
    noMessage();

  if ( SetBestSide + ( sessionRunning() ? 0 : 1 ) > 1 )
    metaData( "Cell" ).selectText( "best side", SearchLeft ? "left" : "right" );

  // Header:
  Options header;
  header.addInteger( "index" );
  header.addText( "session time" );
  header.addLabel( "settings:" );

  // stimulus:
  OutData signal;
  signal.setDelay( PrePause );
  double meanintensity = 0.0;
  NewSignal = true;

  // plot trace:
  plotToggle( true, true, 1.25*Duration, 0.125*Duration );

  timeStamp();

  for ( int count=0;
	( repeats <= 0 || count < repeats ) && softStop() == 0; 
	count++ ) {

    // message:
    if ( repeats == 0 && count%60 == 0 )
      message( "Search ..." );
    else if ( repeats > 0 )
      message( "Search loop <b>" + Str( count ) + "</b> of <b>" + Str( repeats ) + "</b>" );

    // create stimulus:
    if ( NewSignal ) {
      signal.clear();
      signal.setTrace( SearchLeft ? LeftSpeaker[0] : RightSpeaker[0]  );
      applyOutTrace( signal );  // needed for maximum sampling rate!
      if ( Mute ) {
	signal.setSampleInterval( 0.0001 );
	signal.resize( 10 );
	signal = 0.0;
      }
      else {
	if ( Waveform == 1 ) {
	  signal.bandNoiseWave( MinFrequency, Frequency, Duration, 0.3, ramp );
	  ::relacs::clip( -1.0, 1.0, signal );
	  meanintensity = 6.0206; // stdev=0.5
	  meanintensity = 10.458; // stdev = 0.3
	}
	else {
	  signal.sineWave( Frequency, Duration, 1.0, ramp );
	  meanintensity = 3.0103;
	}
      }
      signal.setIntensity( Intensity > 0 && ! Mute ? Intensity + meanintensity : OutData::MuteIntensity );
      //      convert( signal );
      NewSignal = false;
    }
    else {
      signal.setIntensity( Intensity > 0 ? Intensity + meanintensity : OutData::MuteIntensity );
      signal.setTrace( SearchLeft ? LeftSpeaker[0] : RightSpeaker[0]  );
    }

    // output stimulus:
    write( signal );
    if ( signal.error() ) {
      // Attenuator overflow or underflow.
      // Set intensity appropriately and write stimulus again.
      if ( signal.underflow() )
	setIntensity( int( ceil( signal.intensity() - meanintensity )) );
      else
	setIntensity( int( floor( signal.intensity() - meanintensity )) );
      postCustomEvent( 12 );
      write( signal );
    }

    sleepOn( Duration + Pause );
    if ( interrupt() ) {
      writeZero( SearchLeft ? LeftSpeaker[0] : RightSpeaker[0] );
      if ( keepchanges )
	setToDefaults();
      return Aborted;
    }
    timeStamp();

    // adjust gain of daq board:
    if ( adjustgain ) {
      for ( int k=0; k<traces().size(); k++ )
	adjust( trace( k ), trace( k ).signalTime(),
		trace( k ).signalTime() + Duration, 0.8 );
      //      activateGains();
    }

    // save:
    if ( repeats > 0 ) {
      if ( count == 0 ) {
	header.setInteger( "index", totalRuns() );
	header.setText( "session time", sessionTimeStr() );
      }
      for ( int trace=1; trace < events().size(); trace++ ) {
	saveEvents( events( trace ), count, header );
      }
    }

  }

  setMessage();
  writeZero( SearchLeft ? LeftSpeaker[0] : RightSpeaker[0] );
  if ( keepchanges )
    setToDefaults();
  return Completed;
}


void Search::saveEvents( const EventData &events, int count, const Options &header )
{
  // create file:
  ofstream df;
  df.open( addPath( "search-" + Str( events.ident() ).lower() + "-events.dat" ).c_str(),
	   ofstream::out | ofstream::app );
  if ( ! df.good() )
    return;

  // write header and key:
  TableKey spikeskey;
  spikeskey.addNumber( "time", "ms", "%9.2f" );
  if ( count == 0 ) {
    df << '\n' << '\n';
    header.save( df, "# " );
    Options::save( df, "#   ", -1, 0, false, true );
    df << '\n';
    spikeskey.saveKey( df, true, false );
  }

  // write data:
  double t0 = events.signalTime();
  df << '\n';
  df << "# trial: " << count << '\n';
  if ( events.count( t0-PrePause, t0-PrePause+Duration+Pause ) <= 0 ) {
    df << "  -0" << '\n';
  }
  else {
    long jmax = events.previous( t0 + Duration + Pause - PrePause );
    for ( long j=events.next( t0-PrePause ); j<=jmax; j++ ) {
      spikeskey.save( df, 1000.0 * ( events[j] - t0 ), 0 );
      df << '\n';
    }
  }
  
}


void Search::keyPressEvent( QKeyEvent *qke )
{
  switch ( qke->key()) {
  case Qt::Key_Up:
    if ( qke->modifiers() & Qt::AltModifier ) {
      if ( qke->modifiers() & Qt::ShiftModifier )
	setFrequency( int(Frequency + LongFrequencyStep) );
      else
	setFrequency( int(Frequency + ShortFrequencyStep) );
      FLCD->setValue( int(Frequency) );
    }
    else {
      if ( qke->modifiers() & Qt::ShiftModifier )
	setIntensity( int(Intensity + LongIntensityStep) );
      else
	setIntensity( int(Intensity + ShortIntensityStep) );
      ILCD->setValue( int(Intensity) );
    }
    break;
  case Qt::Key_Down:                // arrow down
    if ( qke->modifiers() & Qt::AltModifier ) {
      if ( qke->modifiers() & Qt::ShiftModifier )
	setFrequency( int(Frequency - LongFrequencyStep) );
      else
	setFrequency( int(Frequency - ShortFrequencyStep) );
      FLCD->setValue( int(Frequency) );
    }
    else {
      if ( qke->modifiers() & Qt::ShiftModifier )
	setIntensity( int(Intensity - LongIntensityStep) );
      else
	setIntensity( int(Intensity - ShortIntensityStep) );
      ILCD->setValue( int(Intensity) );
    }
    break;

  case Qt::Key_Left:
    setSpeakerLeft();
    break;
  case Qt::Key_Right:
    setSpeakerRight();
    break;

  case Qt::Key_Pause:
  case Qt::Key_M:
    toggleMute();
    break;

  default:
    RePro::keyPressEvent( qke );

  }
}


void Search::setIntensity( int intensity )
{
  if ( Intensity == intensity ) 
    return;

  Intensity = intensity;
  if ( Intensity < MinIntensity ) 
    Intensity = MinIntensity;
  if ( Intensity > MaxIntensity ) 
    Intensity = MaxIntensity;
  setNumber( "intensity", Intensity );
}


void Search::setDuration( int duration )
{
  double dur = 0.001*duration; // change to sec
 
  if ( Duration == dur ) 
    return;

  Duration = dur;
  if ( Duration < MinDuration ) 
    Duration = MinDuration;
  if ( Duration > MaxDuration ) 
    Duration = MaxDuration;
  setNumber( "duration", Duration );

  // plot trace:
  plotToggle( true, true, 1.25*Duration, 0.125*Duration );

  // new stimulus:
  NewSignal = true;
}


void Search::setPause( int pause )
{
  double pdur = 0.001*pause; // change to sec

  if ( Pause == pdur ) 
    return;

  Pause = pdur;
  if ( Pause < MinPause ) 
    Pause = MinPause;
  if ( Pause > MaxPause ) 
    Pause = MaxPause;
  setNumber( "pause", Pause );
}


void Search::setFrequency( int freq )
{
  double f = freq;

  if ( Frequency == f ) 
    return;

  Frequency = f;
  if ( Frequency < MinFrequency ) 
    Frequency = MinFrequency;
  if ( Frequency > MaxFrequency ) 
    Frequency = MaxFrequency;
  setNumber( "frequency", Frequency );

  // new stimulus:
  NewSignal = true;
}


void Search::setWaveform( int wave )
{
  if ( Waveform == wave ) 
    return;

  Waveform = wave;
  if ( Waveform == 1 )
    selectText( "waveform", "noise" );
  else
    selectText( "waveform", "sine" );

  // new stimulus:
  NewSignal = true;
}


void Search::setWaveformButton( int wave )
{
  if ( wave == 1 )
    NoiseButton->setChecked( true );
  else
    SineButton->setChecked( true );
}


void Search::setSpeaker( bool left )
{
  if ( ! left ) {
    setSpeakerRight();
  }
  else {
    setSpeakerLeft();
  }
}


void Search::setSpeakerLeft( void )
{
  SearchLeft = true;
  if ( SetBestSide + ( sessionRunning() ? 0 : 1 ) > 1 )
    metaData( "Cell" ).selectText( "best side", "left" );
  LeftButton->setChecked( true );
  RightButton->setChecked( false );
  selectText( "side", "left" );
}


void Search::setSpeakerRight( void )
{
  SearchLeft = false;
  if ( SetBestSide + ( sessionRunning() ? 0 : 1 ) > 1 )
    metaData( "Cell" ).selectText( "best side", "right" );
  LeftButton->setChecked( false );
  RightButton->setChecked( true );
  selectText( "side", "right" );
}


void Search::toggleMute( void )
{
  setMute( ! Mute );
}


void Search::setMute( bool mute )
{
  if ( mute != Mute ) {
    if ( ! mute ) {
      Mute = false;
      MuteButton->setDown( false );
      NewSignal = true;
    }
    else {
      Mute = true;
      MuteButton->setDown( true );
      NewSignal = true;
    }
  }
}


void Search::dialogAccepted( void )
{
  ILCD->setValue( int( number( "intensity" ) ) );
  DLCD->setValue( int( number( "duration", "ms" ) ) );
  PLCD->setValue( int( number( "pause", "ms" ) ) );
  FLCD->setValue( int( number( "frequency", "Hz" ) ) );
  setWaveformButton( index( "waveform" ) );
  int side = index( "side" );
  if ( side > 1 )
    side = metaData( "Cell" ).index( "best side" );
  setSpeaker( ( side == 0 ) );
}


void Search::customEvent( QEvent *qce )
{
  ILCD->setValue( int(Intensity) );
  if ( qce->type() - QEvent::User == 11  ) {
    DLCD->setValue( int(1000.0*Duration) );
    PLCD->setValue( int(1000.0*Pause) );
    FLCD->setValue( int(Frequency) );
    setWaveformButton( Waveform );
    int side = index( "side" );
    if ( side > 1 )
      side = metaData( "Cell" ).index( "best side" );
    setSpeaker( ( side == 0 ) );
  }
  else
    RELACSPlugin::customEvent( qce );
}


addRePro( Search );

}; /* namespace auditory */

#include "moc_search.cc"
