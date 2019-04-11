/*
  misc/mirob.cc
  The Mirob module linear robot from MPH.

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

#include <math.h>
#include <limits>
#include <iostream>
#include <stdio.h>
#include <sys/types.h>
#include <unistd.h>
#include <TML_lib.h>
#include <relacs/misc/mirob.h>

using namespace std;


namespace misc {


//Public:

Mirob::Mirob( const string &device )
  : Device( "Mirob" )
{
  Opened = false;
  open( device );
}

Mirob::Mirob( void )
  : Device( "Mirob" )
{
  Opened = false;
}

Mirob::~Mirob( void )
{
  close();
}


int Mirob::open( const string &device )
{
  clearError();
  if ( Opened )
    return 0;

  Info.clear();
  Settings.clear();

  Opened = start();

  if ( Opened ) {
    setDeviceName( "Mirob" );
    setDeviceVendor( "MPH" );
    setDeviceFile( device );
    Device::addInfo();
    return 0;
  }
  else
    return 1;
 }


bool Mirob::start( void ) 
{
  int fd = init_mirob();

  if ( fd <= -1 ) {
    return false;
  }

  FileDescr = fd;

  int si = read_setup();
  if ( si < 0 ) {
    close();
    return false;
  }

  setup_axes(si);
  sleep(1);

  return true;
}


void Mirob::close( void ) 
{
  if ( ! Opened )
    return;

  for ( int k=1; k<=3; k++ ) {
    TS_SelectAxis( k );

    if ( ! TS_Stop() )
      std::cerr << "Failed to stop motion for axis " << k << "! "
		<< TS_GetLastErrorText() << std::endl;
    if ( ! TS_Power( POWER_OFF ) )
      std::cerr << "Failed to power off drive for axis " << k << "! "
		<< TS_GetLastErrorText() << std::endl;
  }

  TS_CloseChannel(FileDescr);
  std::cerr << "Channel closed [OK]" << std::endl;

  Opened = false;
}


  int Mirob::reset( void )
  {
    return 0;
  }

/**
   Move the given axis to an absolute position given in mm relative to the
   zero position.

   axis: is the axis you want to move
   pos: the position in mm
   speed: the speed in "mm/s"
   acc: the acceleration in mm/s^2

 */

int Mirob::move( int axis, double pos, double speed, double acc )
{
  if ( pos < 0.0 )
    return -1;
  double stepsize = get_step_length( axis ); // mm per step
  if ( speed == 0.0 )
    speed = this->Speed;
  if ( acc == 0.0 )
    acc = this->Acc;

  double tml_speed = speed / stepsize * this->tml_time;
  double tml_acc = acc / stepsize * this->tml_time * this->tml_time;
  long steps = (long)round( pos / stepsize );

  TS_SelectAxis( axis + 1 );
  TS_MoveAbsolute( steps, tml_speed, tml_acc, UPDATE_IMMEDIATE, FROM_REFERENCE );
  return 0;
}


int Mirob::step( int axis, double s, double speed, double acc)
{
  if ( s == 0.0 )
    return -1;
  double stepsize = get_step_length( axis ); // mm per step
  if ( speed == 0.0 )
    speed = this->Speed;
  if ( acc == 0.0 )
    acc = this->Acc;

  double tml_speed = speed / stepsize * this->tml_time; // steps per controller time
  double tml_acc = acc / stepsize * this->tml_time * this->tml_time;
  long steps = (long)round( s / stepsize );
  bool additive = false;

  TS_SelectAxis( axis + 1 );
  TS_MoveRelative( steps, tml_speed, tml_acc, additive, UPDATE_IMMEDIATE, FROM_REFERENCE );
  return 0;
}


double Mirob::pos( int axis ) const
{
  long int position = 0;
  TS_SelectAxis( axis+1 );
  TS_GetLongVariable( "APOS", position );
  return position;
}


Point Mirob::pos( void ) const 
{
  Point pos;
  for ( int i=0; i<3; i++ ) {
    long int position = 0;
    TS_SelectAxis( i+1 );
    TS_GetLongVariable("APOS", position);
    //std::cerr << "read position from mirob:"<< position <<"(in steps)" << std::endl;
    //std::cerr << "position from mirob:"<< position*get_step_length(i) <<"(in mm)" << std::endl;
    //std::cerr << "size of position: " << sizeof(position) << std::endl;
    //std::cerr << "MaxValue: " << std::numeric_limits<unsigned int>::max() << std::endl;
    // XXX Why does this differ from pos(int)?
    double pos_with_neg = position;
    double maxValue = std::numeric_limits<unsigned int>::max();

    if ( position * get_step_length(i) > 700 ) {
      pos_with_neg = (double) position - maxValue;
    }
    pos[i] = pos_with_neg*get_step_length( i );

  }
  return pos;
}


int Mirob::stop( int axis )
{
  TS_SelectAxis( axis+1 );
  TS_Stop();
  return 0;
}


int Mirob::stop( void )
{
  int r = 0;
  for ( int k=0; k<3; k++ )
    r |= stop( k );
  return r;
}


int Mirob::wait( void ) const
{
  for( int axis=1; axis<=3; axis++ ) {
    TS_SelectAxis( axis );
    if ( ! TS_SetEventOnMotionComplete(1, 0) ) {
      return 1;
    }
  }
  return 0;
}


//*********************************
//Getter and setter functions:
//*********************************

int Mirob::speed( void ) const
{
  return Speed;
}

bool Mirob::setSpeed( int speed ) 
{
  //too high? (for x and y axis)
  if(speed >= 250) {
    Speed = 250;
    return false;
  // too low?
  } else if (speed <=10) {
    Speed = 10;
    return false;

  } else {
    Speed = speed;
    return true;
  }
}


double Mirob::acceleration( void ) const
{
  return Acc;
}


bool Mirob::setAcceleration( double acc )
{
  if ( acc <= 0 ) {
    Acc = 0.1;
    return false;
  } else if ( acc > 3 ) {
    std::cerr << "Too high acc for z??" << std::endl;
    Acc = acc;
    return true;
  } else {
    Acc = acc;
    return true;
  }
}


void Mirob::set_intern_position( int mirobaxis, long int pos )
{
  TS_SelectAxis( mirobaxis );
  TS_SetPosition( pos );
}


double Mirob::get_step_length( int axis ) const
{
  double ds;
  if ( axis < 0 || axis > 2 )
    return -1;
  else {
    if ( axis == 0 ) {
      ds = 10./10000.;
    }
    if ( axis == 1 ) {
      ds = 10./8192.;
    }
    if ( axis == 2 ) {
      ds = 2./8192.;
    }
  }
  return ds;
}

double Mirob::get_axis_factor( int axis ) const
{
  double dflt = 1000.; // x axis: 1000 steps per mm
  double factor;
  if (axis == 1)
    factor = (8192. / 10.) / dflt; // y axis is a bit faster
  else if (axis == 2)
    factor = (8192. / 2.) / dflt; // z axis much slower but finer resolution
  else
    factor = 1.0;
  return factor;
}

double Mirob::get_max(double a, double b, double c)
{
  double max;
  max = (a > b)   ? a : b;
  max = (c > max) ? c : max;
  return max;
}

/**
   the axis goes to the given limit:
   speed: the speed te axis will use
   positive: true the positive limit/false the negative limit is used
 */
void Mirob::search_home( int mirobaxis, int speed, bool positive )
{
  axis_limits limits;
  limits.mirobaxis = mirobaxis;
  check_limit_switch( limits );

  long position = positive ? 1000 : -1000;
  TS_SelectAxis( mirobaxis );

  //Move with the given speed towards the given limit until you reached it.
  while ( positive ? !limits.positive : !limits.negative ) {
    TS_MoveRelative( position, speed, 0.1, false, UPDATE_IMMEDIATE,
		     FROM_REFERENCE );
    check_limit_switch( limits );
  }
  TS_Stop();

  //Move a tiny bit out of the given limit until it is not active anymore.
  while ( positive ? limits.positive : limits.negative ) {
    TS_MoveRelative( -1*position, speed/2, 0.1, false, UPDATE_IMMEDIATE,
		     FROM_REFERENCE );
    check_limit_switch( limits );
  }
  TS_Stop();

  // Set the position of the axis to 0:
  TS_SetPosition( 0 );
}


void Mirob::go_to_reference( bool positive, int speed )
{
  //all axis go to the given limit:
  // positive- true: the positive limit
  //         - false: the negative limit
  search_home(3, speed, positive);
  search_home(2, speed, positive);
  search_home(1, speed, positive);
  std::cerr << "went to reference position" << std::endl;
}


//*********************************
//Init functions:
//*********************************

int Mirob::init_mirob( void )
{
  BYTE channel = 1;
  BYTE cType = CHANNEL_RS232;
  std::string dev = "/dev/ttyS0";
  DWORD brate = 115200;

  int FileDescr = TS_OpenChannel( dev.c_str(), cType, channel, brate );
  if ( FileDescr < 0 ) {
    std::cerr << "Communication error! " << TS_GetLastErrorText() << std::endl;
    return -1;
  } else {
    std::cerr << "Channel opened [OK]" << std::endl;
  }
  return FileDescr;
}


int Mirob::read_setup( void )
{
  int setupindex = TS_LoadSetup( "mirob2.t.zip" );
  if ( setupindex < 0 ) {
    std::cerr << "Failed to load setup file! " << TS_GetLastErrorText()
	      << std::endl;
    return -1;
  } else {
    std::cerr << "Setup loaded [OK]" << std::endl;
  }
  return setupindex;
}


long Mirob::setup_axes( int setupindex )
{
  long int position;
  for (int k = 1; k <= 3; k++) {
    if ( ! TS_SetupAxis( k, setupindex ) )  {
      std::cerr << "Failed to setup axis " << k << "! "
		<< TS_GetLastErrorText() << std::endl;
    }
    TS_SelectAxis( k );
    TS_GetLongVariable("APOS", position);
    printf("Position %i: %dx\n", k, (int) position);
    TS_DriveInitialisation();
    TS_Power( POWER_ON );
  }
  return position;
}


bool Mirob::switch_on_power( void )
{
  WORD SRL_value = 0;
  for ( int k=1; k <=3; k++ ) {
    TS_SelectAxis( k );
    if ( ! TS_Power( POWER_ON ) )
      return false;

    while ( SRL_value == 0 ) {
      if ( ! TS_ReadStatus( REG_SRL, SRL_value ) )
	return false;
      SRL_value = ( ( SRL_value & 1<<15 ) != 0 ? 1 : 0 );
      //SRL.15 - signals the state of the power stage */
    }
  }
  return true;
}


bool Mirob::power_off( void ) {
  for ( int k=1; k <=3; k++ ) {
    TS_SelectAxis( k );
    if ( ! TS_Power( POWER_OFF ) )
      return false;
  }
  int count = 0;
  while ( checkPowerState() && count < 10000 ) {
    count++;
  }
  if ( count >= 10000 ) {
    cerr << "Axis on timeout!" << endl;
    return false;
  }
  return true;
}


bool Mirob::power_on( void ) {
  for ( int k=1; k <=3; k++ ) {
    TS_SelectAxis( k );
    if ( ! TS_Power( POWER_ON ) )
      return false;
  }
  int count = 0;
  while (!checkPowerState() && count < 10000) {
    count++;
  }
  if ( count >= 10000 ) {
    cerr << "Axis on timeout!" << endl;
    return false;
  }
  return true;
}


void Mirob::check_limit_switch( axis_limits &limits )
{
  WORD MER_value = 0;
  TS_SelectAxis( limits.mirobaxis );
  TS_ReadStatus(REG_MER, MER_value);
  limits.positive = (MER_value & 1<<6) != 0;
  limits.negative =(MER_value & 1<<7) != 0;
}


bool Mirob::checkPowerState( int axis ) {
  WORD msr = 0;
  bool on = true;
  if( axis > 2 )
    axis = 2;
  if( axis < 0 ) {
    for (int i = 0; i < 3; i++) {
      TS_SelectAxis( i );
      TS_ReadStatus(REG_MSR, msr);
      on &= ((msr & ( 1 << 13 )) != 0);
    }
  } else {
    TS_SelectAxis( axis );
    TS_ReadStatus(REG_MSR, msr);
    return msr & ( 1 << 13 );
  }
  return on;
}


bool Mirob::check_pos_limit( int mirobaxis )
{
  WORD MER_value = 0;
  TS_SelectAxis( mirobaxis );
  TS_ReadStatus( REG_MER, MER_value );
  return ( MER_value & 1<<6 );
}


bool Mirob::check_neg_limit( int mirobaxis )
{
  WORD MER_value = 0;
  TS_SelectAxis( mirobaxis );
  TS_ReadStatus( REG_MER, MER_value );
  return ( MER_value & 1<<7 );
}


void Mirob::check_all_reg( int mirobaxis )
{
  TS_SelectAxis( mirobaxis );
  short unsigned int status;

  if(! TS_ReadStatus(REG_MCR, status)) {
    std::cerr<< "error reading status MCR." << std::endl;
    std::cerr<< TS_GetLastErrorText() << std::endl;
  } else {
    std::cerr << "MCR status: " << status << std::endl;
  }

  if(! TS_ReadStatus(REG_MSR, status)) {
    std::cerr<< "error reading status MSR." << std::endl;
    std::cerr<< TS_GetLastErrorText() << std::endl;
  } else {
    std::cerr << "MSR status: " << status << std::endl;
  }

  if(! TS_ReadStatus(REG_ISR, status)) {
    std::cerr<< "error reading status ISR." << std::endl;
    std::cerr<< TS_GetLastErrorText() << std::endl;
  } else {
    std::cerr << "ISR status: " << status << std::endl;
  }

  if(! TS_ReadStatus(REG_SRL, status)) {
    std::cerr<< "error reading status SRL." << std::endl;
    std::cerr<< TS_GetLastErrorText() << std::endl;
  } else {
    std::cerr << "SRL status: " << status << std::endl;
  }

  if(! TS_ReadStatus(REG_SRH, status)) {
    std::cerr<< "error reading status SRH." << std::endl;
    std::cerr<< TS_GetLastErrorText() << std::endl;
  } else {
    std::cerr << "SRH status: " << status << std::endl;
  }

  if(! TS_ReadStatus(REG_MER, status)) {
    std::cerr<< "error reading status MER." << std::endl;
    std::cerr<< TS_GetLastErrorText() << std::endl;
  } else {
    std::cerr << "MER status: " << status << std::endl;
  }
}


bool Mirob::toolClamped()
{
  TS_SelectAxis(3);
  BYTE InValue;
  TS_GetInput(36, InValue);
  return InValue == 1;
}

bool Mirob::toolPresent()
{
  TS_SelectAxis(3);
  BYTE InValue;
  TS_GetInput(37, InValue);
  return InValue == 1;
}

void Mirob::toolRelease()
{
  if (toolClamped()) {
    TS_SelectAxis(3);
    TS_SetOutput(30, 0);
    TS_SetOutput(31, 4);
  }
}

void Mirob::toolFix()
{
  if ( !toolClamped() ) {
    TS_SelectAxis(3);
    TS_SetOutput(30, 2);
    TS_SetOutput(31, 0);
  }
}

}; /* namespace misc */
