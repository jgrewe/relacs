/*
  inlist.cc
  A container for InData

  RELACS - RealTime ELectrophysiological data Acquisition, Control, and Stimulation
  Copyright (C) 2002-2007 Jan Benda <j.benda@biologie.hu-berlin.de>

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

#include <sstream>
#include "inlist.h"
using namespace std;


InList::InList( void )
  : IL()
{
}


InList::InList( InData &data )
  : IL()
{
  push( data );
}


InList::InList( InData *data, bool own )
  : IL()
{
  add( data, own );
}


InList::InList( const InList &il )
  : IL()
{
  IL.resize( il.IL.size() );
  for ( unsigned int k=0; k<IL.size(); k++ ) {
    IL[k].Own = il.IL[k].Own;
    if ( IL[k].Own )
      IL[k].ID = new InData( *(il.IL[k].ID) );
    else
      IL[k].ID = il.IL[k].ID;
  }
}


InList::~InList( void )
{
  clear();
}


void InList::resize( int n, int m, double step )
{
  if ( n <= 0 ) {
    clear();
    return;
  }

  int os = IL.size();
  IL.resize( n, ILE() );
  if ( n > os ) {
    for ( int k=os; k<n; k++ ) {
      IL[k].ID = new InData( m, step );
      IL[k].Own = true;
    }
  }
}


void InList::clear( void )
{
  for ( unsigned int k=0; k<IL.size(); k++ ) {
    if ( IL[k].Own )
      delete IL[k].ID;
  }
  IL.clear();
}


void InList::reserve( int n )
{
  IL.reserve( n );
}


InList &InList::operator=( const InList &il )
{
  if ( &il == this )
    return *this;

  IL.resize( il.IL.size() );
  for ( unsigned int k=0; k<IL.size(); k++ ) {
    IL[k].Own = il.IL[k].Own;
    if ( IL[k].Own )
      IL[k].ID = new InData( *(il.IL[k].ID) );
    else
      IL[k].ID = il.IL[k].ID;
  }

  return *this;
}


const InData &InList::front( void ) const
{
  return *IL.front().ID;
}


InData &InList::front( void )
{
  return *IL.front().ID;
}


const InData &InList::back( void ) const
{
  return *IL.back().ID;
}


InData &InList::back( void )
{
  return *IL.back().ID;
}


const InData &InList::operator[]( const string &ident ) const
{
  for ( unsigned int k=0; k<IL.size(); k++ ) {
    if ( IL[k].ID->ident() == ident ) 
      return *IL[k].ID;
  }
  return front();
}


InData &InList::operator[]( const string &ident )
{
  for ( unsigned int k=0; k<IL.size(); k++ ) {
    if ( IL[k].ID->ident() == ident ) 
      return *IL[k].ID;
  }
  return front();
}


int InList::index( const string &ident ) const
{
  for ( unsigned int k=0; k<IL.size(); k++ ) {
    if ( IL[k].ID->ident() == ident ) 
      return k;
  }
  return -1;
}


void InList::push( InData &data )
{
  IL.push_back( ILE( new InData( data ), true ) );
}


void InList::push( const InList &traces )
{
  IL.reserve( IL.size() + traces.size() );
  for ( int k=0; k<traces.size(); k++ )
    IL.push_back( ILE( new InData( traces[k] ), true ) );
}


void InList::add( InData *data, bool own )
{
  IL.push_back( ILE( data, own ) );
}


void InList::add( const InData *data, bool own )
{
  IL.push_back( ILE( const_cast<InData*>(data), own ) );
}


void InList::add( const InList &traces, bool own )
{
  IL.reserve( IL.size() + traces.size() );
  for ( int k=0; k<traces.size(); k++ )
    IL.push_back( ILE( const_cast<InData*>(&traces[k]), own ) );
}


void InList::erase( int index )
{
  if ( index >= 0 && index < size() ) {
    if ( IL[index].Own )
      delete IL[index].ID;
    IL.erase( IL.begin() + index );
  }
}


bool lessChannelILE( const InList::ILE &a, const InList::ILE &b )
{
  return ( a.ID->channel() < b.ID->channel() );
}


void InList::sortByChannel( void )
{
  sort( IL.begin(), IL.end(), lessChannelILE );
}


bool lessDeviceChannelILE( const InList::ILE &a, const InList::ILE &b )
{
  if ( a.ID->device() == b.ID->device() )
    return ( a.ID->channel() < b.ID->channel() );
  else
    return ( a.ID->device() < b.ID->device() );
}


void InList::sortByDeviceChannel( void )
{
  sort( IL.begin(), IL.end(), lessDeviceChannelILE );
}


void InList::clearBuffer( void )
{
  for ( int k=0; k<size(); k++ )
    operator[]( k ).clearBuffer();
}


void InList::setDevice( int device )
{
  for ( int k=0; k<size(); k++ )
    operator[]( k ).setDevice( device );
}


void InList::setReference( InData::RefType ref )
{
  for ( int k=0; k<size(); k++ )
    operator[]( k ).setReference( ref );
}


void InList::setDither( bool dither )
{
  for ( int k=0; k<size(); k++ )
    operator[]( k ).setDither( dither );
}


void InList::setUnipolar( bool unipolar )
{
  for ( int k=0; k<size(); k++ )
    operator[]( k ).setUnipolar( unipolar );
}


void InList::setStartSource( int startsource )
{
  for ( int k=0; k<size(); k++ )
    operator[]( k ).setStartSource( startsource );
}


void InList::setDelay( double delay )
{
  for ( int k=0; k<size(); k++ )
    operator[]( k ).setDelay( delay );
}


void InList::setPriority( bool priority )
{
  for ( int k=0; k<size(); k++ )
    operator[]( k ).setPriority( priority );
}


void InList::setSampleRate( double rate )
{
  for ( int k=0; k<size(); k++ )
    operator[]( k ).setSampleRate( rate );
}


void InList::setSampleInterval( double step )
{
  for ( int k=0; k<size(); k++ )
    operator[]( k ).setSampleInterval( step );
}


void InList::setContinuous( bool continuous )
{
  for ( int k=0; k<size(); k++ )
    operator[]( k ).setContinuous( continuous );
}


void InList::setScale( double scale )
{
  for ( int k=0; k<size(); k++ )
    operator[]( k ).setScale( scale );
}


void InList::setOffset( double offset )
{
  for ( int k=0; k<size(); k++ )
    operator[]( k ).setOffset( offset );
}


void InList::setUnit( const string &unit )
{
  for ( int k=0; k<size(); k++ )
    operator[]( k ).setUnit( unit );
}


void InList::setUnit( double scale, double offset, const string &unit )
{
  for ( int k=0; k<size(); k++ )
    operator[]( k ).setUnit( scale, offset, unit );
}


void InList::setUpdateTime( double time )
{
  for ( int k=0; k<size(); k++ )
    operator[]( k ).setUpdateTime( time );
}


void InList::clearMode( void )
{
  for ( int k=0; k<size(); k++ )
    operator[]( k ).clearMode();
}


void InList::setMode( int flags )
{
  for ( int k=0; k<size(); k++ )
    operator[]( k ).setMode( flags );
}


void InList::addMode( int flags )
{
  for ( int k=0; k<size(); k++ )
    operator[]( k ).addMode( flags );
}


void InList::delMode( int flags )
{
  for ( int k=0; k<size(); k++ )
    operator[]( k ).delMode( flags );
}


void InList::setSignalIndex( int index )
{
  for ( int k=0; k<size(); k++ )
    operator[]( k ).setSignalIndex( index );
}


void InList::setSignalTime( double time )
{
  for ( int k=0; k<size(); k++ )
    operator[]( k ).setSignalTime( time );
}


void InList::setRestart( void )
{
  for ( int k=0; k<size(); k++ )
    operator[]( k ).setRestart();
}


string InList::errorText( void ) const
{
  ostringstream ss;

  // common errors:
  long long flags = 0xffffffffffffffffLL;
  for ( int k=0; k<size(); k++ )
    flags &= operator[]( k ).error();
  if ( flags > 0 )
    ss << DaqError::errorText( flags ) << ". ";

  // individual errors:
  for ( int k=0; k<size(); k++ ) {
    long long f = operator[]( k ).error() & (~flags);
    if ( f > 0 || !operator[]( k ).errorStr().empty() ) {
      ss << "Channel " << operator[]( k ).channel()
	 << " on device " << operator[]( k ).device() << ": ";
      if ( !operator[]( k ).errorText( f ).empty() )
	ss << operator[]( k ).errorText( f );
      if ( !operator[]( k ).errorText( f ).empty() &&
	   !operator[]( k ).errorStr().empty() )
	ss << ", ";
      if ( !operator[]( k ).errorStr().empty() )
	ss << operator[]( k ).errorStr();
      ss << ". ";
    }
  }
  return ss.str();
}


void InList::clearError( void )
{
  for ( int k=0; k<size(); k++ )
    operator[]( k ).clearError();
}


void InList::setError( long long flags )
{
  for ( int k=0; k<size(); k++ )
    operator[]( k ).setError( flags );
}


void InList::addError( long long flags )
{
  for ( int k=0; k<size(); k++ )
    operator[]( k ).addError( flags );
}


void InList::delError( long long flags )
{
  for ( int k=0; k<size(); k++ )
    operator[]( k ).delError( flags );
}


void InList::addDaqError( int de )
{
  for ( int k=0; k<size(); k++ )
    operator[]( k ).addDaqError( de );
}


void InList::setErrorStr( const string &strg )
{
  for ( int k=0; k<size(); k++ )
    operator[]( k ).setErrorStr( strg );
}


void InList::addErrorStr( const string &strg )
{
  for ( int k=0; k<size(); k++ )
    operator[]( k ).addErrorStr( strg );
}


void InList::setErrorStr( int errnum )
{
  for ( int k=0; k<size(); k++ )
    operator[]( k ).setErrorStr( errnum );
}


void InList::addErrorStr( int errnum )
{
  for ( int k=0; k<size(); k++ )
    operator[]( k ).addErrorStr( errnum );
}


bool InList::success( void ) const
{
  for ( int k=0; k<size(); k++ ) {
    if ( operator[]( k ).failed() )
      return false;
  }
  return true;
}


bool InList::failed( void ) const
{
  for ( int k=0; k<size(); k++ ) {
    if ( operator[]( k ).failed() )
      return true;
  }
  return false;
}


ostream &operator<< ( ostream &str, const InList &data )
{
  for ( int k=0; k<data.size(); k++ ) {
    str << "InData " << k << ":" << '\n';
    str << data[k] << '\n';
  }
  return str;
}


void InList::freeDeviceBuffer( void )
{
  for ( int k=0; k<size(); k++ )
    operator[]( k ).freeDeviceBuffer();
}


void InList::clearDeviceBuffer( void )
{
  for ( int k=0; k<size(); k++ )
    operator[]( k ).clearDeviceBuffer();
}

