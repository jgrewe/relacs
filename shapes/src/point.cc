/*
  point.cc
  A point in 3D space.

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

#include <cmath>
#include <iomanip>
#include <relacs/matrix.h>
#include <relacs/point.h>


namespace relacs {


const Point Point::None = Point( NAN, NAN, NAN );

const Point Point::Origin = Point( 0.0, 0.0, 0.0 );

const Point Point::Ones = Point( 1.0, 1.0, 1.0 );

const Point Point::UnitX = Point( 1.0, 0.0, 0.0 );

const Point Point::UnitY = Point( 0.0, 1.0, 0.0 );

const Point Point::UnitZ = Point( 0.0, 0.0, 1.0 );


Point::Point( void )
{
  for ( int k=0; k<Dim; k++ )
    Coords[k] = 0;
}


Point::Point( const Point &p )
{
  for ( int k=0; k<Dim; k++ )
    Coords[k] = p.Coords[k];
}


Point::Point( double x, double y, double z )
{
  Coords[0] = x;
  Coords[1] = y;
  Coords[2] = z;
}


Point::Point( const double c[3] )
{
  for ( int i=0; i<3; i++ )
    Coords[i] = c[i];
}


Point &Point::assign( double x, double y, double z )
{
  Coords[0] = x;
  Coords[1] = y;
  Coords[2] = z;
  return *this;
}


Point &Point::assign( const Point &p )
{
  for ( int k=0; k<Dim; k++ )
    Coords[k] = p[k];
  return *this;
}


Point &Point::operator=( const Point &p )
{
  for ( int k=0; k<Dim; k++ )
    Coords[k] = p[k];
  return *this;
}


Point Point::operator-( void ) const
{
  Point q( *this );
  for ( int k=0; k<Dim; k++ )
    q[k] = -q[k];
  return q;
}


Point Point::operator+( const Point &p ) const
{
  Point q( *this );
  for ( int k=0; k<Dim; k++ )
    q[k] += p[k];
  return q;
}


Point Point::operator-( const Point &p ) const
{
  Point q( *this );
  for ( int k=0; k<Dim; k++ )
    q[k] -= p[k];
  return q;
}


Point Point::operator*( const Point &p ) const
{
  Point q( *this );
  for ( int k=0; k<Dim; k++ )
    q[k] *= p[k];
  return q;
}


Point Point::operator/( const Point &p ) const
{
  Point q( *this );
  for ( int k=0; k<Dim; k++ )
    q[k] /= p[k];
  return q;
}


Point &Point::operator+=( const Point &p )
{
  for ( int k=0; k<Dim; k++ )
    Coords[k] += p[k];
  return *this;
}


Point &Point::operator-=( const Point &p )
{
  for ( int k=0; k<Dim; k++ )
    Coords[k] -= p[k];
  return *this;
}


Point &Point::operator*=( const Point &p )
{
  for ( int k=0; k<Dim; k++ )
    Coords[k] *= p[k];
  return *this;
}


Point &Point::operator/=( const Point &p )
{
  for ( int k=0; k<Dim; k++ )
    Coords[k] /= p[k];
  return *this;
}


Point Point::operator+( double a ) const
{
  Point p( *this );
  for ( int k=0; k<Dim; k++ )
    p[k] += a;
  return p;
}


Point Point::operator-( double a ) const
{
  Point p( *this );
  for ( int k=0; k<Dim; k++ )
    p[k] -= a;
  return p;
}


Point Point::operator*( double a ) const
{
  Point p( *this );
  for ( int k=0; k<Dim; k++ )
    p[k] *= a;
  return p;
}


Point Point::operator/( double a ) const
{
  Point p( *this );
  for ( int k=0; k<Dim; k++ )
    p[k] /= a;
  return p;
}


Point operator+( double a, const Point &p )
{
  Point pp( p );
  for ( int k=0; k<pp.Dim; k++ )
    pp[k] += a;
  return pp;
}


Point operator-( double a, const Point &p )
{
  Point pp( p );
  for ( int k=0; k<pp.Dim; k++ )
    pp[k] = a - pp[k];
  return pp;
}


Point operator*( double a, const Point &p )
{
  Point pp( p );
  for ( int k=0; k<pp.Dim; k++ )
    pp[k] *= a;
  return pp;
}


Point operator/( double a, const Point &p )
{
  Point pp( p );
  for ( int k=0; k<pp.Dim; k++ )
    pp[k] = a / pp[k];
  return pp;
}


Point &Point::operator+=( double a )
{
  for ( int k=0; k<Dim; k++ )
    Coords[k] += a;
  return *this;
}


Point &Point::operator-=( double a )
{
  for ( int k=0; k<Dim; k++ )
    Coords[k] -= a;
  return *this;
}


Point &Point::operator*=( double a )
{
  for ( int k=0; k<Dim; k++ )
    Coords[k] *= a;
  return *this;
}


Point &Point::operator/=( double a )
{
  for ( int k=0; k<Dim; k++ )
    Coords[k] /= a;
  return *this;
}


double Point::dot( const Point &p ) const
{
  double d = 0.0;
  for ( int k=0; k<Dim; k++ )
    d += Coords[k]*p[k];
  return d;
}


Point Point::cross( const Point &p ) const
{
  return Point( Coords[1]*p.Coords[2] - Coords[2]*p.Coords[1],
		Coords[2]*p.Coords[0] - Coords[0]*p.Coords[2],
		Coords[0]*p.Coords[1] - Coords[1]*p.Coords[0] );
}


double Point::magnitude( void ) const
{
  double d = 0.0;
  for ( int k=0; k<Dim; k++ )
    d += Coords[k]*Coords[k];
  return ::sqrt( d );
}


double Point::angle( const Point &p ) const
{
  double d = dot( p );
  d /= magnitude()*p.magnitude();
  return ::acos( d );
}


Point &Point::normalize( void )
{
  double d = 0.0;
  for ( int k=0; k<Dim; k++ )
    d += Coords[k]*Coords[k];
  double f = 1.0/::sqrt( d );
  for ( int k=0; k<Dim; k++ )
    Coords[k] *= f;
  return *this;
}


Point Point::normalized( void ) const
{
  double d = 0.0;
  for ( int k=0; k<Dim; k++ )
    d += Coords[k]*Coords[k];
  double f = 1.0/::sqrt( d );
  return (*this) * f;
}


bool Point::isNone( void ) const
{
  for ( int k=0; k<Dim; k++ ) {
    if ( ::isnan( Coords[k] ) )
      return true;
  }
  return false;
}


bool Point::operator==( const Point &p ) const
{
  for ( int k=0; k<Dim; k++ ) {
    if ( Coords[k] != p[k] )
      return false;
  }
  return true;
}


bool Point::operator!=( const Point &p ) const
{
  for ( int k=0; k<Dim; k++ ) {
    if ( Coords[k] != p[k] )
      return true;
  }
  return false;
}


bool Point::operator<( const Point &p ) const
{
  for ( int k=0; k<Dim; k++ ) {
    if ( Coords[k] >= p[k] )
      return false;
  }
  return true;
}


bool Point::operator<=( const Point &p ) const
{
  for ( int k=0; k<Dim; k++ ) {
    if ( Coords[k] > p[k] )
      return false;
  }
  return true;
}


bool Point::operator>( const Point &p ) const
{
  for ( int k=0; k<Dim; k++ ) {
    if ( Coords[k] <= p[k] )
      return false;
  }
  return true;
}


bool Point::operator>=( const Point &p ) const
{
  for ( int k=0; k<Dim; k++ ) {
    if ( Coords[k] < p[k] )
      return false;
  }
  return true;
}


double Point::distance( const Point &p ) const
{
  double d = 0.0;
  for ( int k=0; k<Dim; k++ ) {
    double c = Coords[k] - p[k];
    d += c*c;
  }
  return ::sqrt( d );
}


Point Point::center( const Point &p ) const
{
  Point middle;
  for ( int k=0; k<Dim; k++ )
    middle[k] = 0.5*(Coords[k] + p[k]);
  return middle;
}


Point Point::min( const Point &p ) const
{
  Point minp( *this );
  for ( int k=0; k<Dim; k++ ) {
    if ( p[k] < minp[k] )
      minp[k] = p[k];
  }
  return minp;
}


Point Point::max( const Point &p ) const
{
  Point maxp( *this );
  for ( int k=0; k<Dim; k++ ) {
    if ( p[k] > maxp[k] )
      maxp[k] = p[k];
  }
  return maxp;
}


Point min( const deque<Point> &pts )
{
  if ( pts.empty() )
    return Point( 0.0, 0.0, 0.0 );
  Point minp = pts.front();
  for ( auto pi=pts.begin()+1; pi != pts.end(); ++pi )
    minp = minp.min( *pi );
  return minp;
}


Point max( const deque<Point> &pts )
{
  if ( pts.empty() )
    return Point( 0.0, 0.0, 0.0 );
  Point maxp = pts.front();
  for ( auto pi=pts.begin()+1; pi != pts.end(); ++pi )
    maxp = maxp.max( *pi );
  return maxp;
}


Point abs( Point p )
{
  for ( int k=0; k<p.Dim; k++ )
    p[k] = ::abs( p[k] );
  return p;
}


Point &Point::operator*=( const Matrix &m )
{
  Point p( *this );
  for ( int i=0; i<3; i++ ) {
    Coords[i] = 0.0;
    for ( int j=0; j<3; j++ )
      Coords[i] += m(i, j) * p[j];
  }
  return *this;
}


ostream &operator<<( ostream &str, const Point &p ) 
{
  str << "x=" << p[0] << ", ";
  str << "y=" << p[1] << ", ";
  str << "z=" << p[2];
  return str;
}


}; /* namespace relacs */