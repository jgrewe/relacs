/*
  attenuator.h
  Interface for programming an attenuator.

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

#ifndef _RELACS_ATTENUATOR_H_
#define _RELACS_ATTENUATOR_H_ 1


#include <vector>
#include <relacs/device.h>
using namespace std;

namespace relacs {


/*! 
\class Attenuator
\author Jan Benda
\version 1.1
\brief Interface for programming an attenuator.

The Attenuator class defines an interface for accessing attenuator devices
which are used to attenuate output signals.
You have to reimplement open(), isOpen(), close(), lines(), minLevel(),
maxLevel(), levels(), attenuate(), testAttenuate() for a specific attenuator.
The open(), isOpen(), and close() functions are defined in the Device class.

Via a constructor or the open()-function a specific attenuator device
is assigned to the Attenuator-class and opened.
The isOpen()-function checks whether the device driver
for the attenuator is valid and opened.
With close() the device driver for the attenuator is closed.

The number of output lines that can be attenuated by the attenuator
is returned by lines().
The minimum and maximum possible attenuation level is returned by minLevel()
and maxLevel(), respectively.
levels() returns all possible attenuation levels.

The attenuation level of an output line can be set by the attenuate()-function.
Attenuators usually can be set to discrete attenuation levels only.
If you request a specific attenuation level by attenuate(),
than it is very likely that the actually set attenuation level
differs slightly from the requested one.
To make this attenuation level known to the user,
attenuate() sets the \a decibel variable to the actually set attenuation level.

To check whether a requested attenuation level
is possible and what level would be set without actually setting it,
you can use the testAttenuate()-function.

The mute()-function can be used to mute an output line
and testMute() checks whether the output line can be muted.
Alternatively, you can call attenuate() with \a decibel set to MuteAttenuationLevel.

Possible return values of attenuate(), test() and mute() are:
- \c 0: success
- \c NotOpen: The device driver for the attenuator is not open.
- \c InvalidDevice: An invalid device index is requested, i.e.
     the requested output line is not supported by the attenuator device.
- \c ReadError: Failed in reading the attenuation level.
- \c WriteError: Failed in setting the attenuation level.
- \c Underflow: The requested attenuation level is too high,
     i.e. the requested signal amplitude is too small.
- \c Overflow: The requested attenuation level is too low,
     i.e. the requested signal amplitude is too large.

In case you want to use a attenuator device within RELACS, your
Attenuator implementation needs to provide a void default constructor
(i.e. with no parameters) that does not open the device.  Also,
include the header file \c <relacs/relacsplugin.h> and make the
Attenuator device known to RELACS with the \c addAttenuator(
ClassNameOfYourAttenuatorImplementation, PluginSetName ) macro.
*/


class Attenuator : public Device
{

public:

    /*! Return code indicating a too high requested attenuation level,
        i.e. the requested signal amplitude is too small. */
  static const int Underflow = -7;
    /*! Return code indicating a too low requested attenuation level,
        i.e. the requested signal amplitude is too large. */
  static const int Overflow = -8;

    /*! The attenuation level for which the outputline is muted. */
  static const double MuteAttenuationLevel;

    /*! Constructs an Attenuator. */
  Attenuator( void );
    /*! Constructs an Attenuator with device class \a deviceclass. 
        \sa setDeviceClass() */
  Attenuator( const string &deviceclass );
    /*! Destructor. Closes the attenuator device driver. */
  virtual ~Attenuator( void );

    /*! Returns the number of output lines the attenuator device supports.
        Reimplement this  function for a specific attenuator.
        The default implementation returns 1.
        \sa minLevel(), maxLevel(), levels() */
  virtual int lines( void ) const =0;
    /*! Returns the minimum possible attenuation level in decibel.
        This number can be negative, indicating amplification.
        \sa maxLevel(), levels(), lines() */
  virtual double minLevel( void ) const =0;
    /*! Returns the maximum possible attenuation level in decibel.
        \sa minLevel(), levels(), lines() */
  virtual double maxLevel( void ) const =0;
    /*! Returns in \a l all possible attenuation levels
        sorted by increasing attenuation levels (highest last).
        \sa minLevel(), maxLevel(), lines() */
  virtual void levels( vector<double> &l ) const =0;

    /*! Set the attenuation level of the output line specified by its index \a di
        to \a decibel decibel. 
	\a decibel can also be set to MuteAttenuationLevel for requesting 
	to mute the attenuator.
        Returns the actually set level in \a decibel.
	In case the attenuator was muted, MuteAttenuationLevel is returned in \a decibel.
	If the requested attenuation level is too high or too low 
	(Underflow or Overflow), then the maximum or minimum possible 
        attenuation level is set and returned in \a decibel.
	You have to reimplement this  function for a specific attenuator.
	\return
	- 0 on success
	- \c NotOpen: The device driver for the attenuator is not open.
	- \c InvalidDevice: An invalid device index is requested, i.e.
	  the requested output line is not supported by the attenuator device.
	- \c WriteError: Failed in setting the attenuation level.
	- \c Underflow: The requested attenuation level is too high,
	  i.e. the requested signal amplitude is too small.
	- \c Overflow: The requested attenuation level is too low,
	  i.e. the requested signal amplitude is too large.
        \sa testAtenuate(), mute()
    */
  virtual int attenuate( int di, double &decibel ) =0;
    /*! Tests setting the attenuation level of the output line
        specified by its index \a di to \a decibel decibel. 
	\a decibel can also be set to MuteAttenuationLevel for requesting 
	to mute the attenuator.
        Returns the level that would be set in \a decibel. 
	In case the attenuator would be muted, MuteAttenuationLevel is returned in \a decibel.
	If the requested attenuation level is too high or too low 
	(Underflow or Overflow), then the maximum or minimum possible 
        attenuation level is returned in \a decibel.
	You have to reimplement this  function for a specific attenuator.
	\return
	- 0 on success
	- \c NotOpen: The device driver for the attenuator is not open.
	- \c InvalidDevice: An invalid device index is requested, i.e.
	  the requested output line is not supported by the attenuator device.
	- \c Underflow: The requested attenuation level is too high,
	  i.e. the requested signal amplitude is too small.
	- \c Overflow: The requested attenuation level is too low,
	  i.e. the requested signal amplitude is too large.
       \sa attenuate(), testMute()
    */
  virtual int testAttenuate( int di, double &decibel ) =0;

    /*! Mutes the output line specified by its index \a di. 
	The default implementation simply calls 
	attenuate( di, MuteAttenuationLevel ).
	\return
	- 0 on success
	- \c NotOpen: The device driver for the attenuator is not open.
	- \c InvalidDevice: An invalid device index is requested, i.e.
	  the requested output line is not supported by the attenuator device.
	- \c WriteError: Failed in setting the attenuation level.
       \sa attenuate(), testMute()
     */
  virtual int mute( int di );
    /*! Tests muting the output line specified by its index \a di. 
	The default implementation simply calls 
	testAttenuate( di, MuteAttenuationLevel ).
	\return
	- 0 on success
	- \c NotOpen: The device driver for the attenuator is not open.
	- \c InvalidDevice: An invalid device index is requested, i.e.
	  the requested output line is not supported by the attenuator device.
	- \c WriteError: Failed in setting the attenuation level.
       \sa attenuate(), test(), mute()
     */
  virtual int testMute( int di );


 protected:

    /*! Set the device info().
        Call this function from open().
        \sa info() */
  void setInfo( void );

};


}; /* namespace relacs */

#endif /* ! _RELACS_ATTENUATOR_H_ */
