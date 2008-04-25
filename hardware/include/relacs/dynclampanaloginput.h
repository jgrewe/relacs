/*****************************************************************************
 *
 * dynclampanaloginput.h
 * Interface for accessing analog input of a daq-board via a dynamic clamp kernel module.
 *
 * RELACS
 * RealTime ELectrophysiological data Acquisition, Control, and Stimulation
 * Copyright (C) 2002-2008 Jan Benda <j.benda@biologie.hu-berlin.de>
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or
 * (at your option) any later version.
 * 
 * RELACS is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 *****************************************************************************/

#ifndef _DYNCLAMPANALOGINPUT_H_
#define _DYNCLAMPANALOGINPUT_H_

#include <vector>
#include <comedilib.h>
#include <relacs/daqerror.h>
#include <relacs/analoginput.h>
#include "comedianaloginput.h"
#include "dynclampanalogoutput.h"
#include "moduledef.h"

using namespace std;


/*! 
\class DynClampAnalogInput
\author Marco Hackenberg
\version 0.1
\brief Interface for accessing analog input of a daq-board via a dynamic clamp kernel module.
*/


class DynClampAnalogInput : public AnalogInput
{

public:

    /*! Device type id for dynamic clamp analog input device. */
  static const int DynClampAnalogInputType = 5;

    /*! Create a new DynClampAnalogInput without opening a device. */
  DynClampAnalogInput( void );
    /*! Constructs an DynClampAnalogInput with device class name \a device
        and type id \a aotype. */
  DynClampAnalogInput( const string &deviceclass );
    /*! Stop analog input and close the daq driver. */
  ~DynClampAnalogInput( void );

    /*! Open the analog input device specified by \a device.
        Returns zero on success, or InvalidDevice (or any other negative number
	indicating the error).
        \sa isOpen(), close(), reset() */
  int open( const string &devicename, long mode=0 );
    /*! Returns true if dynamic clamp module was succesfully opened.
        \sa open(), close(), reset() */
  bool isOpen( void ) const;
    /*! Stop all activity and close the device.
        \sa open(), isOpen(), reset() */
  void close( void );


    /*! Set the name of the dynamic clamp module file. This has to be done 
        before performing prepareRead() or  startRead().
        Returns zero on success, or InvalidDevice (or any other negative number
	indicating the error).        
	\sa moduleName() \sa setDeviceName() \sa deviceName()  */
  int setModuleName( string modulename );

    /*! Return the name of the dynamic clamp module file.
      \sa setModuleName() \sa setDeviceName() \sa deviceName()  */
  string moduleName( void ) const;

    /*! Returns the pointer to the device file.
      \sa setDeviceName() \sa open() \sa subdevice() */
//  comedi_t* device( void ) const;

    /*! Comedi internal index of analog input subdevice. */
  int subdevice( void ) const;

  int maxBufSize ( void ) const;

    /*! Number of analog input channels. */
  int channels( void ) const;

    /*! Resolution in bits of analog input. */
  int bits( void ) const;

    /*! Maximum sampling rate in Hz of analog input. */
  double maxRate( void ) const;

    /*! Maximum number of analog input ranges. */
  int maxRanges( void ) const;
    /*! Voltage range \a index in Volt for unipolar mode.
        If -1 is returned this range is not supported. */
  double unipolarRange( int index ) const;
    /*! Voltage range \a index in Volt for bipolar mode.
        If -1 is returned this range is not supported. */
  double bipolarRange( int index ) const;

    /*! Device driver specific tests on the settings in \a sigs
        for each input signal.
	Before this function is called, the validity of the settings in 
	\a sigs was already tested by testReadData().
	This function should test whether the settings are really supported
	by the hardware.
	If an error ocurred in any trace, the corresponding errorflags in the
	InData are set and a negative value is returned.
        The channels in \a sigs are not sorted.
        This function is called by testRead(). */
  int testReadDevice( InList &sigs );

    /*! Prepare analog input of the input signals \a sigs on the device.
	If an error ocurred in any signal, the corresponding errorflags in
	InData are set and a negative value is returned.
	This function assumes that \a sigs successfully passed testRead().
        The channels in \a sigs are not sorted. */
  int prepareRead( InList &sigs );

    /*! Convert data of the input signals \a sigs.
	If an error ocurred in any channel, the corresponding errorflags in the
	InData structure are filled and a negative value is returned.
	The input signals are sorted by channel number first
        and are then multiplexed into a buffer of signed short's (2 byte).
        The buffer is attached to the first signal in \a sigs. */
  int convertData( InList &sigs );

    /*! Start analog input of the input signals \a sigs
        after they were prepared by prepareRead().
	If an error ocurred in any signal, the corresponding errorflags in
	InData are set and a negative value is returned.
        The channels in \a sigs are not sorted.
	Also start possible pending acquisition on other devices
	that are known from take(). */
  int startRead( InList &sigs );

    /*! Read data to a running data acquisition.
        Returns the number of data values that were popped from the \a trace- 
	device-buffer (sum over all \a traces).
	If an error ocurred in any channel, the corresponding errorflags in the
	InList structure are filled and a negative value is returned.  */
  int readData( InList &sigs );
  
    /*! Read data to a running data acquisition.
        Returns the number of data values that were popped from the \a trace- 
	device-buffer (sum over all \a traces).
	If an error ocurred in any channel, the corresponding errorflags in the
	InList structure are filled and a negative value is returned.  
	For internal usage!
    */
  int fillReadBuffer( void );

    /*! A template function that is used for the implementation
        of the convertData() function.
	This function first sorts the input signals by channel number
	and then multiplexes the signals in a buffer of type \a T
	after appropriate scaling.
        Data values exceeding the range of the daq board are truncated.
        The buffer is attached to the first signal in \a sigs. */
  template < typename T >
    int convert( InList &sigs );

    /*! Stop any running ananlog input activity on the device.
        Returns zero on success, otherwise one of the flags 
        NotOpen, InvalidDevice, ReadError.
        \sa close(), open(), isOpen() */
  int stop ( void );

    /*! Stop any running ananlog input activity and reset the device.
        Returns zero on success, otherwise one of the flags 
        NotOpen, InvalidDevice, ReadError.
        \sa close(), open(), isOpen() */
  int reset( void );

    /* Reloads the prepared configuration commands of the following acquisition 
       into the registers of the hardware after stop() was performed.
       For internal usage.
       \sa stop(), prepareRead() */
  int reload( void );

    /* True, if configuration command for acquisition is successfully loaded
       into the registers of the hardware.
       For internal usage.
       \sa running(), reload() */
  bool loaded( void ) const;

    /*! True if analog input was prepared using testReadDevice() and prepareRead() */
  bool prepared( void ) const;
  
    /*! True if analog input is running. */
  bool running( void ) const;

    /* Sets the running status and unsets the prepared status. */
  void setRunning( void );

    /*! Get error status of the device. 
        0: no error
	-1: underrun
        other: unknown */
  int error( void ) const;

  
    /*! Index of signal start.
        The defualt implemetation returns -1, indicating that
        no index is available.
        If the analog output driver can return
        an index into the data stream of a running analog input
        where the last analog output started,
        then this function should return the this index.
        You also need to reimplement getAISyncDevice()
        to let the user know about this property. */
  virtual long index( void );

    /*! Check for every analog input and input device in \a ais and \a aos
        whether it can be simultaneously started by startRead()
	from this device (\a syncmode = 0)
	or whether the device driver can read the index of an running
	analog input at the time of starting an analog input (\a syncmode = 1).
	Add the indices of those devices to \a aiinx and \a aoinx. */
  virtual void take( vector< AnalogInput* > &ais, vector< AnalogOutput* > &aos,
		                 vector< int > &aiinx, vector< int > &aoinx );




private:

  ComediAnalogInput *CAI;

  // needed by DynClamp class e.g. for assigning TraceInfo strings to channels
  int SubdeviceID;
  bool IsLoaded;
  bool IsKernelDaqOpened;

  string Modulename;
  int Modulefile;

  unsigned int Subdevice;
  int Channels;
  int Bits;
  double MaxRate;
  int ComediBufferSize;
  unsigned int BufferElemSize;  

  unsigned int ChanList[MAXCHANLIST];

  InList *Sigs;

  int ErrorState;
  mutable bool IsRunning;
  bool IsPrepared;

  vector< DynClampAnalogInput* > DynClampAIs;
  vector< DynClampAnalogOutput* > DynClampAOs;
  vector< int > ComediAIsLink;
  vector< int > ComediAOsLink;

};


#endif
