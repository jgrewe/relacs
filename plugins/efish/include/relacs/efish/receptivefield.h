/*
  efish/receptivefield.h
  Locates the receptive field of a p-unit electrosensory afferent using the mirob roboter.

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

#ifndef _RELACS_EFISH_RECEPTIVEFIELD_H_
#define _RELACS_EFISH_RECEPTIVEFIELD_H_ 1

#include <relacs/repro.h>
#include <relacs/linearrange.h>
#include <relacs/plot.h>
#include <relacs/ephys/traces.h>
#include <relacs/efield/traces.h>
#include <relacs/outdata.h>
#include <relacs/misc/xyzrobot.h>
#include <relacs/point.h>

using namespace relacs;

namespace efish {
/*!
\class ReceptiveField
\brief [RePro] Locates the receptive field of a p-unit electrosensory afferent using the mirob roboter.
\author Jan Grewe
\version 1.0 (Dec 12, 2017)
*/


class ReceptiveField : public RePro,
    public ephys::Traces,
    public efield::Traces
{
  Q_OBJECT

 public:

  ReceptiveField( void );
  virtual int main( void );

 private:
  Plot xPlot, yPlot, posPlot;
  Point safe_pos;
  Point fish_head, fish_tail, reset_position;
  double duration, deltaf, amplitude, pause;
  double best_x, best_y;
  double y_slope;
  int npasses, nfft, nshift;
  misc::XYZRobot *robot = NULL;
  std::vector<int> axis_map, axis_invert;

  int setupRobot( );
  int xScan( OutData &signal );
  int yScan( OutData &signal );
  int scan( OutData &signal, Point &start_pos, Point &end_pos, double speed );
  Point convertToRobotCoords(const Point &fish_point );
  void resetPlots( double xmin, double xmax, double ymin, double ymax );
  //bool rangeSearch( LinearRange &range, double xy_pos, double z_pos,
  //                  std::vector<double> &avg_rates, OutData &signal,
  //		    bool x_search, bool adjust_y, bool simulation );
  void prepareStimulus( OutData &signal, double min, double max, double speed );
  //int presentStimulus( double x_pos, double y_pos, double z_pos,
  //                     int repeat_num, OutData &signal );
  void getSpikes( EventList &spikeTrains );
  void getRate( SampleDataD &rate, const EventData &spike_train, int &start_trial,
                double period, double duration );
  int analyze( EventList &spikeTrains, double &bestt, double &certainty );
  SampleDataD spectrogram( const SampleDataD &rate, int nfft, int nshift );
  bool moveToPosition( double x, double y, double z );
  double getYSlope( void );
};


}; /* namespace efish */

#endif /* ! _RELACS_EFISH_RECEPTIVEFIELD_H_ */
