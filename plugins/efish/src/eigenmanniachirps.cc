/*
  efish/eigenmanniachirps.cc
  Repro for stimulation with the Eigenmannia-like chirps, i.e. incomplete and complete interruptions. To be used for chripchamber as well as ephys experiments.

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

#include <relacs/efish/eigenmanniachirps.h>
#include <numeric>

using namespace relacs;
using namespace std;
namespace efish {

//************************************************************************************
//********                       Eigenamannia EOD                              *******
EigenmanniaEOD::EigenmanniaEOD():
    sampling_interval(1./20000), eod_model(EODModel::REALISTIC)
{}

EigenmanniaEOD::EigenmanniaEOD( const EODModel eod_model ):
    sampling_interval(1./20000), eod_model(eod_model)
{}

EigenmanniaEOD::EigenmanniaEOD( const EODModel eod_model, const double sampling_interval ):
    sampling_interval(sampling_interval), eod_model(eod_model)
{}

SampleDataD EigenmanniaEOD::getEOD( const double eodf, double &duration, const double phase, 
                                    bool full_cycles ) const {
  if ( full_cycles ) {
    duration = ceil(duration * eodf) / eodf; 
  }
  SampleDataD eod( 0.0, duration, 0.0 );
  if ( eod_model == EODModel::REALISTIC ) {
    for ( size_t i = 0; i < harmonic_group_amplitudes.size(); ++i ) {
      eod += ( sin( 0.0, duration, this->sampling_interval, (i+1) * eodf, harmonic_group_phases[i] + (i+1) * phase ) * harmonic_group_amplitudes[i] );
    }
  } else {
    eod += sin(0.0, duration, this->sampling_interval, eodf, phase);
  }
  return eod;
}

double EigenmanniaEOD::phaseShift( const double eodf, double threshold, bool rising_flank ) const {
  if (eod_model == EODModel::SINE) {
    return 0.0;
  }
  double eod_period = 1./eodf;
  double duration = 3. * eod_period;
  SampleDataD eod = getEOD( eodf, duration ); 

  if ( threshold < min(eod) ) {
    threshold = min(eod) + std::numeric_limits<double>::epsilon();
  } else if (threshold > max(eod) ) {
    threshold = max(eod) - std::numeric_limits<double>::epsilon() ;
  }  

  EventData crossings;
  if ( rising_flank ) {
    rising( eod, crossings, threshold );
  } else {
    falling( eod, crossings, threshold );
  }

  double shift = 0.0;
  if ( crossings.size() > 0 ) {
    shift = 2 * pi() * (crossings[0] - sampling_interval)/eod_period ;
  } else {
    cerr << "EigenmanniaEOD: invalid threshold, could not figure out the phase shift!\n";
  }
  return shift;
}

//************************************************************************************
//********                         Type A Chirp                                *******
TypeAChirp::TypeAChirp( const double sampling_interval) :
    sampling_interval(sampling_interval), eod_model(EODModel::REALISTIC) 
    {}

TypeAChirp::TypeAChirp( const double sampling_interval, const EODModel eod_model ):
    sampling_interval( sampling_interval ), eod_model( eod_model )
{}

void TypeAChirp::eodModel( EODModel model ) {
  this->eod_model = model;
}

EODModel TypeAChirp::eodModel( void ) const {
  return this->eod_model;
}

SampleDataD TypeAChirp::getWaveform( const double eodf, const double chirp_duration, SignalContent signal ) const {
  EigenmanniaEOD eod( this->eod_model, this->sampling_interval );
  double eod_period = 1./eodf;
  double begin_eod_duration = eod_period;
  double end_eod_duration = eod_period;
  
  double start_shift = eod.phaseShift( eodf );
  double stop_shift = 0.0;
  if (signal == SignalContent::FULL) {
    stop_shift = eod.phaseShift( eodf, -10.);
    double delta_phi = stop_shift - start_shift; 
    begin_eod_duration = (delta_phi/(2*eod.pi()) * eod_period);
    end_eod_duration = eod_period + (2*eod.pi() - stop_shift + start_shift)/(2*eod.pi()) * eod_period;
  } else if ( signal == SignalContent::NO_DC ) {
    stop_shift = start_shift + 2 * eod.pi();
  }
  SampleDataD start_eod = eod.getEOD( eodf, begin_eod_duration, start_shift, false );
  SampleDataD end_eod = eod.getEOD( eodf, end_eod_duration, stop_shift, false );
  
  double offset = signal == SignalContent::FULL ? start_eod[start_eod.size() - 1] : 0.0; 
  SampleDataD chirp_data( static_cast<int>( floor( chirp_duration / sampling_interval ) ), 0.0, sampling_interval, offset );

  SampleDataD result = start_eod;
  result = result.append( chirp_data );
  result = result.append( end_eod ); 

  return result;
}


//************************************************************************************
//********                         Type B Chirp                                *******
TypeBChirp::TypeBChirp( const double sampling_interval) :
    sampling_interval(sampling_interval), eod_model(EODModel::REALISTIC) 
    {}

TypeBChirp::TypeBChirp( const double sampling_interval, const EODModel eod_model ):
    sampling_interval( sampling_interval ), eod_model( eod_model )
{}

void TypeBChirp::eodModel( EODModel model ) {
  this->eod_model = model;
}

EODModel TypeBChirp::eodModel( void ) const {
  return this->eod_model;
}

SampleDataD TypeBChirp::getWaveform( const double eodf, const double chirp_duration, SignalContent signal ) const {
  EigenmanniaEOD eod( this->eod_model, this->sampling_interval );
  double eod_period = 1./eodf;
  int interruption_count = static_cast<int>( ceil(chirp_duration/(2*eod_period)) );
  double begin_eod_duration = eod_period;
  double end_eod_duration = eod_period;
  
  double start_shift = eod.phaseShift( eodf );
  double stop_shift = 0.0;
  if (signal == SignalContent::FULL) {
    stop_shift = eod.phaseShift( eodf, -10.);
    double delta_phi = stop_shift - start_shift; 
    begin_eod_duration = (delta_phi/(2*eod.pi()) * eod_period);  
    end_eod_duration = eod_period + (2*eod.pi() - stop_shift + start_shift)/(2*eod.pi()) * eod_period;
  } else if ( signal == SignalContent::NO_DC ) {
    stop_shift = start_shift + 2 * eod.pi();
  }
  
  SampleDataD start_eod = eod.getEOD( eodf, begin_eod_duration, start_shift, false );
  double offset = signal == SignalContent::FULL ? start_eod[start_eod.size() - 1] : 0.0; 
  SampleDataD end_eod = eod.getEOD( eodf, end_eod_duration, stop_shift, false );
  SampleDataD intermediate_eod = eod.getEOD( eodf, eod_period, stop_shift, false );  
  SampleDataD interruption_data( static_cast<int>( floor( eod_period / sampling_interval ) ), 0.0, sampling_interval, offset);
  
  SampleDataD result = start_eod;
  for ( int i = 0; i < interruption_count; ++i ) {
    for ( int j = 0; j < interruption_data.size(); ++j) {
      result.push(interruption_data[j]);
    }
    if (i < interruption_count - 1) {
      for ( int j = 0; j < intermediate_eod.size(); ++j)
          result.push(intermediate_eod[j]);    
    }
  }  
  result = result.append(end_eod);
  return result;
}


//************************************************************************************
//********                     EigenmanniaChirps RePro                         *******
EigenmanniaChirps::EigenmanniaChirps( void )
  : RePro( "EigenmanniaChirps", "efish", "Jan Grewe", "1.0", "May 11, 2020" ) {
  // add some options:
  newSection( "Beat parameter" );
  addNumber( "duration", "Total trial duration", 1.0, 0.001, 100000.0, 0.001, "s", "ms" );
  addNumber( "deltaf", "Difference frequency", 20., 0.1, 1000., 1.0, "Hz" );
  addNumber( "fakefish", "Fake a receiver fish with the given frequency, set to zero to use the real one", 0.0, 1., 1500., 10., "Hz" );
  addNumber( "contrast", "Relative strength of the sender.", 20., 0.01, 200.0, 1.0, "%" );
  addSelection( "eodmodel", "Model for EOD creation.", "sinewave|realistic" );
  addInteger( "repeats", "Number of repeated trials with the same conditions.", 10, 0, 100000, 2 ).setStyle( OptWidget::SpecialInfinite );

  newSection( "Chirps" );
  addSelection( "chirptype", "Type of chirp", "TypeA|TypeB" );
  addNumber( "chirpdelay", "Minimum time until the first chrip", 1.0, 0.0, 1000.0, 0.01, "s");
  addNumber( "chirpduration", "Minimum chirp duration, is extended to integer multiple of EOD period", 0.01, 0.001, 0.5, 0.001, "s", "ms" );
  addNumber( "chirprate", "Rate at which the fake fish generates chirps.", 1.0, 0.001, 100.0, 0.1, "Hz" );
  addSelection( "signaltype", "Type of signal, whether it drives all, only ampullary, or only tuberous pathways", "all|tuberous only|ampullary only" );
   
  newSection( "EOD estimation" );
  addSelection( "intrace", "inputTrace" );
  addBoolean( "usepsd", "Use the power spectrum", true );
  addNumber( "mineodfreq", "Minimum expected EOD frequency", 100.0, 0.0, 10000.0, 10.0, "Hz" ).setActivation( "usepsd", "true" );
  addNumber( "maxeodfreq", "Maximum expected EOD frequency", 2000.0, 0.0, 10000.0, 10.0, "Hz" ).setActivation( "usepsd", "true" );
  addNumber( "eodfreqprec", "Precision of EOD frequency measurement", 1.0, 0.0, 100.0, 0.1, "Hz" ).setActivation( "usepsd", "true" );
  addNumber( "averagetime", "Time for computing EOD frequency", 2.0, 0.0, 100000.0, 1.0, "s" );  
}


bool EigenmanniaChirps::estimateEodFrequency( double &fisheodf ) {
  double averagetime = number( "averagetime" );
  fisheodf = number( "fakefish" );
  if ( fisheodf < 0.1 ) {
    if ( !boolean( "usepsd" ) ) {
      fisheodf = events( EODEvents ).frequency( currentTime() - averagetime, currentTime() );
      if ( EODEvents < 0 ) {
	      warning( "need EOD events of the EOD Trace." );
	      fisheodf = number( "fakefisheodf" );
	      return false;
      }
      return true;
    } else {
      double bigeod = 0.0;
      double bigeodf = 0.0; 
      double min_eodf = number( "mineodfreq" );
      double max_eodf = number( "maxeodfreq" );
      double eodf_prec = number( "eodfreqprec" );
      int intrace = index( "inputTrace" );
      if ( intrace == -1 ) {
        return false;
      }
      int nfft = 1;

      nfft = nextPowerOfTwo( (int)::ceil( 1.0/trace( intrace ).stepsize()/eodf_prec ) );
      eodf_prec = 1.0/trace( intrace ).stepsize()/nfft;
      if ( averagetime < 2.0/trace( intrace ).stepsize()/nfft ) {
	      averagetime = 2.0/trace( intrace ).stepsize()/nfft;
	      warning( "averagetime is too small for requested frequency resolution. I set it to " +
		    Str( averagetime ) + "s for now." );
      }

      SampleDataF data( 0.0, averagetime, trace( intrace ).stepsize() );
      trace( intrace ).copy( currentTime() - averagetime, data );
      SampleDataF power( nfft );
      rPSD( data, power );
      double threshold = power.max( min_eodf, max_eodf );
      EventData powerpeaks( 1000, true );
      peaks( power, powerpeaks, 0.2*threshold );
      double maxpower = 0.0;
      double maxfreq = 0.0;
      for ( int i=0; i<powerpeaks.size(); i++ ) {
	      if ( powerpeaks[i] >= min_eodf && powerpeaks[i] <= max_eodf ) {
	        if ( powerpeaks.eventSize( i ) > maxpower ) {
	          maxpower = powerpeaks.eventSize( i );
	          maxfreq = powerpeaks[i];
	        }
	      }
      }
      if ( bigeod < maxpower ) {
	      bigeod = maxpower;
	      bigeodf = maxfreq;
      }
      fisheodf = bigeodf;
      return true;
    }
  }
  return true;
}


bool EigenmanniaChirps::createStimulus( void ) {
  stimData.clear();
  double sender_eodf = eodf + deltaf;
  EigenmanniaEOD eod( eod_model );
  SampleDataD eod_waveform;
  SampleDataD first_eod_waveform;
  SampleDataD chirp_waveform;
  int chirp_count = static_cast<int>( floor( stimulus_duration * chirp_rate ) );
  double ici = stimulus_duration / chirp_count;
  if ( ici < chirp_duration ) {
    return false; 
  }
  if ( chirp_count * chirp_duration >= stimulus_duration ) {
    return false;
  }
  if ( chirp_type == ChirpType::TYPE_A ){
    TypeAChirp chirp( sampling_interval, eod_model);
    chirp_waveform = chirp.getWaveform( sender_eodf, chirp_duration, signal_content );
  } else {
    TypeBChirp chirp( sampling_interval, eod_model);
    chirp_waveform = chirp.getWaveform( sender_eodf, chirp_duration, signal_content );
  }
  chirp_waveform.save("ChirpWaveform.dat");
  double shift = eod.phaseShift( sender_eodf );
  first_eod_waveform = eod.getEOD( sender_eodf, chirp_delay, shift, false );
  eod_waveform = eod.getEOD( sender_eodf, ici, shift, false );

  stimData = first_eod_waveform;
  for ( int i = 0; i < chirp_count; ++i ){
    stimData.append( chirp_waveform );
    stimData.append( eod_waveform );
  }

  stimData.save("stim_data.dat");
  stimData.setTrace( GlobalEField );
  stimData.setSampleInterval( sampling_interval );
  // stimData.description().addText( "Type", "unrewarded" ).addFlags( OutData::Mutable );
  // stimData.description()["Frequency"].addFlags( OutData::Mutable );
  //stimData.setIdent( ident + "_unrewarded" );
  outList.push( stimData );
  
  //stimData.sineWave();
  return true;
}

void EigenmanniaChirps::readOptions( void ) {
  stimulus_duration = number( "duration", 0.0, "s" );
  deltaf = number( "deltaf", 0.0, "Hz" );
  chirp_duration = number( "chirpduration", 0.0, "s" );
  chirp_rate = number( "chirprate", 0.0, "Hz" );
  chirp_delay = number( "chirpdelay", 0.0, "s" );
  sampling_interval = 1./20000;
  
  string model_selection = text( "eodmodel" );
  if (model_selection == "sinewave") {
    eod_model = EODModel::SINE;
  } else {
    eod_model = EODModel::REALISTIC;
  }
  
  string chirp_selection = text( "chirptype" );
  if ( chirp_selection == "TypeA") {
    chirp_type = ChirpType::TYPE_A;
  } else {
    chirp_type = ChirpType::TYPE_B;
  }
  
  string chirp_location = text( "signaltype" );
  if ( chirp_location == "all" ) {
    signal_content = SignalContent::FULL;
  } else if ( chirp_location == "tuberous only" ) {
    signal_content = SignalContent::NO_DC;
  } else{
    cerr << "ampullary only is not supported yet!" << endl;
    signal_content = SignalContent::FULL;
  }
  
  bool success = estimateEodFrequency( eodf );
  if (!success) {
    cerr << "Could not estimate the fish frequency!" << endl;
  }
}


int EigenmanniaChirps::main( void ) {
  // get options:
  readOptions();
  bool stimulus_ok = createStimulus();
  if (!stimulus_ok) {
    return Failed;
  }


  return Completed;
}


addRePro( EigenmanniaChirps, efish );

}; /* namespace efish */

#include "moc_eigenmanniachirps.cc"
