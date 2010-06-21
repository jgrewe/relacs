/*****************************************************************************
 *
 * optwidgetbase.cc
 * 
 *
 * RELACS
 * Relaxed ELectrophysiological data Acquisition, Control, and Stimulation
 * Copyright (C) 2002-2010 Jan Benda <benda@bio.lmu.de>
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

#include <cmath>
#include <QHBoxLayout>
#include <QFileDialog>
#include <relacs/optwidgetbase.h>


namespace relacs {


OptWidgetBase::OptWidgetBase( Options::iterator param, QWidget *label,
			      Options *oo, OptWidget *ow, QMutex *mutex )
  : QObject(),
    Param( param ),
    OO( oo ),
    OW( ow ),
    LabelW( label ),
    W( 0 ),
    UnitLabel( 0 ),
    UnitBrowseW( 0 ),
    OMutex( mutex ),
    Editable( true ),
    ContUpdate( ow->continuousUpdate() ),
    InternChanged( false ),
    InternRead( false )
{
  if ( OW->readOnlyMask() < 0 ||
       ( OW->readOnlyMask() > 0 && ( (*param).flags() & OW->readOnlyMask() ) ) )
    Editable = false;
  else
    (*Param).delFlags( OW->changedFlag() );
  ow->addWidget( this );
}


OptWidgetBase::~OptWidgetBase( void )
{
}

void OptWidgetBase::get( void )
{
}


void OptWidgetBase::reset( void )
{
}


void OptWidgetBase::resetDefault( void )
{
}


void OptWidgetBase::update( void )
{
}


void OptWidgetBase::setMutex( QMutex *mutex )
{
  OMutex = mutex;
}


QWidget *OptWidgetBase::valueWidget( void )
{
  return W;
}


bool OptWidgetBase::editable( void ) const
{
  return Editable;
}


Parameter &OptWidgetBase::param( void )
{
  return *Param;
}


void OptWidgetBase::setUnitLabel( QLabel *l )
{
  if ( l != 0 ) {
    UnitLabel = l;
    UnitBrowseW = 0;
  }
}


void OptWidgetBase::initActivation( void )
{
}


void OptWidgetBase::addActivation( OptWidgetBase *w )
{
  Widgets.push_back( w );
  initActivation();
}


void OptWidgetBase::activateOption( bool eq )
{
  bool ac = (*Param).activation() ? eq : !eq;

  if ( LabelW != 0 )
    LabelW->setEnabled( ac );
  if ( W != 0 )
    W->setEnabled( ac );
  if ( UnitLabel != 0 )
    UnitLabel->setEnabled( ac );
  if ( UnitBrowseW != 0 )
    UnitBrowseW->setEnabled( ac );
}


OptWidgetText::OptWidgetText( Options::iterator param, QWidget *label,
			      Options *oo, OptWidget *ow,
			      QMutex *mutex, QWidget *parent )
  : OptWidgetBase( param, label, oo, ow, mutex ),
    EW( 0 ),
    Value( "" ),
    LW( 0 ),
    BrowseButton( 0 )
{
  if ( Editable ) {
    W = EW = new QLineEdit( (*Param).text( "%s" ).c_str(), parent );
    OptWidget::setValueStyle( W, (*Param).style(), false, true );
    Value = EW->text().toLatin1().data();
    connect( EW, SIGNAL( textChanged( const QString& ) ),
	     this, SLOT( textChanged( const QString& ) ) );
    if ( (*Param).style() & OptWidget::Browse ) {
      BrowseButton = new QPushButton( "Browse...", parent );
      UnitBrowseW = BrowseButton;
      connect( BrowseButton, SIGNAL( clicked( void ) ),
	       this, SLOT( browse( void ) ) );
    }
  }
  else {
    LW = new QLabel( (*Param).text( "%s" ).c_str(), parent );
    OptWidget::setValueStyle( LW, (*Param).style() );
    //    LW->setBackgroundMode( QWidget::PaletteMid );
    LW->setFrameStyle( QFrame::Panel | QFrame::Sunken );
    LW->setLineWidth( 2 );
    W = LW;
  }
}


void OptWidgetText::get( void )
{
  if ( Editable ) {
    InternRead = true;
    bool cn = OO->notifying();
    OO->unsetNotify();
    (*Param).setText( EW->text().toLatin1().data() );
    if ( (*Param).text( 0, "%s" ) != Value )
      (*Param).addFlags( OW->changedFlag() );
    Value = EW->text().toLatin1().data();
    OO->setNotify( cn );
    InternRead = false;
  }
}


void OptWidgetText::reset( void )
{
  InternChanged = true;
  if ( Editable ) {
    EW->setText( (*Param).text( 0, "%s" ).c_str() );
  }
  else {
    LW->setText( (*Param).text( 0, "%s" ).c_str() );
  }
  InternChanged = false;
}


void OptWidgetText::resetDefault( void )
{
  if ( Editable ) {
    InternChanged = true;
    EW->setText( (*Param).defaultText( "%s" ).c_str() );
    InternChanged = false;
  }
}


void OptWidgetText::update( void )
{
  if ( UnitLabel != 0 ) {
    InternChanged = true;
    UnitLabel->setText( (*Param).outUnit().c_str() );
    InternChanged = false;
  }
}


void OptWidgetText::textChanged( const QString &s )
{
  if ( InternRead )
    return;

  if ( ContUpdate && Editable ) {
    if ( InternChanged ) {
      Value = EW->text().toLatin1().data();
      bool cn = OO->notifying();
      OO->unsetNotify();
      (*Param).setText( Value );
      (*Param).delFlags( OW->changedFlag() );
      OO->setNotify( cn );
    }
    else {
      if ( OMutex != 0 )
	OMutex->lock();
      OW->disableUpdate();
      bool cn = OO->notifying();
      OO->unsetNotify();
      (*Param).setText( s.toLatin1().data() );
      if ( (*Param).text( 0, "%s" ) != Value )
	(*Param).addFlags( OW->changedFlag() );
      Value = EW->text().toLatin1().data();
      if ( cn )
	OO->notify();
      (*Param).delFlags( OW->changedFlag() );
      OO->setNotify( cn );
      OW->enableUpdate();
      if ( OMutex != 0 )
	OMutex->unlock();
    }
  }
  for ( unsigned int k=0; k<Widgets.size(); k++ ) {
    Widgets[k]->activateOption( Widgets[k]->param().testActivation( s.toLatin1().data() ) );
  }
}


void OptWidgetText::initActivation( void )
{
  string s = "";
  if ( EW != 0 )
    s = EW->text().toLatin1().data();
  else
    s = LW->text().toLatin1().data();
  Widgets.back()->activateOption( Widgets.back()->param().testActivation( s ) );
}


void OptWidgetText::browse( void )
{
  if ( OMutex != 0 )
    OMutex->lock();
  QFileDialog* fd = new QFileDialog( 0 );
  if ( (*Param).style() & OptWidget::BrowseExisting ) {
    fd->setFileMode( QFileDialog::ExistingFile );
    fd->setWindowTitle( "Open File" );
    fd->setDirectory( Str( (*Param).text( 0 ) ).dir().c_str() );
  }
  else if ( (*Param).style() & OptWidget::BrowseAny ) {
    fd->setFileMode( QFileDialog::AnyFile );
    fd->setWindowTitle( "Save File" );
    fd->setDirectory( Str( (*Param).text( 0 ) ).dir().c_str() );
  }
  else if ( (*Param).style() & OptWidget::BrowseDirectory ) {
    fd->setFileMode( QFileDialog::Directory );
    fd->setWindowTitle( "Choose directory" );
    fd->setDirectory( Str( (*Param).text( 0 ) ).preventSlash().dir().c_str() );
  }
  fd->setFilter( "All (*)" );
  fd->setViewMode( QFileDialog::List );
  if ( OMutex != 0 )
    OMutex->unlock();
  if ( fd->exec() == QDialog::Accepted ) {
    if ( OMutex != 0 )
      OMutex->lock();
    Str filename = "";
    QStringList qsl = fd->selectedFiles();
    if ( qsl.size() > 0 )
      Str filename = qsl[0].toLatin1().data();
    if ( ( (*Param).style() & OptWidget::BrowseAbsolute ) == 0 )
      filename.stripWorkingPath( 3 );
    if ( ( (*Param).style() & OptWidget::BrowseDirectory ) )
      filename.provideSlash();
    OW->disableUpdate();
    bool cn = OO->notifying();
    OO->unsetNotify();
    (*Param).setText( filename );
    if ( (*Param).text( 0 ) != Value )
      (*Param).addFlags( OW->changedFlag() );
    EW->setText( (*Param).text( 0, "%s" ).c_str() );
    Value = EW->text().toLatin1().data();
    if ( cn )
      OO->notify();
    (*Param).delFlags( OW->changedFlag() );
    OO->setNotify( cn );
    OW->enableUpdate();
    if ( OMutex != 0 )
      OMutex->unlock();
  }
}


void OptWidgetText::setUnitLabel( QLabel *l )
{
  if ( l != 0 ) {
    UnitLabel = l;
    UnitBrowseW = l;
  }
}


QPushButton *OptWidgetText::browseButton( void )
{
  return BrowseButton;
}


OptWidgetMultiText::OptWidgetMultiText( Options::iterator param, QWidget *label,
					Options *oo, OptWidget *ow,
					QMutex *mutex, QWidget *parent )
  : OptWidgetBase( param, label, oo, ow, mutex ),
    EW( 0 ),
    CI( 0 ),
    Inserted( false ),
    Update( true ),
    Value( "" ),
    LW( 0 )
{
  if ( Editable ) {
    W = EW = new QComboBox( parent );
    EW->setEditable( ((*Param).style() & OptWidget::SelectText) == 0 );
    if ( ((*Param).style() & OptWidget::SelectText) > 0 )
      OptWidget::setValueStyle( W, (*Param).style(), false, true, true );
    else
      OptWidget::setValueStyle( W, (*Param).style(), false, true );
    EW->setInsertPolicy( QComboBox::InsertAtTop );
    EW->setDuplicatesEnabled( false );
    if ( (*Param).style() & OptWidget::ComboAutoCompletion == 0 )
      EW->setCompleter( 0 );
    reset();
    connect( EW, SIGNAL( editTextChanged( const QString& ) ),
	     this, SLOT( insertText( const QString & ) ) );
    Value = EW->itemText( 0 ).toLatin1().data();
    connect( EW, SIGNAL( currentIndexChanged( const QString& ) ),
	     this, SLOT( textChanged( const QString& ) ) );
    connect( EW, SIGNAL( activated( const QString& ) ),
	     this, SLOT( textChanged( const QString& ) ) );
  }
  else {
    LW = new QLabel( (*Param).text().c_str(), parent );
    OptWidget::setValueStyle( LW, (*Param).style() );
    //    LW->setBackgroundMode( QWidget::PaletteMid );
    LW->setFrameStyle( QFrame::Panel | QFrame::Sunken );
    LW->setLineWidth( 2 );
    W = LW;
  }
}


void OptWidgetMultiText::get( void )
{
  if ( Editable ) {
    InternRead = true;
    bool cn = OO->notifying();
    OO->unsetNotify();
    (*Param).setText( EW->currentText().toLatin1().data() );
    for ( int k=0; k<EW->count(); k++ )
      (*Param).addText( EW->itemText( k ).toLatin1().data() );
    if ( (*Param).text( 0 ) != Value )
      (*Param).addFlags( OW->changedFlag() );
    Value = EW->itemText( 0 ).toLatin1().data();
    OO->setNotify( cn );
    InternRead = false;
  }
}


void OptWidgetMultiText::reset( void )
{
  InternChanged = true;
  if ( Editable ) {
    Update = false;
    EW->clear();
    if ( (*Param).size() > 0 ) {
      string first = (*Param).text( 0 );
      int firstindex = 0;
      for ( int k=0; k<(*Param).size(); k++ ) {
	string s = (*Param).text( k );
	EW->addItem( s.c_str() );
	if ( s == first )
	  firstindex = k;
      }
      if ( firstindex > 0 ) {
	EW->removeItem( 0 );
	EW->setCurrentIndex( firstindex-1 );
      }
      else
	EW->setCurrentIndex( 0 );
    }
    CI = 0;
    Inserted = false;
    Update = true;
  }
  else {
    LW->setText( (*Param).text( 0 ).c_str() );
  }
  InternChanged = false;
}


void OptWidgetMultiText::resetDefault( void )
{
  if ( Editable ) {
    InternChanged = true;
    Update = false;
    EW->setEditText( (*Param).defaultText().c_str() );
    Update = true;
    InternChanged = false;
  }
}


void OptWidgetMultiText::update( void )
{
  if ( UnitLabel != 0 ) {
    InternChanged = true;
    UnitLabel->setText( (*Param).outUnit().c_str() );
    InternChanged = false;
  }
}


void OptWidgetMultiText::textChanged( const QString &s )
{
  if ( InternRead )
    return;

  if ( ContUpdate && Editable && Update) {
    if ( InternChanged ) {
      Value = EW->itemText( 0 ).toLatin1().data();
      bool cn = OO->notifying();
      OO->unsetNotify();
      (*Param).setText( Value );
      for ( int k=0; k<EW->count(); k++ )
	(*Param).addText( EW->itemText( k ).toLatin1().data() );
      (*Param).delFlags( OW->changedFlag() );
      OO->setNotify( cn );
    }
    else {
      if ( OMutex != 0 )
	OMutex->lock();
      OW->disableUpdate();
      bool cn = OO->notifying();
      OO->unsetNotify();
      (*Param).setText( s.toLatin1().data() );
      for ( int k=0; k<EW->count(); k++ )
	(*Param).addText( EW->itemText( k ).toLatin1().data() );
      if ( (*Param).text( 0 ) != Value )
	(*Param).addFlags( OW->changedFlag() );
      Value = EW->itemText( 0 ).toLatin1().data();
      if ( cn )
	OO->notify();
      (*Param).delFlags( OW->changedFlag() );
      OO->setNotify( cn );
      OW->enableUpdate();
      if ( OMutex != 0 )
	OMutex->unlock();
    }
  }
  for ( unsigned int k=0; k<Widgets.size(); k++ ) {
    Widgets[k]->activateOption( Widgets[k]->param().testActivation( s.toLatin1().data() ) );
  }
}


void OptWidgetMultiText::initActivation( void )
{
  string s = "";
  if ( EW != 0 )
    s = EW->currentText().toLatin1().data();
  else
    s = LW->text().toLatin1().data();
  Widgets.back()->activateOption( Widgets.back()->param().testActivation( s ) );
}


void OptWidgetMultiText::insertText( const QString &text )
{
  if ( ! Update )
    return;

  if ( OMutex != 0 )
    OMutex->lock();
  if ( CI == EW->currentIndex() &&
       ( CI > 0 || EW->currentText() != EW->itemText( 0 ) ) ) {
    if ( ! Inserted ) {
      EW->insertItem( 0, EW->currentText() );
      EW->setCurrentIndex( 0 );
      CI = 0;
      Inserted = true;
    }
    else
      EW->setItemText( 0, EW->currentText() );
  }
  else {
    CI = EW->currentIndex();
    if ( Inserted ) {
      if ( EW->count() > 0 ) {
	EW->removeItem( 0 );
	EW->setCurrentIndex( CI-1 );
	CI = EW->currentIndex();
      }
      Inserted = false;
    }
  }
  if ( OMutex != 0 )
    OMutex->unlock();
}


void OptWidgetMultiText::setUnitLabel( QLabel *l )
{
  if ( l != 0 ) {
    UnitLabel = l;
    UnitBrowseW = l;
  }
}


OptWidgetNumber::OptWidgetNumber( Options::iterator param, QWidget *label,
				  Options *oo, OptWidget *ow,
				  QMutex *mutex, QWidget *parent )
  : OptWidgetBase( param, label, oo, ow, mutex ),
    EW( 0 ),
    Value( 0.0 ),
    LW( 0 ),
    LCDW( 0 )
{
  if ( Editable ) {
    double val = (*Param).number( (*Param).outUnit() );
    double min = (*Param).minimum( (*Param).outUnit() );
    double max = (*Param).maximum( (*Param).outUnit() );
    double step = (*Param).step( (*Param).outUnit() );
    W = EW = new DoubleSpinBox( parent );
    EW->setRange( min, max );
    EW->setSingleStep( step );
    if ( (*Param).isNumber() )
      EW->setFormat( (*Param).format() );
    else
      EW->setFormat( "%.0f" );
    EW->setValue( val );
    OptWidget::setValueStyle( W, (*Param).style(), false, true );
    Value = EW->value();
    connect( EW, SIGNAL( valueChanged( double ) ),
	     this, SLOT( valueChanged( double ) ) );
  }
  else {
    if ( (*Param).style() & OptWidget::ValueLCD ) {
      LCDW = new QLCDNumber( parent );
      LCDW->setSegmentStyle( QLCDNumber::Filled );
      LCDW->setSmallDecimalPoint( true );
      LCDW->display( (*Param).text().c_str() );	
      OptWidget::setValueStyle( LCDW, (*Param).style(), true );
      // size:
      if ( ( (*Param).style() & OptWidget::ValueHuge ) == OptWidget::ValueHuge )
	LCDW->setFixedHeight( 16 * LCDW->sizeHint().height() / 10 );
      else if ( (*Param).style() & OptWidget::ValueLarge )
	LCDW->setFixedHeight( 13 * LCDW->sizeHint().height() / 10 );
      else if ( (*Param).style() & OptWidget::ValueSmall )
	LCDW->setFixedHeight( 8 * LCDW->sizeHint().height() / 10 );
      W = LCDW;
    }
    else {
      LW = new QLabel( (*Param).text().c_str(), parent );
      //    LW->setBackgroundMode( QWidget::PaletteMid );
      LW->setAlignment( Qt::AlignRight | Qt::AlignVCenter );
      LW->setFrameStyle( QFrame::Panel | QFrame::Sunken );
      LW->setLineWidth( 2 );
      OptWidget::setValueStyle( LW, (*Param).style() );
      LW->setFixedHeight( LW->sizeHint().height() );
      W = LW;
    }
  }
}


void OptWidgetNumber::get( void )
{
  if ( Editable ) {
    InternRead = true;
    bool cn = OO->notifying();
    OO->unsetNotify();
    (*Param).setNumber( EW->value(), (*Param).outUnit() );
    if ( fabs( (*Param).number( (*Param).outUnit() ) - Value ) > 0.0001*(*Param).step() )
      (*Param).addFlags( OW->changedFlag() );
    Value = EW->value();
    OO->setNotify( cn );
    InternRead = false;
  }
}


void OptWidgetNumber::reset( void )
{
  InternChanged = true;
  if ( Editable ) {
    EW->setValue( (*Param).number( (*Param).outUnit() ) );
  }
  else {
    if ( (*Param).style() & OptWidget::ValueLCD )
      LCDW->display( (*Param).text( "", (*Param).outUnit() ).c_str() );
    else
      LW->setText( (*Param).text( "", (*Param).outUnit() ).c_str() );
  }
  InternChanged = false;
}


void OptWidgetNumber::resetDefault( void )
{
  if ( Editable ) {
    InternChanged = true;
    EW->setValue( (*Param).defaultNumber( (*Param).outUnit() ) );
    InternChanged = false;
  }
}


void OptWidgetNumber::update( void )
{
  InternChanged = true;
  if ( UnitLabel != 0 )
    UnitLabel->setText( (*Param).outUnit().c_str() );
  if ( Editable ) {
    double val = (*Param).number( (*Param).outUnit() );
    double min = (*Param).minimum( (*Param).outUnit() );
    double max = (*Param).maximum( (*Param).outUnit() );
    double step = (*Param).step( (*Param).outUnit() );
    EW->setRange( min, max );
    EW->setSingleStep( step );
    if ( (*Param).isNumber() )
      EW->setFormat( (*Param).format() );
    else
      EW->setFormat( "%.0f" );
    EW->setValue( val );
  }
  InternChanged = false;
}


void OptWidgetNumber::valueChanged( double v )
{
  if ( InternRead )
    return;

  if ( ContUpdate && Editable ) {
    if ( InternChanged ) {
      Value = EW->value();
      bool cn = OO->notifying();
      OO->unsetNotify();
      (*Param).setNumber( Value, (*Param).outUnit() );
      (*Param).delFlags( OW->changedFlag() );
      OO->setNotify( cn );
    }
    else {
      if ( OMutex != 0 )
	OMutex->lock();
      OW->disableUpdate();
      bool cn = OO->notifying();
      OO->unsetNotify();
      (*Param).setNumber( v, (*Param).outUnit() );
      if ( fabs( v - Value ) > 0.0001*(*Param).step() )
	(*Param).addFlags( OW->changedFlag() );
      Value = EW->value();
      if ( cn )
	OO->notify();
      (*Param).delFlags( OW->changedFlag() );
      OO->setNotify( cn );
      OW->enableUpdate();
      if ( OMutex != 0 )
	OMutex->unlock();
    }
  }

  double tol = 0.2;
  if ( (*Param).isNumber() )
    tol = 0.01*(*Param).step();
  for ( unsigned int k=0; k<Widgets.size(); k++ ) {
    Widgets[k]->activateOption( Widgets[k]->param().testActivation( v, tol ) );
  }
}


void OptWidgetNumber::initActivation( void )
{
  double v = 0.0;
  if ( EW == 0 && LCDW == 0 )
    v = (*Param).number( (*Param).outUnit() );
  else {
    if ( EW != 0 )
      v = EW->value();
    else if ( LCDW != 0 )
      v = LCDW->value();
  }
  double tol = 0.2;
  if ( (*Param).isNumber() )
    tol = 0.01*(*Param).step();
  Widgets.back()->activateOption( Widgets.back()->param().testActivation( v, tol ) );
}


void OptWidgetNumber::setUnitLabel( QLabel *l )
{
  if ( l != 0 ) {
    UnitLabel = l;
    UnitBrowseW = l;
  }
}


OptWidgetBoolean::OptWidgetBoolean( Options::iterator param, Options *oo,
				    OptWidget *ow, const string &request,
				    QMutex *mutex, QWidget *parent )
  : OptWidgetBase( param, 0, oo, ow, mutex ),
    EW( 0 ),
    Value( false )
{
  W = new QWidget;
  QHBoxLayout *hb = new QHBoxLayout;
  hb->setContentsMargins( 0, 0, 0, 0 );
  W->setLayout( hb );
  EW = new QCheckBox;
  hb->addWidget( EW );
  OptWidget::setValueStyle( EW, (*Param).style(), false, true );
  QLabel *label = new QLabel( request.c_str() );
  OptWidget::setLabelStyle( label, (*Param).style() );
  hb->addWidget( label );
  reset();
  if ( Editable ) {
    EW->setEnabled( true );
    Value = EW->isChecked();
    connect( EW, SIGNAL( toggled( bool ) ),
	     this, SLOT( valueChanged( bool ) ) );
  }
  else {
    EW->setEnabled( false );
  }
}


void OptWidgetBoolean::get( void )
{
  if ( Editable ) {
    InternRead = true;
    bool cn = OO->notifying();
    OO->unsetNotify();
    (*Param).setBoolean( EW->isChecked() );
    if ( (*Param).boolean( 0 ) != Value )
      (*Param).addFlags( OW->changedFlag() );
    Value = EW->isChecked();
    OO->setNotify( cn );
    InternRead = false;
  }
}


void OptWidgetBoolean::reset( void )
{
  InternChanged = true;
  EW->setChecked( (*Param).boolean() );
  InternChanged = false;
}


void OptWidgetBoolean::resetDefault( void )
{
  if ( Editable ) {
    InternChanged = true;
    EW->setChecked( (*Param).defaultBoolean() );
    InternChanged = false;
  }
}


void OptWidgetBoolean::valueChanged( bool v )
{
  if ( InternRead )
    return;

  if ( ContUpdate && Editable ) {
    if ( InternChanged ) {
      Value = EW->isChecked();
      bool cn = OO->notifying();
      OO->unsetNotify();
      (*Param).setBoolean( Value );
      (*Param).delFlags( OW->changedFlag() );
      OO->setNotify( cn );
    }
    else {
      if ( OMutex != 0 )
	OMutex->lock();
      OW->disableUpdate();
      bool cn = OO->notifying();
      OO->unsetNotify();
      (*Param).setBoolean( v );
      if ( (*Param).boolean( 0 ) != Value )
	(*Param).addFlags( OW->changedFlag() );
      Value = EW->isChecked();
      if ( cn )
	OO->notify();
      (*Param).delFlags( OW->changedFlag() );
      OO->setNotify( cn );
      OW->enableUpdate();
      if ( OMutex != 0 )
	OMutex->unlock();
    }
  }
  string b( v ? "true" : "false" );
  for ( unsigned int k=0; k<Widgets.size(); k++ ) {
    Widgets[k]->activateOption( Widgets[k]->param().testActivation( b ) );
  }
}


void OptWidgetBoolean::initActivation( void )
{
  string b( EW->isChecked() ? "true" : "false" );
  Widgets.back()->activateOption( Widgets.back()->param().testActivation( b ) );
}


OptWidgetDate::OptWidgetDate( Options::iterator param, Options *oo,
			      OptWidget *ow, QMutex *mutex, QWidget *parent )
  : OptWidgetBase( param, 0, oo, ow, mutex ),
    DE( 0 ),
    LW( 0 ),
    Year( 0 ),
    Month( 0 ),
    Day( 0 )
{
  if ( Editable ) {
    Year = (*Param).year( 0 );
    Month = (*Param).month( 0 );
    Day = (*Param).day( 0 );
    W = DE = new QDateEdit( QDate( Year, Month, Day ), parent );
    OptWidget::setValueStyle( W, (*Param).style(), false, true );
    DE->setDisplayFormat( "yyyy-MM-dd" );
    Year = DE->date().year();
    Month = DE->date().month();
    Day = DE->date().day();
    connect( DE, SIGNAL( dateChanged( const QDate& ) ),
	     this, SLOT( valueChanged( const QDate& ) ) );
  }
  else {
    LW = new QLabel( (*Param).text().c_str(), parent );
    OptWidget::setValueStyle( LW, (*Param).style() );
    //    LW->setBackgroundMode( QWidget::PaletteMid );
    LW->setFrameStyle( QFrame::Panel | QFrame::Sunken );
    LW->setLineWidth( 2 );
    W = LW;
  }
}


void OptWidgetDate::get( void )
{
  if ( Editable ) {
    InternRead = true;
    bool cn = OO->notifying();
    OO->unsetNotify();
    (*Param).setDate( DE->date().year(), DE->date().month(), DE->date().day() );
    if ( (*Param).year( 0 ) != Year ||
	 (*Param).month( 0 ) != Month ||
	 (*Param).day( 0 ) != Day )
      (*Param).addFlags( OW->changedFlag() );
    Year = DE->date().year();
    Month = DE->date().month();
    Day = DE->date().day();
    OO->setNotify( cn );
    InternRead = false;
  }
}


void OptWidgetDate::reset( void )
{
  InternChanged = true;
  DE->setDate( QDate( (*Param).year(), (*Param).month(), (*Param).day() ) );
  InternChanged = false;
}


void OptWidgetDate::resetDefault( void )
{
  if ( Editable ) {
    InternChanged = true;
    DE->setDate( QDate( (*Param).defaultYear(),
			(*Param).defaultMonth(),
			(*Param).defaultDay() ) );
    InternChanged = false;
  }
}


void OptWidgetDate::valueChanged( const QDate &date )
{
  if ( InternRead )
    return;

  if ( ContUpdate && Editable ) {
    if ( InternChanged ) {
      Year = DE->date().year();
      Month = DE->date().month();
      Day = DE->date().day();
      bool cn = OO->notifying();
      OO->unsetNotify();
      (*Param).setDate( Year, Month, Day );
      (*Param).delFlags( OW->changedFlag() );
      OO->setNotify( cn );
    }
    else {
      if ( OMutex != 0 )
	OMutex->lock();
      OW->disableUpdate();
      bool cn = OO->notifying();
      OO->unsetNotify();
      (*Param).setDate( date.year(), date.month(), date.day() );
      if ( (*Param).year( 0 ) != Year ||
	   (*Param).month( 0 ) != Month ||
	   (*Param).day( 0 ) != Day )
	(*Param).addFlags( OW->changedFlag() );
      Year = (*Param).year( 0 );
      Month = (*Param).month( 0 );
      Day = (*Param).day( 0 );
      if ( cn )
	OO->notify();
      (*Param).delFlags( OW->changedFlag() );
      OO->setNotify( cn );
      OW->enableUpdate();
      if ( OMutex != 0 )
	OMutex->unlock();
    }
  }
  string s = "";
  if ( DE != 0 )
    s = DE->date().toString( Qt::ISODate ).toLatin1().data();
  else
    s = LW->text().toLatin1().data();
  for ( unsigned int k=0; k<Widgets.size(); k++ ) {
    Widgets[k]->activateOption( Widgets[k]->param().testActivation( s ) );
  }
}


void OptWidgetDate::initActivation( void )
{
  string s = "";
  if ( DE != 0 )
    s = DE->date().toString( Qt::ISODate ).toLatin1().data();
  else
    s = LW->text().toLatin1().data();
  Widgets.back()->activateOption( Widgets.back()->param().testActivation( s ) );
}


OptWidgetTime::OptWidgetTime( Options::iterator param, Options *oo,
			      OptWidget *ow, QMutex *mutex, QWidget *parent )
  : OptWidgetBase( param, 0, oo, ow, mutex ),
    TE( 0 ),
    LW( 0 ),
    Hour( 0 ),
    Minutes( 0 ),
    Seconds( 0 )
{
  if ( Editable ) {
    Hour = (*Param).hour( 0 );
    Minutes = (*Param).minutes( 0 );
    Seconds = (*Param).seconds( 0 );
    W = TE = new QTimeEdit( QTime( Hour, Minutes, Seconds ), parent );
    OptWidget::setValueStyle( W, (*Param).style(), false, true );
    TE->setDisplayFormat( "hh:mm:ss" );
    Hour = TE->time().hour();
    Minutes = TE->time().minute();
    Seconds = TE->time().second();
    connect( TE, SIGNAL( timeChanged( const QTime& ) ),
	     this, SLOT( valueChanged( const QTime& ) ) );
  }
  else {
    LW = new QLabel( (*Param).text().c_str(), parent );
    OptWidget::setValueStyle( LW, (*Param).style() );
    //    LW->setBackgroundMode( QWidget::PaletteMid );
    LW->setFrameStyle( QFrame::Panel | QFrame::Sunken );
    LW->setLineWidth( 2 );
    W = LW;
  }
}


void OptWidgetTime::get( void )
{
  if ( Editable ) {
    InternRead = true;
    bool cn = OO->notifying();
    OO->unsetNotify();
    (*Param).setTime( TE->time().hour(), TE->time().minute(), TE->time().second() );
    if ( (*Param).hour( 0 ) != Hour ||
	 (*Param).minutes( 0 ) != Minutes ||
	 (*Param).seconds( 0 ) != Seconds )
      (*Param).addFlags( OW->changedFlag() );
    Hour = TE->time().hour();
    Minutes = TE->time().minute();
    Seconds = TE->time().second();
    OO->setNotify( cn );
    InternRead = false;
  }
}


void OptWidgetTime::reset( void )
{
  InternChanged = true;
  TE->setTime( QTime( (*Param).hour(), (*Param).minutes(), (*Param).seconds() ) );
  InternChanged = false;
}


void OptWidgetTime::resetDefault( void )
{
  if ( Editable ) {
    InternChanged = true;
    TE->setTime( QTime( (*Param).defaultHour(),
			(*Param).defaultMinutes(),
			(*Param).defaultSeconds() ) );
    InternChanged = false;
  }
}


void OptWidgetTime::valueChanged( const QTime &time )
{
  if ( InternRead )
    return;

  if ( ContUpdate && Editable ) {
    if ( InternChanged ) {
      Hour = TE->time().hour();
      Minutes = TE->time().minute();
      Seconds = TE->time().second();
      bool cn = OO->notifying();
      OO->unsetNotify();
      (*Param).setTime( Hour, Minutes, Seconds );
      (*Param).delFlags( OW->changedFlag() );
      OO->setNotify( cn );
    }
    else {
      if ( OMutex != 0 )
	OMutex->lock();
      OW->disableUpdate();
      bool cn = OO->notifying();
      OO->unsetNotify();
      (*Param).setTime( time.hour(), time.minute(), time.second() );
      if ( (*Param).hour( 0 ) != Hour ||
	   (*Param).minutes( 0 ) != Minutes ||
	   (*Param).seconds( 0 ) != Seconds )
	(*Param).addFlags( OW->changedFlag() );
      Hour = (*Param).hour( 0 );
      Minutes = (*Param).minutes( 0 );
      Seconds = (*Param).seconds( 0 );
      if ( cn )
	OO->notify();
      (*Param).delFlags( OW->changedFlag() );
      OO->setNotify( cn );
      OW->enableUpdate();
      if ( OMutex != 0 )
	OMutex->unlock();
    }
  }
  string s = "";
  if ( TE != 0 )
    s = TE->time().toString( Qt::ISODate ).toLatin1().data();
  else
    s = LW->text().toLatin1().data();
  for ( unsigned int k=0; k<Widgets.size(); k++ ) {
    Widgets[k]->activateOption( Widgets[k]->param().testActivation( s ) );
  }
}


void OptWidgetTime::initActivation( void )
{
  string s = "";
  if ( TE != 0 )
    s = TE->time().toString( Qt::ISODate ).toLatin1().data();
  else
    s = LW->text().toLatin1().data();
  Widgets.back()->activateOption( Widgets.back()->param().testActivation( s ) );
}


OptWidgetLabel::OptWidgetLabel( Options::iterator param, Options *oo,
				OptWidget *ow, QMutex *mutex, QWidget *parent )
  : OptWidgetBase( param, 0, oo, ow, mutex )
{
  string id = ( (*Param).style() & OptWidget::MathLabel ) ?
    (*Param).ident().htmlUnit() : (*Param).ident().html();
  QLabel *l = new QLabel( id.c_str(), parent );
  l->setTextFormat( Qt::RichText );
  l->setAlignment( Qt::AlignLeft );
  l->setWordWrap( false );
  W = l;
  OptWidget::setLabelStyle( W, (*Param).style() );
}


OptWidgetSeparator::OptWidgetSeparator( Options::iterator param, Options *oo,
					OptWidget *ow, QMutex *mutex, QWidget *parent )
  : OptWidgetBase( param, 0, oo, ow, mutex )
{
  QWidget *h = new QWidget( parent );
  QHBoxLayout *layout = new QHBoxLayout;
  layout->setContentsMargins( 0, 0, 0, 0 );
  h->setLayout( layout );
  W = h;
  if ( (*Param).ident().empty() ) {
    QLabel *line = new QLabel;
    layout->addWidget( line );
    line->setFrameStyle( QFrame::HLine | QFrame::Sunken );
    line->setLineWidth( 1 );
    line->setMidLineWidth( 0 );
    line->setFixedHeight( 4 );
  }
  else {
    QLabel *line;
    line = new QLabel;
    layout->addWidget( line );
    line->setFrameStyle( QFrame::HLine | QFrame::Sunken );
    line->setLineWidth( 1 );
    line->setMidLineWidth( 0 );
    line->setFixedHeight( 4 );
    layout->setStretchFactor( line, 10 );
    QLabel *label = new QLabel( (*Param).ident().c_str() );
    layout->addWidget( label );
    OptWidget::setLabelStyle( label, (*Param).style() );
    label->setAlignment( Qt::AlignHCenter | Qt::AlignVCenter );
    int w = label->sizeHint().width();
    label->setFixedWidth( w + 4*w/(*Param).ident().size() );
    layout->setStretchFactor( label, 1 );
    line = new QLabel;
    layout->addWidget( line );
    line->setFrameStyle( QFrame::HLine | QFrame::Sunken );
    line->setLineWidth( 1 );
    line->setMidLineWidth( 0 );
    line->setFixedHeight( 4 );
    layout->setStretchFactor( line, 30 );
  }
}


}; /* namespace relacs */

#include "moc_optwidgetbase.cc"

