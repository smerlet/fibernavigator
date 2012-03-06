#include "TrackingWindow.h"
#include "SelectionBox.h"
#include "SelectionEllipsoid.h"
#include "../dataset/Anatomy.h"
#include "../dataset/Fibers.h"
#include "../dataset/ODFs.h"
#include "../dataset/SplinePoint.h"
#include "../dataset/Surface.h"
#include "../dataset/Tensors.h"
#include "../misc/IsoSurface/CIsoSurface.h"
#include "../main.h"

IMPLEMENT_DYNAMIC_CLASS( TrackingWindow, wxScrolledWindow )

BEGIN_EVENT_TABLE( TrackingWindow, wxScrolledWindow )
EVT_PAINT( TrackingWindow::OnPaint )
EVT_SIZE( TrackingWindow::OnSize )
END_EVENT_TABLE()


TrackingWindow::TrackingWindow( wxWindow *pParent, MainFrame *pMf, wxWindowID id, const wxPoint &pos, const wxSize &size )
: wxScrolledWindow( pParent, id, pos, size, wxBORDER_NONE, _T("RTT Canvas") ),
m_pMainFrame( pMf ),
m_pNoteBook( pParent )
{
    SetBackgroundColour( *wxLIGHT_GREY );
    SetCursor( wxCursor( wxCURSOR_HAND ) );
    m_pTrackingSizer = new wxBoxSizer( wxVERTICAL );
    SetSizer( m_pTrackingSizer );
    SetAutoLayout( true );

	//Content of RTT panel
    /********************************/

	m_pBtnSelectFile = new wxButton( this, wxID_ANY,wxT("Tensor not selected"), wxDefaultPosition, wxSize(200, -1) );
    Connect( m_pBtnSelectFile->GetId(), wxEVT_COMMAND_BUTTON_CLICKED, wxCommandEventHandler(TrackingWindow::OnSelectFile) );
	
    m_pBtnStart = new wxToggleButton( this, wxID_ANY,wxT("Start tracking"), wxPoint(0,30), wxSize(140, -1) );
    Connect( m_pBtnStart->GetId(), wxEVT_COMMAND_TOGGLEBUTTON_CLICKED, wxCommandEventHandler(TrackingWindow::OnStartTracking) );
    m_pBtnStart->Enable(false);

    wxImage bmpDelete(MyApp::iconsPath+ wxT("delete.png" ), wxBITMAP_TYPE_PNG);
    wxBitmapButton *m_pbtnDelete = new wxBitmapButton(this, wxID_ANY, bmpDelete, wxPoint(140,30), wxSize(60,-1));
    Connect(m_pbtnDelete->GetId(),wxEVT_COMMAND_BUTTON_CLICKED, wxTreeEventHandler(TrackingWindow::OnClearBox));

    m_pTextFA = new wxStaticText( this, wxID_ANY, wxT("Min FA"), wxPoint(0,60), wxSize(60, -1), wxALIGN_RIGHT );
    m_pSliderFA = new MySlider( this, wxID_ANY, 0, 1, 50, wxPoint(60,60), wxSize(100, -1), wxSL_HORIZONTAL | wxSL_AUTOTICKS );
    m_pSliderFA->SetValue( 10 );
    Connect( m_pSliderFA->GetId(), wxEVT_COMMAND_SLIDER_UPDATED, wxCommandEventHandler(TrackingWindow::OnSliderFAMoved) );
    m_pTxtFABox = new wxTextCtrl( this, wxID_ANY, wxT("0.10"), wxPoint(160,60), wxSize(55, -1), wxTE_CENTRE | wxTE_READONLY );

    m_pTextAngle = new wxStaticText( this, wxID_ANY, wxT("Max angle"), wxPoint(0,90), wxSize(60, -1), wxALIGN_RIGHT );
    m_pSliderAngle = new MySlider( this, wxID_ANY, 0, 1, 90, wxPoint(60,90), wxSize(100, -1), wxSL_HORIZONTAL | wxSL_AUTOTICKS );
    m_pSliderAngle->SetValue( 60 );
    Connect( m_pSliderAngle->GetId(), wxEVT_COMMAND_SLIDER_UPDATED, wxCommandEventHandler(TrackingWindow::OnSliderAngleMoved) );
    m_pTxtAngleBox = new wxTextCtrl( this, wxID_ANY, wxT("60.0 "), wxPoint(160,90), wxSize(55, -1), wxTE_CENTRE | wxTE_READONLY );

    m_pTextStep = new wxStaticText( this, wxID_ANY, wxT("Step"), wxPoint(0,120), wxSize(60, -1), wxALIGN_RIGHT );
    m_pSliderStep = new MySlider( this, wxID_ANY, 0, 5, 20, wxPoint(60,120), wxSize(100, -1), wxSL_HORIZONTAL | wxSL_AUTOTICKS );
    m_pSliderStep->SetValue( 10 );
    Connect( m_pSliderStep->GetId(), wxEVT_COMMAND_SLIDER_UPDATED, wxCommandEventHandler(TrackingWindow::OnSliderStepMoved) );
    m_pTxtStepBox = new wxTextCtrl( this, wxID_ANY, wxT("1.0 mm"), wxPoint(160,120), wxSize(55, -1), wxTE_CENTRE | wxTE_READONLY );

    m_pTextPuncture = new wxStaticText( this, wxID_ANY, wxT("Puncture"), wxPoint(0,150), wxSize(60, -1), wxALIGN_RIGHT );
    m_pSliderPuncture = new MySlider( this, wxID_ANY, 0, 0, 10, wxPoint(60,150), wxSize(100, -1), wxSL_HORIZONTAL | wxSL_AUTOTICKS );
    m_pSliderPuncture->SetValue( 2 );
    Connect( m_pSliderPuncture->GetId(), wxEVT_COMMAND_SLIDER_UPDATED, wxCommandEventHandler(TrackingWindow::OnSliderPunctureMoved) );
    m_pTxtPunctureBox = new wxTextCtrl( this, wxID_ANY, wxT("0.2"), wxPoint(160,150), wxSize(55, -1), wxTE_CENTRE | wxTE_READONLY );

	m_pTextMinLength = new wxStaticText( this, wxID_ANY, wxT("Min length"), wxPoint(0,180), wxSize(60, -1), wxALIGN_RIGHT );
    m_pSliderMinLength = new MySlider( this, wxID_ANY, 0, 0, 400, wxPoint(60,180), wxSize(100, -1), wxSL_HORIZONTAL | wxSL_AUTOTICKS );
    m_pSliderMinLength->SetValue( 10 );
    Connect( m_pSliderMinLength->GetId(), wxEVT_COMMAND_SLIDER_UPDATED, wxCommandEventHandler(TrackingWindow::OnSliderMinLengthMoved) );
    m_pTxtMinLengthBox = new wxTextCtrl( this, wxID_ANY, wxT("10 mm"), wxPoint(160,180), wxSize(55, -1), wxTE_CENTRE | wxTE_READONLY );

	m_pTextMaxLength = new wxStaticText( this, wxID_ANY, wxT("Max length"), wxPoint(0,210), wxSize(60, -1), wxALIGN_RIGHT );
    m_pSliderMaxLength = new MySlider( this, wxID_ANY, 0, 0, 300, wxPoint(60,210), wxSize(100, -1), wxSL_HORIZONTAL | wxSL_AUTOTICKS );
    m_pSliderMaxLength->SetValue( 200 );
    Connect( m_pSliderMaxLength->GetId(), wxEVT_COMMAND_SLIDER_UPDATED, wxCommandEventHandler(TrackingWindow::OnSliderMaxLengthMoved) );
    m_pTxtMaxLengthBox = new wxTextCtrl( this, wxID_ANY, wxT("200 mm"), wxPoint(160,210), wxSize(55, -1), wxTE_CENTRE | wxTE_READONLY );

    wxToggleButton *m_pToggleRandom = new wxToggleButton( this, wxID_ANY,wxT("Use random seeds"), wxPoint(0,240), wxSize(140, -1) );
    Connect( m_pToggleRandom->GetId(), wxEVT_COMMAND_TOGGLEBUTTON_CLICKED, wxCommandEventHandler(TrackingWindow::OnRandomSeeding) );

    //wxToggleButton *m_pToggleInterp = new wxToggleButton( m_pTrackingWindow, wxID_ANY,wxT("Interpolation"), wxPoint(0,270), wxSize(140, -1) );
    //m_pTrackingWindow->Connect( m_pToggleInterp->GetId(), wxEVT_COMMAND_TOGGLEBUTTON_CLICKED, wxCommandEventHandler(TrackingWindow::OnInterpolate) );

}

void TrackingWindow::OnSize( wxSizeEvent &WXUNUSED(event) )
{

}

void TrackingWindow::OnPaint( wxPaintEvent &WXUNUSED(event) )
{
    wxPaintDC dc( this );
}

wxSizer* TrackingWindow::getWindowSizer()
{
    return m_pTrackingSizer;
}


void TrackingWindow::OnStartTracking( wxCommandEvent& WXUNUSED(event) )
{
    m_pMainFrame->m_pDatasetHelper->m_isRTTReady = !m_pMainFrame->m_pDatasetHelper->m_isRTTReady;
    m_pMainFrame->m_pDatasetHelper->m_isRTTDirty = true;

    if( !m_pMainFrame->m_pDatasetHelper->m_isRTTReady )
    {
        m_pMainFrame->m_pMainGL->m_pRealTimeFibers->clearFibersRTT();
        m_pMainFrame->m_pMainGL->m_pRealTimeFibers->clearColorsRTT();
        m_pMainFrame->m_pDatasetHelper->m_isRTTDirty = false;
    }
}

void TrackingWindow::OnClearBox( wxTreeEvent&    event )
{
    m_pMainFrame->onDeleteTreeItem( event );
    m_pMainFrame->m_pMainGL->m_pRealTimeFibers->clearFibersRTT();
    m_pMainFrame->m_pMainGL->m_pRealTimeFibers->clearColorsRTT();
    m_pMainFrame->m_pDatasetHelper->m_isRTTDirty = false;
    m_pMainFrame->m_pDatasetHelper->m_isRTTReady = false;
    m_pBtnStart->SetValue( false );
}

void TrackingWindow::OnSliderFAMoved(wxCommandEvent& WXUNUSED(event))
{
    float sliderValue = m_pSliderFA->GetValue() / 100.0f;
    m_pTxtFABox->SetValue( wxString::Format( wxT( "%.2f"), sliderValue) );
    m_pMainFrame->m_pMainGL->m_pRealTimeFibers->setFAThreshold( sliderValue );
	m_pMainFrame->m_pDatasetHelper->m_isRTTDirty = true;
}

void TrackingWindow::OnSliderAngleMoved( wxCommandEvent& WXUNUSED(event) )
{
    float sliderValue = m_pSliderAngle->GetValue();
    m_pTxtAngleBox->SetValue(wxString::Format( wxT( "%.1f "), sliderValue) );
    m_pMainFrame->m_pMainGL->m_pRealTimeFibers->setAngleThreshold( sliderValue );
	m_pMainFrame->m_pDatasetHelper->m_isRTTDirty = true;
}

void TrackingWindow::OnSliderStepMoved( wxCommandEvent& WXUNUSED(event) )
{
    float sliderValue = m_pSliderStep->GetValue() / 10.0f;
    m_pTxtStepBox->SetValue(wxString::Format( wxT( "%.1f mm"), sliderValue) );
    m_pMainFrame->m_pMainGL->m_pRealTimeFibers->setStep( sliderValue );
	m_pMainFrame->m_pDatasetHelper->m_isRTTDirty = true;
}

void TrackingWindow::OnSliderPunctureMoved( wxCommandEvent& WXUNUSED(event) )
{
    float sliderValue = m_pSliderPuncture->GetValue() / 10.0f;
    m_pTxtPunctureBox->SetValue(wxString::Format( wxT( "%.1f"), sliderValue) );
    m_pMainFrame->m_pMainGL->m_pRealTimeFibers->setPuncture( sliderValue );
	m_pMainFrame->m_pDatasetHelper->m_isRTTDirty = true;
}

void TrackingWindow::OnSliderMinLengthMoved( wxCommandEvent& WXUNUSED(event) )
{
    float sliderValue = m_pSliderMinLength->GetValue();
    m_pTxtMinLengthBox->SetValue(wxString::Format( wxT( "%.1f mm"), sliderValue) );
    m_pMainFrame->m_pMainGL->m_pRealTimeFibers->setMinFiberLength( sliderValue );
	m_pMainFrame->m_pDatasetHelper->m_isRTTDirty = true;
}

void TrackingWindow::OnSliderMaxLengthMoved( wxCommandEvent& WXUNUSED(event) )
{
    float sliderValue = m_pSliderMaxLength->GetValue();
    m_pTxtMaxLengthBox->SetValue(wxString::Format( wxT( "%.1f mm"), sliderValue) );
    m_pMainFrame->m_pMainGL->m_pRealTimeFibers->setMaxFiberLength( sliderValue );
	m_pMainFrame->m_pDatasetHelper->m_isRTTDirty = true;
}

void TrackingWindow::OnSelectFile( wxCommandEvent& WXUNUSED(event) )
{
    //Tensor data
    long item = m_pMainFrame->m_pListCtrl->GetNextItem( -1, wxLIST_NEXT_ALL, wxLIST_STATE_SELECTED );
    Tensors* pTensorInfo = dynamic_cast<Tensors*>((DatasetInfo*)m_pMainFrame->m_pListCtrl->GetItemData( item ));
    
	if( pTensorInfo != NULL )
	{
		m_pBtnSelectFile->SetLabel( pTensorInfo->getName() );

		//Set Step
		float step = m_pMainFrame->m_pDatasetHelper->m_xVoxel / 2.0f;
		m_pSliderStep->SetValue( step * 10.0f );
		m_pTxtStepBox->SetValue( wxString::Format( wxT( "%.1f mm"), step) );
		m_pMainFrame->m_pMainGL->m_pRealTimeFibers->setStep( step );

		//Copy useful matrices
		m_pMainFrame->m_pMainGL->m_pRealTimeFibers->setTensorsMatrix( pTensorInfo->getTensorsMatrix() );
		m_pMainFrame->m_pMainGL->m_pRealTimeFibers->setTensorsFA( pTensorInfo->getTensorsFA() );
		m_pMainFrame->m_pMainGL->m_pRealTimeFibers->setTensorsEV( pTensorInfo->getTensorsEV() );
	}
}

void TrackingWindow::OnRandomSeeding( wxCommandEvent& WXUNUSED(event) )
{
    m_pMainFrame->m_pDatasetHelper->m_isRandomSeeds = !m_pMainFrame->m_pDatasetHelper->m_isRandomSeeds;
    m_pMainFrame->m_pDatasetHelper->m_isRTTDirty = true;
}

void TrackingWindow::OnInterpolate( wxCommandEvent& WXUNUSED(event) )
{
    m_pMainFrame->m_pDatasetHelper->m_interpolateTensors = !m_pMainFrame->m_pDatasetHelper->m_interpolateTensors;
    m_pMainFrame->m_pDatasetHelper->m_isRTTDirty = true;
}


