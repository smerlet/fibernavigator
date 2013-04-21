/////////////////////////////////////////////////////////////////////////////
// Name:            EAPs.cpp
// Author:          athena team-project -> MERLET Sylvain <-
// Creation Date:   10/07/2009
//
// Description: This is the implementation file for EAPs class.
//
// Last modifications:
//      by : SMerlet - 20/04/2013
/////////////////////////////////////////////////////////////////////////////

#include "EAPs.h"

#include "DatasetManager.h"
#include "../Logger.h"
#include "../gui/MyListCtrl.h"
#include "../misc/nifti/nifti1_io.h"
#include "../misc/Fantom/FMatrix.h"

#include <wx/math.h>
#include <wx/xml/xml.h>

#include <algorithm>
#include <complex>
using std::complex;

#include <fstream>
#include <limits>
using std::numeric_limits;

#include <map>
using std::map;
using std::pair;

#include <vector>
using std::vector;


#include <math.h>
// #include <gsl.h>
#include <gsl/gsl_sf_laguerre.h>
#include <boost/math/special_functions/laguerre.hpp>
#include <boost/math/special_functions/gamma.hpp>
#include <boost/math/special_functions/factorials.hpp>


#define DEF_POS   wxDefaultPosition
#define DEF_SIZE  wxDefaultSize



EAPs::EAPs( const wxString &filename )
:   ODFs( filename )
{
    m_order = 0;
    m_sh_basis = SH_BASIS_DIPY;

#ifdef __WXMSW__
    m_name = filename.AfterLast( '\\' );
#else
    m_name = filename.AfterLast( '/' );
#endif
    
    m_type = EAPS;
}

EAPs::~EAPs()
{
    Logger::getInstance()->print( wxT( "Executing EAPs destructor..." ), LOGLEVEL_DEBUG );
// 	shoreDataAranged.clear();

    Logger::getInstance()->print( wxT( "EAPs destructor done." ), LOGLEVEL_DEBUG );
}

bool EAPs::load( nifti_image *pHeader, nifti_image *pBody )
{
    m_columns = pHeader->dim[1];
    m_rows    = pHeader->dim[2];
    m_frames  = pHeader->dim[3];
    m_bands_EAP   = pHeader->dim[4];


    m_voxelSizeX = pHeader->dx;
    m_voxelSizeY = pHeader->dy;
    m_voxelSizeZ = pHeader->dz;
        
    float voxelX = DatasetManager::getInstance()->getVoxelX();
    float voxelY = DatasetManager::getInstance()->getVoxelY();
    float voxelZ = DatasetManager::getInstance()->getVoxelZ();

    if( m_voxelSizeX != voxelX || m_voxelSizeY != voxelY || m_voxelSizeZ != voxelZ )
    {
        Logger::getInstance()->print( wxT( "Voxel size different from anatomy." ), LOGLEVEL_ERROR );
        return false;
    }

    m_type = EAPS;

    

    float *eapData = (float*)pBody->data; //Reange ici une bonne fois pour ttoute

    int nVoxels = m_columns * m_rows * m_frames;
    
	shoreDataAranged.resize(nVoxels * m_bands_EAP, 0);
	 
    // We need to do a bit of moving around with the data in order to have it like we want.
	// This step and the next one could merge together to decrease the computation cost
    for( int i( 0 ); i < nVoxels; ++i )
    {
        for( int j( 0 ); j < m_bands_EAP; ++j )
        {
            shoreDataAranged[i * m_bands_EAP + j] = eapData[j * nVoxels + i];
        }
    }
    
    
	radialOrder_EAP=4;
	angularOrder_EAP=0;
	if ((radialOrder_EAP%2)==0)
	{
		angularOrder_EAP=radialOrder_EAP;
	}
	else{
		angularOrder_EAP=radialOrder_EAP+1;
	}

	
    //double radius=1e-3;
    
    // TODO Sylvain you should use the value that you want. The current value
    // is in m_displayRadius. You should set a default value in the constructor.

    double radius = 1e-3;
	std::vector< float > odfFloatData=shoreToSh(radius);
	std::cout << "odfFloatData.size(): " << odfFloatData.size() << std::endl;

	std::cout << "avant createStructure" << std::endl;
    // Once the file has been read successfully, we need to create the structure 
    // that will contain all the sphere points representing the ODFs.
	
	// TODO Sylvain you are here.
	m_bands=(angularOrder_EAP+1)*(angularOrder_EAP+2)/2;

    createStructure( odfFloatData );

    m_isLoaded = true;

    return true;
}


///////////////////////////////////////////////////////////////////////////
// This function fills up a vector (m_Points) that contains all the point 
// information for all the ODFs. 
//
// i_fileFloatData  : The SHORE coefficients read from the loaded nifti file.
// Returns true if successful, false otherwise.
///////////////////////////////////////////////////////////////////////////
bool EAPs::createStructure( std::vector< float >& shore_data )
{

    return ODFs::createStructure( shore_data );
}

void EAPs::drawGlyph(int zVoxel, int yVoxel, int xVoxel, AxisType axis)
{
    ODFs::drawGlyph(zVoxel, yVoxel, xVoxel, axis);
}


std::vector< float > EAPs::shoreToSh(double radius)
{
	
	double zeta(700); //A inserer dans la fonction
	zeta=1/(4 * pow(M_PI,2) * zeta);
	double x(pow(radius,2)/zeta);

	int nVoxels = m_columns * m_rows * m_frames;

 
	//EAP modeling in Spherical Harmonic basis at a radius R
	unsigned ShNumber=(angularOrder_EAP+1)*(angularOrder_EAP+2)/2; //comparer SH_number et m_band de odf
	std::vector< float > odfFloatData( nVoxels * ShNumber );
	
	//j is the SH index
	unsigned j(0);
	unsigned counter(0);
	
	//SH coefficients computation
	for( int i( 0 ); i < nVoxels; ++i )
	{		
		counter=0;
		for( int n( 0 ); n <= radialOrder_EAP; ++n )
		{
			j=0;
			for( int l( 0 ); l <= n; l+=2 )
			{					
				for( int m( -l ); m <= l; ++m )
				{
					odfFloatData[i * ShNumber + j]+=shoreDataAranged[i * m_bands_EAP + counter]*shoreFunction(n,l,zeta,x);
					counter++;
					j++;
				}
				
			}
		}
	}

    return odfFloatData;
}


double EAPs::shoreFunction(unsigned n, unsigned l, double zeta, double x)
{
	
	
	double res(1);
	res*=pow(-1,n-l/2);
	res*=boost::math::laguerre( n  - l , l + 0.5 , x );
	res*=exp(-x/2);
	res*=kappa(n,l,zeta);
	res*=pow(x,l/2);

	return res;
}


double EAPs::kappa(unsigned n, unsigned l, double zeta)
{
	double res(0);
	if(n-l<0)
	{
		res= sqrt( (2* 1) / (pow(zeta,1.5) * boost::math::tgamma(n + 1.5)) );
	}
	else{
		return sqrt((2* boost::math::factorial<double>(n -l)) / (pow(zeta,1.5) * boost::math::tgamma(n + 1.5)) );
	}
	return res;
}





bool EAPs::updateDisplayRadius()
{
    float newRad = m_pSliderRadius->GetValue() / 10000.0f;

    if( newRad != m_displayRadius )
    {
        m_displayRadius = newRad;
// 		m_displayRadius = 1e-3;
        std::vector< float > odfFloatData=shoreToSh(m_displayRadius);
		createStructure( odfFloatData );
        
        return true;
    }
    
    return false;
}

void EAPs::createPropertiesSizer( PropertiesWindow *pParent )
{
    ODFs::createPropertiesSizer( pParent );
    
    setColorWithPosition( true );
    
    wxBoxSizer *pBoxMain = new wxBoxSizer( wxVERTICAL );
    
    //////////////////////////////////////////////////////////////////////////

    m_pSliderRadius = new MySlider( pParent, wxID_ANY, 50, 0, 100, DEF_POS, wxSize( 150, -1 ), wxSL_HORIZONTAL  );
    
    wxFlexGridSizer *pGridSliders = new wxFlexGridSizer( 2 );

    
    pGridSliders->Add( new wxStaticText( pParent, wxID_ANY, wxT( "EAP Rad" ) ), 0, wxALIGN_RIGHT | wxALIGN_CENTER_VERTICAL | wxALL, 1 );
    pGridSliders->Add( m_pSliderRadius, 0, wxALIGN_LEFT | wxEXPAND | wxALL, 1 );
    
    pBoxMain->Add( pGridSliders, 0, wxEXPAND, 0 );

    //////////////////////////////////////////////////////////////////////////
        
    m_pSliderLightAttenuation->SetValue( m_pSliderLightAttenuation->GetMin() );
    m_pSliderLightXPosition->SetValue( m_pSliderLightXPosition->GetMin() );
    m_pSliderLightYPosition->SetValue( m_pSliderLightYPosition->GetMin() );
    m_pSliderLightZPosition->SetValue( m_pSliderLightZPosition->GetMin() );
    
    //////////////////////////////////////////////////////////////////////////
    
    m_pPropertiesSizer->Add( pBoxMain, 0, wxFIXED_MINSIZE | wxEXPAND, 0 );
    
    pParent->Connect( m_pSliderRadius->GetId(), wxEVT_COMMAND_SLIDER_UPDATED,       wxCommandEventHandler( PropertiesWindow::OnEAPRadiusSliderMoved ) );
}

void EAPs::updatePropertiesSizer()
{

    ODFs::updatePropertiesSizer();
    
    //TODO JC GUI THIS IS NOT NEEDED.
    /*m_pSliderLightAttenuation->Enable( false );
    m_pSliderLightXPosition->Enable( false );
    m_pSliderLightYPosition->Enable( false );
    m_pSliderLightZPosition->Enable( false );
    m_pBtnFlipX->Enable( false );
    m_pBtnFlipY->Enable( false );
    m_pBtnFlipZ->Enable( false );
    
    m_pSliderMinHue->SetValue(     getColor( MIN_HUE )    * 100 );
    m_pSliderMaxHue->SetValue(     getColor( MAX_HUE )    * 100 );
    m_pSliderSaturation->SetValue( getColor( SATURATION ) * 100 );
    m_pSliderLuminance->SetValue(  getColor( LUMINANCE )  * 100 );
    m_pSliderLOD->SetValue(        (int)getLOD() );
    m_pSliderDisplay->SetValue(    getDisplayFactor() );
    m_pSliderScalingFactor->SetValue( getScalingFactor() * 10.0f );
    
    m_pToggleAxisFlipX->SetValue( isAxisFlipped( X_AXIS ) );
    m_pToggleAxisFlipY->SetValue( isAxisFlipped( Y_AXIS ) );
    m_pToggleAxisFlipZ->SetValue( isAxisFlipped( Z_AXIS ) );
    m_pToggleColorWithPosition->SetValue( getColorWithPosition() );*/
    
    //m_psliderScalingFactor->SetValue(m_psliderScalingFactor->GetMin());
    
    /*if( !isDisplayShape( AXIS ) )
    {
        m_pLblThres->Hide();
        m_pSliderFlood->Hide();
        m_pTxtThres->Hide();
        m_pBtnMainDir->Hide();
    }
    else
    {
        m_pLblThres->Show();
        m_pSliderFlood->Show();
        m_pTxtThres->Show();
        m_pBtnMainDir->Show();
    }*/
}
