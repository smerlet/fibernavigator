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
#include "../gfx/ShaderHelper.h"
#include "../gui/MyListCtrl.h"
#include "../gui/SceneManager.h"
#include "../misc/nifti/nifti1_io.h"
#include "../misc/Fantom/FMatrix.h"

#include <GL/glew.h>
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

#define DEF_POS   wxDefaultPosition
#define DEF_SIZE  wxDefaultSize



EAPs::EAPs( const wxString &filename )
:   Glyph(), 
    m_isMaximasSet   ( false ),
    m_axisThreshold  ( 0.5f ),
    m_order          ( 0 ),
    m_radiusAttribLoc( 0 ),
    m_radiusBuffer   ( NULL ),    
    m_nbors          ( NULL ),
	m_sh_basis       ( SH_BASIS_DIPY )
{
    m_scalingFactor = 5.0f;
    m_fullPath = filename;

#ifdef __WXMSW__
    m_name = filename.AfterLast( '\\' );
#else
    m_name = filename.AfterLast( '/' );
#endif
    
    m_type = EAPS;

    // Generating hemispheres
    generateSpherePoints( m_scalingFactor );
// 	
 	odfs=new ODFs(filename); // Ne pas oublier de delete l'odfs dans le destructeur
}

EAPs::~EAPs()
{
    Logger::getInstance()->print( wxT( "Executing EAPs destructor..." ), LOGLEVEL_DEBUG );
    if( m_radiusBuffer )
    {
        glDeleteBuffers( 1, m_radiusBuffer );
        delete [] m_radiusBuffer;
    }

    if( m_nbors != NULL )
    {
        delete [] m_nbors;
        m_nbors = NULL;
    }
    Logger::getInstance()->print( wxT( "EAPs destructor done." ), LOGLEVEL_DEBUG );
}

bool EAPs::load( nifti_image *pHeader, nifti_image *pBody )
{
    m_columns = pHeader->dim[1];
    m_rows    = pHeader->dim[2];
    m_frames  = pHeader->dim[3];
    m_bands   = pHeader->dim[4];

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

    int nVoxels = m_columns * m_rows * m_frames;

    float* eapData = (float*)pBody->data;

	
	std::vector< float > odfFloatData( nVoxels * m_bands );
	
	double radius=1e-3;
	odfFloatData=shoreToSh( eapData, radius, nVoxels, m_bands);

    // Once the file has been read successfully, we need to create the structure 
    // that will contain all the sphere points representing the ODFs.
	
// 	creer un pointeur sur un objet ODF dans le constructeur EAPs et ensuite s'en sevir pour appeler createStrucure
	// TODO Sylvain you are here.
    //odfs->createStructure( odfFloatData ); 
//     createStructure( odfFloatData );
// 	createStructure( l_data);


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
// // 	Generer les ODFs a partir des coeffient SHORE et utiliser un objet ODFs pour la vizualisation.
// // Problem : createStructure de la classe ODFs est private.. :( Demander a JC pour trouver un autre moyen ou simplement 
// // creer une fonction qui prend les SH et fait create structure dans le style de la fonction lODFs::load (une sorte de wrapper)
// 
// draft de la fonction:
// 
// 	//On transforme les coefficnet shore vers des coefficnets sh
// 	sh_data=shoreToSh(shore_data, 0.5)
// 
// 	//Ensuite on met en forme les sh_data pour pouvoir utilser createStructure de la class ODF. Attention createStructure est private => faire un wrapper!
//     // We need to do a bit of moving around with the data in order to have it like we want.
// 	std::vector< float > l_fileFloatData( l_nSize * m_bands );
//     for( int i( 0 ); i < l_nSize; ++i )
//     {
//         for( int j( 0 ); j < m_bands; ++j )
//         {
//             l_fileFloatData[i * m_bands + j] = sh_data[j * l_nSize + i];
//         }
//     }
// 
//     // Once the file has been read successfully, we need to create the structure 
//     // that will contain all the sphere points representing the ODFs.
//     ODFs::createStructure( l_fileFloatData );
// 	
// 	
// 
// 
    return true;
}

void EAPs::drawGlyph(int zVoxel, int yVoxel, int xVoxel, AxisType axis)
{
    odfs->drawGlyph(zVoxel, yVoxel, xVoxel, axis);
}



std::vector< float > EAPs::shoreToSh( float* shoreData, double radius, int nVoxels, int m_bands)
{
// 	Converti les coefficient SHORE vers des coefficients SH

// see https://fr.wikipedia.org/wiki/C%2B%2B11
// dans <math.h>:
//  	double assoc_laguerre( unsigned n, unsigned m, double x ) ;
// 	double laguerre( unsigned n, double x ) ;
// 	exponentielles : exp
// (regarder aussi cmath.h)	


    std::vector< float > odfFloatData( nVoxels * m_bands );

    // We need to do a bit of moving around with the data in order to have it like we want.
    for( int i( 0 ); i < nVoxels; ++i )
    {
        for( int j( 0 ); j < m_bands; ++j )
        {
            odfFloatData[i * m_bands + j] = shoreData[j * nVoxels + i]*radius;
        }
    }
    
    return odfFloatData;
}


// def EAPmatrix(self,r, theta, phi):
// 	"Compute the matrix function used to model the diffusion propagator"
// 	radialOrder = self.radialOrder
// 	zeta = 1/(4 * pi**2 * self.zeta)
// 
// 	if mod(radialOrder,2)==0:
// 		angularOrder=radialOrder
// 	else:
// 		angularOrder=radialOrder+1
// 	
// 	M = zeros((r.shape[0],(radialOrder+1)*((radialOrder+1)/2)*(2*radialOrder+1)))
// 	Y = SphericalHarmonics.matrix(theta, phi, order=angularOrder)
// 
// 	counter=0;
// 	for n in range(radialOrder+1):
// 		for l in range(0,n+1,2):
// 			for m in range(-l,l+1):
// 				j = SphericalHarmonics.j(l,m)
// 				#print(counter)
// 				#print "(n,l,m) = (%d,%d,%d)" % (n,l,m)
// 				#print(counter)
// 				M[:,counter] = (-1)**(n - l/2) * \
// 					Y[:,j] * \
// 					genlaguerre(n - l,l + 0.5)(r ** 2 / zeta) * \
// 					exp(- r ** 2 / (2 * zeta)) * \
// 					HermitePolynomial.kappa(zeta, n, l) * \
// 					(r ** 2 / zeta)**(l/2)
// 				#print(sum(genlaguerre(n - l/2,l + 0.5)(r ** 2 / zeta)))
// 				counter+=1
// 	return M[:,0:counter]


// double kappa(double zeta, int n, int l)
// {
// 	double res=0
// 	if (n-l<0)
// 	{
// 		res= sqrt()
// double Helper::getLegendrePlm( int i_m, double i_x )		
// Helper::getFactorial( i_l - l_absm )
// }
	
// @staticmethod
// def kappa(zeta, n , l):
// 	if n-l<0 :
// 		return sqrt( (2* 1) / (zeta ** 1.5 * sp.gamma(n + 1.5)) )
// 	else :
// 		return sqrt( (2* math.factorial(n -l)) / (zeta ** 1.5 * sp.gamma(n + 1.5)) )

void EAPs::createPropertiesSizer( PropertiesWindow *pParent )
{
    Glyph::createPropertiesSizer( pParent );
    
    setColorWithPosition( true );
    
    wxBoxSizer *pBoxMain = new wxBoxSizer( wxVERTICAL );
    
    //////////////////////////////////////////////////////////////////////////
    
    // TODO radius slider goes here
    //m_pSliderRadius = new MySlider( pParent, wxID_ANY, 5, 0, 10,    DEF_POS, wxSize( 100, -1 ), wxSL_HORIZONTAL | wxSL_AUTOTICKS );
    //m_pTxtThres    = new wxTextCtrl(   pParent, wxID_ANY, wxT( "0.5"), DEF_POS, wxSize(  40, -1 ), wxTE_READONLY);
    //m_pLblThres    = new wxStaticText( pParent, wxID_ANY, wxT( "Threshold" ) );
    //m_pBtnMainDir  = new wxButton(     pParent, wxID_ANY, wxT( "Recalculate" ), DEF_POS, wxSize( 140, -1 ) );
    //wxRadioButton *pRadDescoteauxBasis = new wxRadioButton( pParent, wxID_ANY, wxT( "Descoteaux" ), DEF_POS, DEF_SIZE, wxRB_GROUP );
    //wxRadioButton *pRadTournierBasis   = new wxRadioButton( pParent, wxID_ANY, wxT( "MRtrix" ) );
    //wxRadioButton *pRadDipyBasis   = new wxRadioButton( pParent, wxID_ANY, wxT( "Dipy" ) );
    ////     wxRadioButton *pRadOriginalBasis   = new wxRadioButton( pParent, wxID_ANY, wxT( "RR5768" ) );
    ////     wxRadioButton *pRadPTKBasis        = new wxRadioButton( pParent, wxID_ANY, wxT( "PTK" ) );
    
    //////////////////////////////////////////////////////////////////////////
    
    /*wxBoxSizer *pBoxFlood = new wxBoxSizer( wxHORIZONTAL );
    pBoxFlood->Add( m_pLblThres,   0, wxALIGN_CENTER_VERTICAL | wxALL, 1 );
    pBoxFlood->Add( m_pSliderFlood, 1, wxALIGN_CENTER_VERTICAL | wxALL, 1 );
    pBoxFlood->Add( m_pTxtThres, 0, wxFIXED_MINSIZE | wxALIGN_RIGHT | wxALIGN_CENTER_VERTICAL | wxALL, 1 );
    pBoxMain->Add( pBoxFlood, 0, wxEXPAND, 0 );*/
    
    //////////////////////////////////////////////////////////////////////////
    
    //pBoxMain->Add( m_pBtnMainDir, 0, wxALIGN_CENTER | wxEXPAND | wxRIGHT | wxLEFT, 24 );
    
    /*wxBoxSizer *pBoxShBasis = new wxBoxSizer( wxVERTICAL );
    pBoxShBasis->Add( new wxStaticText( pParent, wxID_ANY, wxT( "Sh Basis:" ) ), 0, wxALIGN_LEFT | wxALL, 1 );
    
    wxBoxSizer *pBoxShBasisRadios = new wxBoxSizer( wxVERTICAL );
    //     pBoxShBasisRadios->Add( pRadOriginalBasis,   0, wxALIGN_LEFT | wxALL, 1 );
    pBoxShBasisRadios->Add( pRadDescoteauxBasis, 0, wxALIGN_LEFT | wxALL, 1 );
    pBoxShBasisRadios->Add( pRadTournierBasis,   0, wxALIGN_LEFT | wxALL, 1 );
    pBoxShBasisRadios->Add( pRadDipyBasis,   0, wxALIGN_LEFT | wxALL, 1 );
    //     pBoxShBasisRadios->Add( pRadPTKBasis,        0, wxALIGN_LEFT | wxALL, 1 );
    pBoxShBasis->Add( pBoxShBasisRadios, 0, wxALIGN_LEFT | wxLEFT, 32 );
    
    pBoxMain->Add( pBoxShBasis, 0, wxFIXED_MINSIZE | wxEXPAND, 0 );*/
    
    //////////////////////////////////////////////////////////////////////////
    
    //pParent->Connect( m_pSliderFlood->GetId(),      wxEVT_COMMAND_SLIDER_UPDATED,       wxCommandEventHandler( PropertiesWindow::OnSliderAxisMoved ) );
    //pParent->Connect( m_pBtnMainDir->GetId(),       wxEVT_COMMAND_BUTTON_CLICKED,       wxCommandEventHandler( PropertiesWindow::OnRecalcMainDir ) );
    ////     pParent->Connect( pRadOriginalBasis->GetId(),   wxEVT_COMMAND_RADIOBUTTON_SELECTED, wxCommandEventHandler( PropertiesWindow::OnOriginalShBasis ) );
    //pParent->Connect( pRadDescoteauxBasis->GetId(), wxEVT_COMMAND_RADIOBUTTON_SELECTED, wxCommandEventHandler( PropertiesWindow::OnDescoteauxShBasis ) );
    //pParent->Connect( pRadTournierBasis->GetId(),   wxEVT_COMMAND_RADIOBUTTON_SELECTED, wxCommandEventHandler( PropertiesWindow::OnTournierShBasis ) );
    //pParent->Connect( pRadDipyBasis->GetId(),   wxEVT_COMMAND_RADIOBUTTON_SELECTED, wxCommandEventHandler( PropertiesWindow::OnDipyShBasis ) );
    ////     pParent->Connect( pRadPTKBasis->GetId(),        wxEVT_COMMAND_RADIOBUTTON_SELECTED, wxCommandEventHandler( PropertiesWindow::OnPTKShBasis ) );
    
    //////////////////////////////////////////////////////////////////////////
    
    //     pRadOriginalBasis->SetValue(   isShBasis( SH_BASIS_RR5768 ) );
    //pRadDescoteauxBasis->SetValue( isShBasis( SH_BASIS_DESCOTEAUX ) );
    //pRadTournierBasis->SetValue(   isShBasis( SH_BASIS_TOURNIER ) );
    //pRadDipyBasis->SetValue(   isShBasis( SH_BASIS_DIPY ) );
    //     pRadPTKBasis->SetValue(        isShBasis( SH_BASIS_PTK ) );
    
    m_pSliderLightAttenuation->SetValue( m_pSliderLightAttenuation->GetMin() );
    m_pSliderLightXPosition->SetValue( m_pSliderLightXPosition->GetMin() );
    m_pSliderLightYPosition->SetValue( m_pSliderLightYPosition->GetMin() );
    m_pSliderLightZPosition->SetValue( m_pSliderLightZPosition->GetMin() );
    
    //////////////////////////////////////////////////////////////////////////
    
    m_pPropertiesSizer->Add( pBoxMain, 0, wxFIXED_MINSIZE | wxEXPAND, 0 );
}

void EAPs::updatePropertiesSizer()
{
    //     Glyph::updatePropertiesSizer();
    DatasetInfo::updatePropertiesSizer();
    
    m_pSliderLightAttenuation->Enable( false );
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
    m_pToggleColorWithPosition->SetValue( getColorWithPosition() );
    
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
