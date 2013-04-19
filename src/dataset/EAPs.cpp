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
/*

#include <math.h>
#include <gsl.h>*/
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
	
    //double radius=1e-3;
    
    // TODO Sylvain you should use the value that you want. The current value
    // is in m_displayRadius. You should set a default value in the constructor.
    double radius = 1.0;
	odfFloatData=shoreToSh( eapData, radius, nVoxels, m_bands);
	std::cout << "avant createStructure" << std::endl;
    // Once the file has been read successfully, we need to create the structure 
    // that will contain all the sphere points representing the ODFs.
	
// 	creer un pointeur sur un objet ODF dans le constructeur EAPs et ensuite s'en sevir pour appeler createStrucure
	// TODO Sylvain you are here.
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
    return ODFs::createStructure( shore_data );
}

void EAPs::drawGlyph(int zVoxel, int yVoxel, int xVoxel, AxisType axis)
{
    ODFs::drawGlyph(zVoxel, yVoxel, xVoxel, axis);
}

std::vector< float > EAPs::shoreToSh( float* shoreData,  double radius, int nVoxels, int m_bands)
{
// 	Converti les coefficient SHORE vers des coefficients SH

// see https://fr.wikipedia.org/wiki/C%2B%2B11
// dans <math.h>:
//  	double assoc_laguerre( unsigned n, unsigned m, double x ) ;
// 	double laguerre( unsigned n, double x ) ;
// 	exponentielles : exp
// (regarder aussi cmath.h)	
// 	genlaguerre(n - l,l + 0.5)(r ** 2 / zeta)
	unsigned n(3);
	unsigned l(2);
	double x(1e-2);
	double lagf(0.0);
	double zeta(700); //A inserer dans la fonction

	//Il me faut  definir quels sont les coefficient d'harmonic spherique pour un r fixe
	double res(shoreFunction(n,l,zeta,x));

	std::cout << "res: " << res << std::endl;
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


double EAPs::shoreFunction(unsigned n, unsigned l, double zeta, double x)
{
		double res(1);
		res*=boost::math::laguerre( n  , l , x );
		std::cout << "res: " << res << std::endl;
		res*=exp(-x/2);
		std::cout << "res: " << res << std::endl;
		res*=kappa(n,l,zeta);
		std::cout << "res: " << res << std::endl;
		res*=pow(x,2);

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


bool EAPs::updateDisplayRadius()
{
    float newRad = m_pSliderRadius->GetValue() / 100.0f;

    if( newRad != m_displayRadius )
    {
        m_displayRadius = newRad;
        
        // Sylvain TODO here update your coefficients and call createStructure.
        
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
