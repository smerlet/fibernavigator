/////////////////////////////////////////////////////////////////////////////
// Name:            EAPs.cpp
// Author:          athena team-project -> MERLET Sylvain <-
// Creation Date:   20/04/2013
//
// Description: EAPs class.
//
// Last modifications:
//      by : SMerlet - 21/04/2013
/////////////////////////////////////////////////////////////////////////////
#ifndef EAPS_H_ 
#define EAPS_H_

#include "DatasetInfo.h"
#include "Glyph.h"
#include "ODFs.h"
#include "../misc/nifti/nifti1_io.h"
#include "../misc/Fantom/FMatrix.h"

#include <complex>
#include <map>

class MySlider;

// enum SH_BASIS { SH_BASIS_RR5768, SH_BASIS_DESCOTEAUX, SH_BASIS_TOURNIER, SH_BASIS_PTK, SH_BASIS_DIPY };

class EAPs : public Glyph
{
	
	


public:
    // Constructor/Destructor
    EAPs( const wxString &filename );
    virtual ~EAPs();

    // From DatasetInfo
    bool load( nifti_image *pHeader, nifti_image *pBody );
 

    //Vars
    bool   m_isMaximasSet;
    float  m_axisThreshold;
	ODFs * odfs; // Ne pas oublier de delete l'odfs dans le destructeur



private:
    // From Glyph
//     bool createStructure  ( std::vector< float > &i_fileFloatData );


    // Variables
    int     m_order;
    GLuint  m_radiusAttribLoc;
    GLuint* m_radiusBuffer; 

    std::vector<std::pair<float,int> >* m_nbors;
	std::vector< float > shoreToSh( float* shoreData, double radius, int nVoxels, int m_bands);
	
	
    SH_BASIS                            m_sh_basis;
};

#endif /* EAPS_H */
