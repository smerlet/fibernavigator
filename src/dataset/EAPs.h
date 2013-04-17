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
    // TODO implement
    //void draw();
    bool load( nifti_image *pHeader, nifti_image *pBody );
    // TODO Sylvain
    // bool save( wxXmlNode *pNode ) const;
    
    // Methods
    virtual void createPropertiesSizer( PropertiesWindow *parent );
    virtual void updatePropertiesSizer();
    
    // Inherited from Glyph
    

    //Vars
    // TODO those should all be private.
    bool   m_isMaximasSet;
    float  m_axisThreshold;
	ODFs * odfs; // Ne pas oublier de delete l'odfs dans le destructeur

protected:
    // Inherited from Glyph
    virtual bool createStructure( std::vector< float >& fileFloatData );
    virtual void drawGlyph      ( int      i_zVoxel, 
                                  int      i_yVoxel, 
                                  int      i_xVoxel, 
                                  AxisType i_axis );


private:
    // Variables
    int     m_order;
    float   m_displayRadius;
    GLuint  m_radiusAttribLoc;
    GLuint* m_radiusBuffer;
    
    // GUI elements
    MySlider *m_pSliderRadius;

    std::vector<std::pair<float,int> >* m_nbors;
	std::vector< float > shoreToSh( float* shoreData, double radius, int nVoxels, int m_bands);
	
	
    SH_BASIS                            m_sh_basis;
};

#endif /* EAPS_H */
