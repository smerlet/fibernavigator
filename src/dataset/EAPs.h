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

#include <complex>
#include <map>

class MySlider;

class EAPs : public ODFs
{
	
	


public:
    // Constructor/Destructor
    EAPs( const wxString &filename );
    virtual ~EAPs();

    // From DatasetInfo
    bool load( nifti_image *pHeader, nifti_image *pBody );
    // TODO Sylvain
    // bool save( wxXmlNode *pNode ) const;
    
    // Reaction to events.
    bool updateDisplayRadius();
    
    // GUI Methods
    virtual void createPropertiesSizer( PropertiesWindow *parent );
    virtual void updatePropertiesSizer();    

    //Vars

protected:
    // Inherited from Glyph
    virtual bool createStructure( std::vector< float >& fileFloatData );
    virtual void drawGlyph      ( int      i_zVoxel, 
                                  int      i_yVoxel, 
                                  int      i_xVoxel, 
                                  AxisType i_axis );

private:
    // Variables
    float   m_displayRadius;
    
    // GUI elements
    MySlider *m_pSliderRadius;


	std::vector< float > shoreToSh(double radius);
	double shoreFunction(unsigned n, unsigned l, double zeta, double x)	;
	double kappa(unsigned n, unsigned l, double zeta);
	unsigned m_bands_EAP;
	unsigned radialOrder_EAP;
	unsigned angularOrder_EAP;
// 	float* eapData;
	std::vector< float > shoreDataAranged;

};

#endif /* EAPS_H */
