/////////////////////////////////////////////////////////////////////////////
// Name:            ODFs.h
// Author:          Imagicien ->LAMIRANDE-NADEAU Julien & NAZRATI R�da<-
// Creation Date:   10/07/2009
//
// Description: ODFs class.
//
// Last modifications:
//      by : Imagicien - 12/11/2009
/////////////////////////////////////////////////////////////////////////////
#ifndef ODFS_H_
#define ODFS_H_

#include "DatasetInfo.h"
#include "Glyph.h"

#include <complex>
#include <map>


using namespace std;

class ODFs : public Glyph
{
public:
    // Constructor/Destructor
    ODFs( DatasetHelper* datasetHelper );
    virtual ~ODFs();

    // From DatasetInfo
    void draw       ();
    bool load       ( wxString i_fileName );

    // Functions
    bool loadNifti  ( wxString i_fileName );
    virtual void createPropertiesSizer(PropertiesWindow *parent);
    virtual void updatePropertiesSizer();
    
    bool isShBasis( int i_sh_basis ) {return m_sh_basis == i_sh_basis;};
    vector< vector < float > > getCoeffs() {return m_coefficients;}
    vector< FMatrix > getShMatrix() {return m_shMatrix;}
    vector< FMatrix > getPhiTheta() {return m_phiThetaDirection;}

    void setShBasis(int value){m_sh_basis = value;}
    void changeShBasis(ODFs*,DatasetHelper*, int);
    
	void			  extractMaximas();

    //Vars
    wxString    m_lastODF_path;

    struct direction_value { double x,y,z,v; };
    struct direction { double x,y,z; };

	bool   m_isMaximasSet;
	float  m_axisThreshold;

    MySlider        *m_pSliderFlood;
    wxStaticText    *m_pTextThres;
    wxTextCtrl      *m_pTxtThresBox;
    wxButton        *m_pbtnMainDir;

private:
    // From Glyph
    bool createStructure  ( vector< float >& i_fileFloatData );
    void drawGlyph        ( int i_zVoxel, int i_yVoxel, int i_xVoxel, AxisType i_axis );
    void loadBuffer       ();
    void sliderPosChanged ( AxisType i_axis );  
    

    // Functions
    void             computeXRadiusSlice();
    void             computeYRadiusSlice();
    void             computeZRadiusSlice();

    void             computeRadiiArray          ( const FMatrix &i_B, vector< float > &i_C, 
                                                  vector< float >& o_radius, 
                                                  pair< float, float > &o_minMax );

    void             getODFpoints               ( FMatrix &i_phiThetaDirection, 
                                                  vector < float > &o_deformedMeshPts );    
    complex< float > getSphericalHarmonic       ( int i_l, int i_m, float i_theta, float i_phi );
    void             getSphericalHarmonicMatrix ( const vector< float > &i_meshPts, 
                                                  FMatrix &o_phiThetaDirection, 
                                                  FMatrix &o_shMatrix );
    void             getSphericalHarmonicMatrixRR5768 ( const vector< float > &i_meshPts, 
                                                           FMatrix &o_phiThetaDirection, 
                                                           FMatrix &o_shMatrix );
    void getSphericalHarmonicMatrixDescoteauxThesis ( const vector< float > &i_meshPts, 
                                                      FMatrix &o_phiThetaDirection, 
                                                      FMatrix &o_shMatrix );
    void       getSphericalHarmonicMatrixPTK ( const vector< float > &i_meshPts, 
                                               FMatrix &o_phiThetaDirection, 
                                               FMatrix &o_shMatrix );
    void             getSphericalHarmonicMatrixTournier ( const vector< float > &i_meshPts, 
                                                          FMatrix &o_phiThetaDirection, 
                                                          FMatrix &o_shMatrix );
    void             loadRadiusBuffer           ( AxisType i_axis );
    void             reloadRadiusBuffer         ( AxisType i_axis );

    
	
    std::vector<Vector> getODFmax(vector < float >  coefs,const FMatrix & SHmatrix, 
                                  const FMatrix & grad,
								  const float & max_thresh);
    void			    set_nbors(FMatrix i_phiThetaDirection);
    float               get_min_angle();
    void setScalingFactor( float i_scalingFactor );
    


    // Variables
    int     m_order;
    GLuint  m_radiusAttribLoc;
    GLuint* m_radiusBuffer; 
    wxRadioButton *m_pRadiobtnOriginalBasis;
    wxRadioButton *m_pRadiobtnDescoteauxBasis;
    wxRadioButton *m_pRadiobtnTournierBasis;
    wxRadioButton *m_pRadiobtnPTKBasis;
    


    vector< vector < float > >          m_coefficients;
    vector< vector < float > >          m_radius;
    vector< FMatrix >                   m_shMatrix;
    vector< FMatrix >                   m_phiThetaDirection;    
    vector< float >                     m_meshPts;    
    map< int, pair< float, float > >    m_radiiMinMaxMap;
  
    float								m_angle_min;
    std::vector<std::pair<float,int> >* m_nbors;
    std::vector<std::vector<Vector> >   m_mainDirections;

	int                                 m_sh_basis;
};

#endif /* ODFS_H_ */
