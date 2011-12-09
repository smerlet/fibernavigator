/*
 *  The Fibers class declaration.
 *
 */

#ifndef FIBERS_H_
#define FIBERS_H_

#include "wx/wxprec.h"

#ifndef WX_PRECOMP
#include "wx/wx.h"
#endif

#include <vector>

#include "DatasetInfo.h"
#include <GL/glew.h>
#include "KdTree.h"
#include "SplinePoint.h"
#include "../gui/SelectionObject.h"
#include "../misc/Fantom/FVector.h"
#include "Octree.h"

enum FiberFileType
{
    ASCII_FIBER = 0,
    ASCII_VTK   = 1,
    BINARY_VTK  = 2,
};

using namespace std;

/**
 * This class represents a set of fibers.
 * It supports loading different fibers file types.
 * It holds the fibers data, and can also be used to
 * compute statistics and selected fibers data.
 */
class Fibers : public DatasetInfo, public wxTreeItemData
{
public:
    Fibers( DatasetHelper *pDatasetHelper );
    // TODO copy constructors: should they be allowed?
    virtual ~Fibers();

    // Fibers loading methods
    bool    load(             wxString  filename );
    bool    loadTRK(    const wxString &filename );
    bool    loadCamino( const wxString &filename );
    bool    loadMRtrix( const wxString &filename );
    bool    loadPTK(    const wxString &filename );
    bool    loadVTK(    const wxString &filename );
    bool    loadDmri(   const wxString &filename );
    void    loadTestFibers();
    
    void    updateFibersColors();
    
    void    generateFiberVolume();
    
    void    save( wxString filename );
    void    saveDMRI( wxString filename );
    
    int     getPointsPerLine(     const int lineId );
    int     getStartIndexForLine( const int lineId );
    
    int     getLineForPoint( const int pointIdx );
    
    void    resetColorArray();
    
    void    updateLinesShown();

    void    generateKdTree();
    bool    getBarycenter( SplinePoint *pPoint );
    
    void    initializeBuffer();

    void    draw();    
    void    switchNormals( bool positive );

    float   getPointValue( int  ptIndex );
    int     getLineCount();
    int     getPointCount();
    bool    isSelected(    int  fiberId );
    
    void    setFibersLength();
    
    void    updateFibersFilters();

    vector< bool >  getFilteredFibers();

    void    flipAxis( AxisType i_axe );
    
    int     getFibersCount() const
    { 
        return m_countLines; 
    }
    
    // TODO check if we can set const
    Octree* getOctree() const
    {
        return m_pOctree;
    }
    
    vector< int > getReverseIdx() const
    {
        return m_reverse;
    }
    
    virtual void createPropertiesSizer( PropertiesWindow *pParent );
    virtual void updatePropertiesSizer();


    GLuint  getGLuint( )
    {
        return 0;
    };
    
    float   getMaxFibersLength()
    {
        return m_maxLength;
    }
    
    float   getMinFibersLength()
    {
        return m_minLength;
    }
    
    // Empty derived methods
    void    activateLIC()      {};
    void    clean()            {};
    
    void    drawVectors()      {};
    void    generateTexture()  {};
    void    generateGeometry() {};
    void    smooth()           {};

private:
    void            colorWithTorsion(     float *pColorData );
    void            colorWithCurvature(   float *pColorData );
    void            colorWithDistance(    float *pColorData );
    void            colorWithMinDistance( float *pColorData );
    
    string          intToString( const int number );
    void            toggleEndianess();
    
    void            calculateLinePointers();
    void            createColorArray( const bool colorsLoadedFromFile );
    
    void            resetLinesShown();
    vector< bool >  getLinesShown( SelectionObject *pSelectionObject );
    void            objectTest(    SelectionObject *pSelectionObject );
    
    void            barycenterTest( int left, int right, int axis );

    void            drawFakeTubes();
    void            drawSortedLines();

    void            freeArrays();

    bool            getFiberCoordValues( int fiberIndex, vector< Vector > &fiberPoints );

    // Variables
    bool            m_isSpecialFiberDisplay;
    Vector          m_barycenter;
    vector< float > m_boxMax;
    vector< float > m_boxMin;
    vector< float > m_colorArray;
    int             m_count;
    int             m_countLines;
    int             m_countPoints;
    bool            m_isInitialized;
    vector< int >   m_lineArray;
    vector< int >   m_linePointers;
    vector< float > m_pointArray;
    vector< float > m_normalArray;
    bool            m_normalsPositive;
    vector< int >   m_reverse;
    vector< bool >  m_selected;
    vector< bool >  m_filtered;
    vector< float > m_length;
    float           m_maxLength;
    float           m_minLength;
    vector<float  > m_localizedAlpha;
    float           m_cachedThreshold;

    KdTree          *m_pKdTree;
    Octree          *m_pOctree;
    
    // GUI members
    wxButton       *m_pGeneratesFibersDensityVolume;
    wxSlider       *m_pSliderFibersFilterMin;
    wxSlider       *m_pSliderFibersFilterMax;
    wxSlider       *m_pSliderFibersSampling;
    wxToggleButton *m_pToggleLocalColoring;
    wxToggleButton *m_pToggleNormalColoring;
    wxRadioButton  *m_pRadioNormalColoring;
    wxRadioButton  *m_pRadioDistanceAnchoring;
    wxRadioButton  *m_pRadioMinDistanceAnchoring;
    wxRadioButton  *m_pRadioCurvature;
    wxRadioButton  *m_pRadioTorsion;
};

#endif /* FIBERS_H_ */
