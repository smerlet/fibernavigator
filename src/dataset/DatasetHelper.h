/*
 * DatasetHelper.h
 *
 *  Created on: 27.07.2008
 *      Author: ralph
 */
#ifndef DATASETHELPER_H_
#define DATASETHELPER_H_

#include "wx/wxprec.h"

#ifndef WX_PRECOMP
#include "wx/wx.h"
#endif

#include "wx/textfile.h"
#include "wx/file.h"
#include "wx/image.h"
#include "wx/wfstream.h"
#include "wx/datstrm.h"
#include "wx/txtstrm.h"
#include "wx/xml/xml.h"

#include <vector>

#include "DatasetInfo.h"

#include "../gui/SelectionObject.h"

#include "../gui/MainFrame.h"
#include "../gfx/TheScene.h"
#include "SplinePoint.h"

#include "AnatomyHelper.h"
#include "../gfx/ShaderHelper.h"

#include "../misc/lic/TensorField.h"

#include "../misc/Fantom/FMatrix.h"

#include "../Logger.h"

class MainFrame;
class DatasetInfo;
class TheScene;
class AnatomyHelper;
class ShaderHelper;
class SplinePoint;
class SelectionObject;
class SelectionTree;
class Fibers;
class FibersGroup;
class TensorField;
class Surface;

// TODO to be removed
typedef std::vector< std::vector< SelectionObject* > > SelectionObjectList;
typedef std::vector< SelectionObject* > SelectionObjectVector;

class DatasetHelper 
{
public:
    // Constructor/destructor
    DatasetHelper( MainFrame* i_mainFrame );
    virtual ~DatasetHelper();

    // Functions
    //void out_of_memory(); 
    bool load( const int i_index );
    bool load( wxString    i_filename,
               int         i_index		= -1, 
               const float i_threshold  = 0.0f, 
               const bool  i_active		= true,
               const bool  i_showFS		= true, 
               const bool  i_useTex		= true, 
               const float i_alpha		= 1.0f,
			   wxString    i_name		= _T( ""),
			   int		   i_version	= 1,
			   const bool  i_isFiberGroup = false,
			   const bool  i_isScene = false );
    void finishLoading ( DatasetInfo*, bool isChild = false );

    bool loadFmriClusters();
    bool loadScene     ( const wxString i_filename );
    bool loadTextFile  ( wxString* i_string, const wxString i_filename );
    bool fileNameExists( const wxString i_filename );

    //! Saves the current scene to an xml file
    void save( const wxString filename );

    SelectionObjectList getSelectionObjects();
    void   deleteAllPoints();
    void   deleteAllSelectionObjects();
    void   updateAllSelectionObjects();
    Vector mapMouse2World( const int i_x, const int i_y,GLdouble i_projection[16], GLint i_viewport[4], GLdouble i_modelview[16]);
    Vector mapMouse2WorldBack( const int i_x, const int i_y,GLdouble i_projection[16], GLint i_viewport[4], GLdouble i_modelview[16]);   

    void createIsoSurface();
    void createDistanceMapAndIso();
    void createDistanceMap();
    void createCutDataset();
    /*
     * Called from MainFrame when a kdTree thread signals it's finished
     */
    void treeFinished();

    void updateView( const float i_x, const float i_y, const float i_z );

    void changeZoom( const int i_z );
    void moveScene ( int i_x, int i_y );

    void doMatrixManipulation();

	bool getFibersGroupDataset( FibersGroup* &i_fiberGroup );
	bool getSelectedFiberDataset ( Fibers* &i_fiber );
    bool getSurfaceDataset( Surface* &i_surface );
    bool getTextureDataset( std::vector< DatasetInfo* > &o_types ); 
    std::vector< float >* getVectorDataset();
    TensorField* getTensorField();

    void toggleShowAllSelectionObjects() { m_showObjects = ! m_showObjects; };
    void toggleActivateAllSelectionObjects() { m_activateObjects = ! m_activateObjects; };
    bool togglePointMode()        { return m_pointMode = ! m_pointMode; };
    bool getPointMode()           { return m_pointMode; };

    void updateLoadStatus();
	void updateItemsId();
	void updateItemsPosition();

    void doLicMovie       ( int i_mode );
    void createLicSliceSag( int i_slize );
    void createLicSliceCor( int i_slize );
    void createLicSliceAxi( int i_slize );
    void licMovieHelper();

    void increaseAnimationStep();
    int  m_debugLevel;

    /////////////////////////////////////////////////////////////////////////////////
    // general info about the datasets
    /////////////////////////////////////////////////////////////////////////////////
    int                   m_rows;
    int                   m_columns;
    int                   m_frames;
    std::vector< float >* m_floatDataset;
    std::vector<Vector>   m_rulerPts;
    bool                  m_isRulerToolActive;
    bool                  m_isDrawerToolActive;
    double                m_rulerFullLength;
    double                m_rulerPartialLength;
    int                   m_fibersSamplingFactor;
    bool                  m_isSegmentActive;
    int                   m_SegmentMethod;
    bool                  m_isFloodfillActive;
    bool                  m_isSelectBckActive;
    bool                  m_isSelectObjActive;
    bool                  m_isObjfilled;
    bool                  m_isBckfilled;
    bool                  m_isObjCreated;
    bool                  m_isBckCreated;
    bool                  m_isBoxCreated;
    bool                  m_thresSliderMoved;
    bool graphcutReady()  { return (m_isObjfilled && m_isBckfilled && m_isBoxCreated); };
    bool m_showCrosshair;
        
    /////////////////////////////////////////////////////////////////////////////////
    //RTT vars
    /////////////////////////////////////////////////////////////////////////////////
    bool m_isRTTReady;
    bool m_isRTTDirty;
    
	bool m_isRTTActive;
    bool m_isRandomSeeds;
	bool m_interpolateTensors;
	bool m_isFileSelected;


    float m_xVoxel;
    float m_yVoxel;
    float m_zVoxel;
    
    FMatrix m_niftiTransform;

    unsigned int m_countFibers;

    bool m_loadDataset;
    bool m_scnFileLoaded;
    bool m_anatomyLoaded;
    bool m_meshLoaded;
    bool m_fibersGroupLoaded;
	bool m_fibersLoaded;
    bool m_vectorsLoaded;
    bool m_tensorsFieldLoaded;
    bool m_tensorsLoaded;
    bool m_ODFsLoaded;
    bool m_surfaceLoaded;
    bool m_surfaceIsDirty;

    /////////////////////////////////////////////////////////////////////////////////
    // state variables for rendering
    /////////////////////////////////////////////////////////////////////////////////
    Matrix4fT m_transform;
    bool      m_useVBO;
    int       m_quadrant;
    int       m_textures;
    //! if set the shaders will be reloaded during next render() call
    bool      m_scheduledReloadShaders;
    // the screenshot button has been pressed, next render pass it will be executed
    bool      m_scheduledScreenshot;

    bool      m_showObjects;
    bool      m_activateObjects;
    bool      m_blendAlpha;
    bool      m_pointMode;
    bool      m_isShowAxes;

    int       m_animationStep;

    

    float     m_frustum[6][4]; // Contains the information of the planes forming the frustum.
    /////////////////////////////////////////////////////////////////////////////////
    // state variables for menu entries
    /////////////////////////////////////////////////////////////////////////////////
    bool m_showSagittal;
    bool m_showCoronal;
    bool m_showAxial;

    
    wxImage m_drawColorIcon;
    float m_xSlize;
    float m_ySlize;
    float m_zSlize;

    bool  m_lighting;         // m_lighting for fibers.
    bool  m_blendTexOnMesh;   // Ignore threshold for textures on meshes.
    bool  m_useLic;           // Show the lic texture on spline surface.
    bool  m_drawVectors;      // Draw vectors as small lines on spline surface.
    float m_normalDirection;  // Normal direction of the spline surface.
    bool  m_geometryShadersSupported;
    bool  m_clearToBlack;
    bool  m_useFibersGeometryShader;
    bool  m_filterIsoSurf;
    int   m_colorMap;
    bool  m_showColorMapLegend;
    bool  m_displayMinMaxCrossSection;
    bool  m_displayGlyphOptions;

	
	enum  DrawMode
	{
		DRAWMODE_PEN = 0,
		DRAWMODE_ERASER = 1,
        DRAWMODE_INVALID
	};
	DrawMode m_drawMode;
	int     m_drawSize;
	bool    m_drawRound;
	bool    m_draw3d;
    bool    m_canUseColorPicker;
	wxColor m_drawColor;
	

    bool  m_morphing;

    bool  m_boxLockIsOn;
    bool  m_semaphore;
    int   m_threadsActive;

   /////////////////////////////////////////////////////////////////////////////////
    // variables for mouse interaction
    /////////////////////////////////////////////////////////////////////////////////
    bool  m_isDragging;
    bool  m_isrDragging;
    bool  m_ismDragging;
    float m_zoom;
    float m_xMove;
    float m_yMove;

    bool  m_texAssigned;
    bool  m_selBoxChanged;
    bool  m_guiBlocked;

    int   m_geforceLevel;

    wxString m_lastError;
    wxString m_lastPath;
    wxString m_scenePath;
    wxString m_scnFileName;
    wxString m_screenshotPath;
    wxString m_screenshotName;
    
    // TODO place in the correct place
    SelectionTree *m_pSelectionTree;

    /////////////////////////////////////////////////////////////////////////////////
    // pointers to often used objects
    /////////////////////////////////////////////////////////////////////////////////
    AnatomyHelper*   m_anatomyHelper;
    SelectionObject* m_boxAtCrosshair;
    SplinePoint*     m_lastSelectedPoint;
    SelectionObject* m_lastSelectedObject;
    TheScene*        m_theScene;
    MainFrame*       m_mainFrame;
    ShaderHelper*    m_shaderHelper;
    
    
};

#define ID_KDTREE_FINISHED    50

#endif /* DATASETHELPER_H_ */
