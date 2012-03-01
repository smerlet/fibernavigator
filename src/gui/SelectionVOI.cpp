// TODO header syntax

#include "SelectionVOI.h"

#include "../misc/IsoSurface/CBoolIsoSurface.h"

#include <algorithm>
#include <functional>
#include <limits>

///////////////////////////////////////////////////////////////////////////
// Constructor
// i_center             : The center of the ellipsoid.
// i_size               : The size of the ellipsoid.
// i_dataHelper         : The datasetHelper associated with this ellipsoid.
///////////////////////////////////////////////////////////////////////////
/*SelectionEllipsoid::SelectionEllipsoid( Vector i_center, Vector i_size, DatasetHelper* i_datasetHelper ) :
SelectionObject( i_center, i_size, i_datasetHelper )
{
    m_gfxDirty   = true;
    m_name       = wxT( "ellipsoid" );
    m_objectType = ELLIPSOID_TYPE;
    
    update();
}*/

///////////////////////////////////////////////////////////////////////////
// Destructor
///////////////////////////////////////////////////////////////////////////
/*SelectionEllipsoid::~SelectionEllipsoid()
{
    
}*/

SelectionVOI::SelectionVOI( DatasetHelper *pDH, Anatomy *pSourceAnatomy, const float threshold, const ThresholdingOperationType opType )
    : SelectionObject( Vector(0.0, 0.0, 0.0), Vector(0.0, 0.0, 0.0), pDH ),
      m_voiSize( 0 )
{
    m_nbRows   = pSourceAnatomy->getRows();
    m_nbCols   = pSourceAnatomy->getColumns();
    m_nbFrames = pSourceAnatomy->getFrames();
    
    vector< float > * const pAnatDataset = pSourceAnatomy->getFloatDataset();
    
    m_includedVoxels.assign( pAnatDataset->size(), false );
    
    if( opType == THRESHOLD_EQUAL )
    {
        std::transform( pAnatDataset->begin(), pAnatDataset->end(),
                       m_includedVoxels.begin(), bind2nd( std::equal_to< float >(), threshold ) );
    }
    else if( opType == THRESHOLD_GREATER )
    {
        std::transform( pAnatDataset->begin(), pAnatDataset->end(),
                       m_includedVoxels.begin(), bind2nd( std::greater< float >(), threshold ) );
    }
    else if( opType == THRESHOLD_GREATER_EQUAL )
    {
        std::transform( pAnatDataset->begin(), pAnatDataset->end(),
                       m_includedVoxels.begin(), bind2nd( std::greater_equal< float >(), threshold ) );
    }
    else if( opType == THRESHOLD_SMALLER )
    {
        std::transform( pAnatDataset->begin(), pAnatDataset->end(),
                       m_includedVoxels.begin(), bind2nd( std::less< float >(), threshold ) );
    }
    else if( opType == THRESHOLD_SMALLER_EQUAL )
    {
        std::transform( pAnatDataset->begin(), pAnatDataset->end(),
                       m_includedVoxels.begin(), bind2nd( std::less_equal< float >(), threshold ) );

    }

    //m_isosurface    = new CIsoSurface( m_datasetHelper, pSourceAnatomy );
    m_pIsoSurface = new CBoolIsoSurface( m_datasetHelper, m_includedVoxels );
    //m_isosurface->setThreshold( threshold );
    //m_isosurface->GenerateWithThreshold();
    m_pIsoSurface->GenerateSurface();
    
    wxString mystring(wxT("[VOI] - ") + pSourceAnatomy->getName());
    m_name          =  mystring;
    
    /*m_objectType    = CISO_SURFACE_TYPE;
    m_sourceAnatomy = i_anatomy;*/
    // TODO HERE was going to compute the value of the 
    // m_center and m_size params.
    unsigned int xIdxMin( m_nbCols );
    unsigned int yIdxMin( m_nbRows );
    unsigned int zIdxMin( m_nbFrames );
    unsigned int xIdxMax( 0 );
    unsigned int yIdxMax( 0 );
    unsigned int zIdxMax( 0 );
    
    unsigned int dataIdx( 0 );
    for( unsigned int zPos( 0 ); zPos < m_nbFrames; ++zPos )
    {
        for( unsigned int yPos( 0 ); yPos < m_nbRows; ++yPos )
        {
            for( unsigned int xPos( 0 ); xPos < m_nbCols; ++xPos, ++dataIdx )
            {
                if( m_includedVoxels[dataIdx] )
                {
                    xIdxMin = std::min( xIdxMin, xPos );
                    yIdxMin = std::min( yIdxMin, yPos );
                    zIdxMin = std::min( zIdxMin, zPos );
                    xIdxMax = std::max( xIdxMax, xPos );
                    yIdxMax = std::max( yIdxMax, yPos );
                    zIdxMax = std::max( zIdxMax, zPos );
                }
            }
        }
    }
    
    // Convert to space coordinates
    float spaceXMin( xIdxMin * pDH->m_xVoxel );
    float spaceYMin( yIdxMin * pDH->m_yVoxel );
    float spaceZMin( zIdxMin * pDH->m_zVoxel );
    
    // For the max value, since we want to grab all of the voxel,
    // and the idx * the size of the voxel gives the begininng of the voxel
    // in space, we adjust with the +1, then the - 0.1 * voxelSize to
    // make sure not to fall in an inexisting coordinate.
    float spaceXMax( ( xIdxMax + 1 ) * pDH->m_xVoxel - 0.1 * pDH->m_xVoxel );
    float spaceYMax( ( yIdxMax + 1 ) * pDH->m_yVoxel - 0.1 * pDH->m_yVoxel );
    float spaceZMax( ( zIdxMax + 1 ) * pDH->m_zVoxel - 0.1 * pDH->m_zVoxel );
    //float spaceXMax( ( xIdxMax  ) * pDH->m_xVoxel  );
    //float spaceYMax( ( yIdxMax  ) * pDH->m_yVoxel  );
    //float spaceZMax( ( zIdxMax  ) * pDH->m_zVoxel  );
    
    setCenter( ( spaceXMax + spaceXMin ) / 2.0f, 
               ( spaceYMax + spaceYMin ) / 2.0f, 
               ( spaceZMax + spaceZMin ) / 2.0f );
    
    setSize( spaceXMax - spaceXMin, spaceYMax - spaceYMin, spaceZMax - spaceZMin );
    
    m_voiSize = std::count( m_includedVoxels.begin(), m_includedVoxels.end(), true );
}

///////////////////////////////////////////////////////////////////////////
// This is the specific implementation to draw a selectionEllipsoid object.
//
// i_color      : The color of the ellipsoid to draw.
///////////////////////////////////////////////////////////////////////////
void SelectionVOI::drawObject( GLfloat * pColor )
{
    // Do not use the suggested color.
    // For the moment a selection VOI will always be red.
    glColor4f( 1.0f, 0.0f, 0.0f, pColor[3] );
    
    //m_isosurface->draw();
    m_pIsoSurface->draw();
    /*glColor4f( i_color[0], i_color[1], i_color[2], i_color[3] );
    
    glDepthMask(GL_FALSE);
    glEnable( GL_BLEND );
    glBlendFunc( GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA );    
    
    glPushMatrix();
    glTranslatef( m_center.x, m_center.y, m_center.z );
    GLUquadricObj* l_quadric = gluNewQuadric();
    gluQuadricNormals( l_quadric, GLU_SMOOTH );
    glScalef( m_xRadius, m_yRadius, m_zRadius );
    gluSphere( l_quadric, 1.0f, 32, 32 );
    glPopMatrix();
    
    glDisable( GL_BLEND );
    glDepthMask(GL_TRUE);*/
}

///////////////////////////////////////////////////////////////////////////
// COMMENT
//
// i_ray        :
///////////////////////////////////////////////////////////////////////////
hitResult SelectionVOI::hitTest( Ray* i_ray )
{    
    hitResult hr = { false, 0.0f, 0, NULL };
    
    /*if( m_isVisible && m_isActive && m_objectType == ELLIPSOID_TYPE ) 
    {
        int   picked  = 0;
        float tpicked = 0;
        float cx = m_center.x;
        float cy = m_center.y;
        float cz = m_center.z;
        float sx = m_size.x * m_datasetHelper->m_xVoxel;
        float sy = m_size.y * m_datasetHelper->m_yVoxel;
        float sz = m_size.z * m_datasetHelper->m_zVoxel;
        
        if( wxGetKeyState( WXK_CONTROL ) )
        {
            BoundingBox *bb = new BoundingBox( cx, cy, cz, sx, sy, sz );
            
            bb->setCenter( m_minX , cy, cz );
            bb->setSize( sx, sy, sz );
            bb->setSizeX( m_datasetHelper->m_xVoxel );
            hr = bb->hitTest( i_ray );
            if( hr.hit )
            {
                
                if( picked == 0 )
                {
                    picked = 11;
                    tpicked = hr.tmin;
                }
                else 
                {
                    if( hr.tmin < tpicked )
                    {
                        picked = 11;
                        tpicked = hr.tmin;
                    }
                }
            }
            bb->setCenter( m_maxX, cy, cz );
            hr = bb->hitTest( i_ray );
            if( hr.hit ) 
            {
                if( picked == 0 ) 
                {
                    picked = 12;
                    tpicked = hr.tmin;
                }
                else 
                {
                    if( hr.tmin < tpicked ) 
                    {
                        picked  = 12;
                        tpicked = hr.tmin;
                    }
                }
            }
            bb->setCenter( cx, m_minY, cz );
            bb->setSize( sx, sy, sz);
            bb->setSizeY( m_datasetHelper->m_yVoxel );
            hr = bb->hitTest( i_ray );
            if( hr.hit )
            {
                if( picked == 0 )
                {
                    picked = 13;
                    tpicked = hr.tmin;
                }
                else 
                {
                    if( hr.tmin < tpicked ) 
                    {
                        picked = 13;
                        tpicked = hr.tmin;
                    }
                }
            }
            bb->setCenter( cx, m_maxY, cz );
            hr = bb->hitTest( i_ray );
            if( hr.hit)
            {
                if( picked == 0 )
                {
                    picked = 14;
                    tpicked = hr.tmin;
                }
                else 
                {
                    if( hr.tmin < tpicked ) 
                    {
                        picked = 14;
                        tpicked = hr.tmin;
                    }
                }
            }
            bb->setCenter( cx, cy, m_minZ );
            bb->setSize( sx, sy, sz );
            bb->setSizeZ( m_datasetHelper->m_zVoxel );
            hr = bb->hitTest( i_ray );
            if( hr.hit ) 
            {
                if( picked == 0 )
                {
                    picked  = 15;
                    tpicked = hr.tmin;
                }
                else 
                {
                    if( hr.tmin < tpicked )
                    {
                        picked = 15;
                        tpicked = hr.tmin;
                    }
                }
            }
            bb->setCenter( cx, cy, m_maxZ );
            hr = bb->hitTest( i_ray );
            if( hr.hit )
            {
                if( picked == 0 ) 
                {
                    picked = 16;
                    tpicked = hr.tmin;
                }
                else 
                {
                    if( hr.tmin < tpicked ) 
                    {
                        picked = 16;
                        tpicked = hr.tmin;
                    }
                }
            }
            
        }
        else  // if (wxGetKeyState(WXK_CONTROL))
        {
            BoundingBox *bb = new BoundingBox( cx, cy, cz, sx, sy, sz );
            hr = bb->hitTest( i_ray );
            if( hr.hit )
            {
                if( picked == 0 )
                {
                    picked = 10;
                    tpicked = hr.tmin;
                    
                }
                else 
                {
                    if( hr.tmin < tpicked )
                    {
                        picked = 10;
                        tpicked = hr.tmin;
                    }
                }
            }
        }
        if( picked != 0 )
        {
            hr.hit = true;
            hr.tmin = tpicked;
            hr.picked = picked;
            hr.object = this;
        }
        
    }*/
    
    m_hitResult = hr;
    
    return hr;
}

///////////////////////////////////////////////////////////////////////////
// This function set the proper radius of the ellipsoid.
///////////////////////////////////////////////////////////////////////////
void SelectionVOI::objectUpdate()
{
    /*m_xRadius = ( m_maxX - m_minX ) / 2.0f;
    m_yRadius = ( m_maxY - m_minY ) / 2.0f;
    m_zRadius = ( m_maxZ - m_minZ ) / 2.0f;*/
}

bool SelectionVOI::isPointInside( const float xPos, const float yPos, const float zPos ) const
{
    unsigned int xVoxelCoord( static_cast< unsigned int >( ( xPos / m_datasetHelper->m_xVoxel ) ) );
    unsigned int yVoxelCoord( static_cast< unsigned int >( ( yPos / m_datasetHelper->m_yVoxel ) ) );
    unsigned int zVoxelCoord( static_cast< unsigned int >( ( zPos / m_datasetHelper->m_zVoxel )  ) );
    
    unsigned int dataCoord( zVoxelCoord * m_nbCols * m_nbRows + yVoxelCoord * m_nbCols + xVoxelCoord );
    
    return m_includedVoxels.at( dataCoord );
}

void SelectionVOI::createPropertiesSizer( PropertiesWindow *pParent )
{
    SelectionObject::createPropertiesSizer( pParent );
    
    m_propertiesSizer->AddSpacer( 8 );
    
    wxSizer *pSizer = new wxBoxSizer( wxHORIZONTAL );
    
    m_pVOISize = new wxTextCtrl( pParent, wxID_ANY, wxString::Format( wxT("%d"), m_voiSize ), wxDefaultPosition, wxDefaultSize , wxTE_CENTRE | wxTE_READONLY);    
    m_pVOISize->SetBackgroundColour( *wxLIGHT_GREY );
    
    pSizer->Add( new wxStaticText( pParent, wxID_ANY, wxT( "Nb. of voxels: " ),wxDefaultPosition, wxDefaultSize, wxALIGN_CENTER ), 0, wxALIGN_CENTER );
    pSizer->Add( m_pVOISize, 0, wxALIGN_CENTER );
    
    m_propertiesSizer->Add( pSizer, 0, wxALIGN_CENTER );
}

void SelectionVOI::updatePropertiesSizer()
{
    SelectionObject::updatePropertiesSizer();
}

SelectionVOI::~SelectionVOI()
{
    delete m_pIsoSurface;
    m_pIsoSurface = NULL;
}