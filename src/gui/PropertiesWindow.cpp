#include <wx/colordlg.h>

#include "PropertiesWindow.h"
#include "SelectionBox.h"
#include "SelectionEllipsoid.h"
#include "SelectionTree.h"
#include "TrackingWindow.h"
#include "../dataset/Anatomy.h"
#include "../dataset/Fibers.h"
#include "../dataset/FibersGroup.h"
#include "../dataset/ODFs.h"
#include "../dataset/SplinePoint.h"
#include "../dataset/Surface.h"
#include "../dataset/Tensors.h"
#include "../gui/SelectionObject.h"
#include "../gui/SelectionVOI.h"
#include "../misc/IsoSurface/CIsoSurface.h"
#include "../Logger.h"

IMPLEMENT_DYNAMIC_CLASS(PropertiesWindow, wxScrolledWindow)

BEGIN_EVENT_TABLE(PropertiesWindow, wxScrolledWindow)
EVT_PAINT(      PropertiesWindow::OnPaint)
EVT_SIZE(       PropertiesWindow::OnSize)
END_EVENT_TABLE()

PropertiesWindow::PropertiesWindow( wxWindow *parent, MainFrame *mf, wxWindowID id,
                    const wxPoint &pos, const wxSize &size )
        : wxScrolledWindow( parent, id, pos, size, wxBORDER_NONE, _T("test canvas") )
{
    m_noteBook = parent;
    m_mainFrame = mf;
    SetBackgroundColour( *wxLIGHT_GREY );
    SetCursor( wxCursor( wxCURSOR_HAND ) );
    propertiesSizer = new wxBoxSizer( wxVERTICAL );
    SetSizer( propertiesSizer );
    SetAutoLayout(true);
}


void PropertiesWindow::OnSize( wxSizeEvent &WXUNUSED(event) )
{

}

void PropertiesWindow::OnPaint( wxPaintEvent &WXUNUSED(event) )
{
    wxPaintDC dc( this );
}

void PropertiesWindow::OnListItemUp(wxCommandEvent& WXUNUSED(event))
{
	DatasetInfo* pDatasetInfo = ((DatasetInfo*)m_mainFrame->m_pListCtrl->GetItemData(m_mainFrame->m_currentListItem));
	if(pDatasetInfo != NULL)
	{
		if( pDatasetInfo->getType() == FIBERSGROUP )
		{
			FibersGroup* pFibersGroup = (FibersGroup*)pDatasetInfo;
			pFibersGroup->OnMoveUp();
			m_mainFrame->refreshAllGLWidgets();
			return;
		}
	}
	
	long prevItemId = m_mainFrame->m_currentListItem - 1;
	if(prevItemId > -1)
	{
		DatasetInfo* pPrevDatasetInfo = ((DatasetInfo*)m_mainFrame->m_pListCtrl->GetItemData(prevItemId));
		if( pPrevDatasetInfo != NULL)
		{
			if( pPrevDatasetInfo->getType() == FIBERS && pDatasetInfo->getType() != FIBERS )
			{
				FibersGroup* pFibersGroup;
				m_mainFrame->m_pDatasetHelper->getFibersGroupDataset(pFibersGroup);
				if(pFibersGroup != NULL)
				{
					int nbChilds = pFibersGroup->getFibersCount();
					m_mainFrame->m_pListCtrl->moveItemAt(m_mainFrame->m_currentListItem, prevItemId - nbChilds);
					m_mainFrame->m_pListCtrl->EnsureVisible(prevItemId - nbChilds);
					m_mainFrame->refreshAllGLWidgets();
					return;
				}
			}
		}
	}
    m_mainFrame->m_pListCtrl->moveItemUp(m_mainFrame->m_currentListItem);
    m_mainFrame->m_pListCtrl->EnsureVisible(m_mainFrame->m_currentListItem);   
    m_mainFrame->refreshAllGLWidgets();
}

void PropertiesWindow::OnListItemDown( wxCommandEvent& WXUNUSED(event) )
{
	DatasetInfo* pDatasetInfo = ((DatasetInfo*)m_mainFrame->m_pListCtrl->GetItemData(m_mainFrame->m_currentListItem));
	if(pDatasetInfo != NULL)
	{
		if( pDatasetInfo->getType() == FIBERSGROUP )
		{
			DatasetInfo* pDataset = (DatasetInfo*)m_mainFrame->m_pListCtrl->GetItemData(m_mainFrame->m_pListCtrl->GetItemCount() - 1);
			if(pDataset->getType() != FIBERS ) // if the last item in the list is not a fiber (child), then move down
			{
				FibersGroup* pFibersGroup = (FibersGroup*)pDatasetInfo;
				pFibersGroup->OnMoveDown();
				m_mainFrame->refreshAllGLWidgets();
			}
			return;
		}
	}
	
	long nextItemId = m_mainFrame->m_currentListItem + 1;
	if(nextItemId > -1)
	{
		DatasetInfo* pNextDatasetInfo = ((DatasetInfo*)m_mainFrame->m_pListCtrl->GetItemData(nextItemId));
		if( pNextDatasetInfo != NULL)
		{
			if( pNextDatasetInfo->getType() == FIBERSGROUP )
			{
				FibersGroup* pFibersGroup = (FibersGroup*)pNextDatasetInfo;
				int nbChilds = pFibersGroup->getFibersCount();
				m_mainFrame->m_pListCtrl->moveItemAt(m_mainFrame->m_currentListItem, nextItemId + nbChilds);
				m_mainFrame->m_pListCtrl->EnsureVisible(nextItemId + nbChilds);
				m_mainFrame->refreshAllGLWidgets();
				return;
			}
		}
	}
	m_mainFrame->m_pListCtrl->moveItemDown(m_mainFrame->m_currentListItem);
    m_mainFrame->m_pListCtrl->EnsureVisible(m_mainFrame->m_currentListItem);
    m_mainFrame->refreshAllGLWidgets();
}

void PropertiesWindow::OnDeleteListItem( wxEvent& WXUNUSED(event) )
{
    m_mainFrame->deleteListItem();
}

void PropertiesWindow::OnToggleIntensityBtn( wxEvent& WXUNUSED(event) )
{
	if (m_mainFrame->m_pCurrentSceneObject != NULL && m_mainFrame->m_currentListItem != -1)
    {
		DatasetInfo* pDatasetInfo = ((DatasetInfo*)m_mainFrame->m_pCurrentSceneObject);
        if( pDatasetInfo != NULL)
		{
			DatasetHelper* pDatasetHelper = m_mainFrame->m_pDatasetHelper;
			if(pDatasetHelper)
			{
				FibersGroup* pFibersGroup;
				pDatasetHelper->getFibersGroupDataset(pFibersGroup);
				if(pFibersGroup)
				{
					pFibersGroup->OnToggleIntensityBtn();
				}
			}
		}
	}
	m_mainFrame->refreshAllGLWidgets();
}

void PropertiesWindow::OnToggleOpacityBtn( wxEvent& WXUNUSED(event) )
{
	if (m_mainFrame->m_pCurrentSceneObject != NULL && m_mainFrame->m_currentListItem != -1)
    {
		DatasetInfo* pDatasetInfo = ((DatasetInfo*)m_mainFrame->m_pCurrentSceneObject);
        if( pDatasetInfo != NULL)
		{
			DatasetHelper* pDatasetHelper = m_mainFrame->m_pDatasetHelper;
			if(pDatasetHelper)
			{
				FibersGroup* pFibersGroup;
				pDatasetHelper->getFibersGroupDataset(pFibersGroup);
				if(pFibersGroup)
				{
					pFibersGroup->OnToggleOpacityBtn();
				}
			}
		}
	}
	m_mainFrame->refreshAllGLWidgets();
}

void PropertiesWindow::OnToggleMinMaxLengthBtn( wxEvent& WXUNUSED(event) )
{
	if (m_mainFrame->m_pCurrentSceneObject != NULL && m_mainFrame->m_currentListItem != -1)
    {
		DatasetInfo* pDatasetInfo = ((DatasetInfo*)m_mainFrame->m_pCurrentSceneObject);
        if( pDatasetInfo != NULL)
		{
			DatasetHelper* pDatasetHelper = m_mainFrame->m_pDatasetHelper;
			if(pDatasetHelper)
			{
				FibersGroup* pFibersGroup;
				pDatasetHelper->getFibersGroupDataset(pFibersGroup);
				if(pFibersGroup)
				{
					pFibersGroup->OnToggleMinMaxLengthBtn();
				}
			}
		}
	}
	m_mainFrame->refreshAllGLWidgets();
}

void PropertiesWindow::OnToggleSubsamplingBtn( wxEvent& WXUNUSED(event) )
{
	if (m_mainFrame->m_pCurrentSceneObject != NULL && m_mainFrame->m_currentListItem != -1)
    {
		DatasetInfo* pDatasetInfo = ((DatasetInfo*)m_mainFrame->m_pCurrentSceneObject);
        if( pDatasetInfo != NULL)
		{
			DatasetHelper* pDatasetHelper = m_mainFrame->m_pDatasetHelper;
			if(pDatasetHelper)
			{
				FibersGroup* pFibersGroup;
				pDatasetHelper->getFibersGroupDataset(pFibersGroup);
				if(pFibersGroup)
				{
					pFibersGroup->OnToggleSubsamplingBtn();
				}
			}
		}
	}
	m_mainFrame->refreshAllGLWidgets();
}

void PropertiesWindow::OnToggleCrossingFibersBtn( wxEvent& WXUNUSED(event) )
{
	if (m_mainFrame->m_pCurrentSceneObject != NULL && m_mainFrame->m_currentListItem != -1)
    {
		DatasetInfo* pDatasetInfo = ((DatasetInfo*)m_mainFrame->m_pCurrentSceneObject);
        if( pDatasetInfo != NULL)
		{
			DatasetHelper* pDatasetHelper = m_mainFrame->m_pDatasetHelper;
			if(pDatasetHelper)
			{
				FibersGroup* pFibersGroup;
				pDatasetHelper->getFibersGroupDataset(pFibersGroup);
				if(pFibersGroup)
				{
					pFibersGroup->OnToggleCrossingFibersBtn();
				}
			}
		}
	}
}

void PropertiesWindow::OnToggleColorModeBtn( wxEvent& WXUNUSED(event) )
{
	if (m_mainFrame->m_pCurrentSceneObject != NULL && m_mainFrame->m_currentListItem != -1)
    {
		DatasetInfo* pDatasetInfo = ((DatasetInfo*)m_mainFrame->m_pCurrentSceneObject);
        if( pDatasetInfo != NULL)
		{
			DatasetHelper* pDatasetHelper = m_mainFrame->m_pDatasetHelper;
			if(pDatasetHelper)
			{
				FibersGroup* pFibersGroup;
				pDatasetHelper->getFibersGroupDataset(pFibersGroup);
				if(pFibersGroup)
				{
					pFibersGroup->OnToggleColorModeBtn();
				}
			}
		}
	}
	m_mainFrame->refreshAllGLWidgets();
}

void PropertiesWindow::OnToggleLocalColoringBtn( wxEvent& WXUNUSED(event) )
{
	if (m_mainFrame->m_pCurrentSceneObject != NULL && m_mainFrame->m_currentListItem != -1)
    {
		DatasetInfo* pDatasetInfo = ((DatasetInfo*)m_mainFrame->m_pCurrentSceneObject);
        if( pDatasetInfo != NULL)
		{
			DatasetHelper* pDatasetHelper = m_mainFrame->m_pDatasetHelper;
			if(pDatasetHelper)
			{
				FibersGroup* pFibersGroup;
				pDatasetHelper->getFibersGroupDataset(pFibersGroup);
				if(pFibersGroup)
				{
					pFibersGroup->OnToggleLocalColoring();
				}
			}
		}
	}
	m_mainFrame->refreshAllGLWidgets();
}

void PropertiesWindow::OnToggleNormalColoringBtn( wxEvent& WXUNUSED(event) )
{
	if (m_mainFrame->m_pCurrentSceneObject != NULL && m_mainFrame->m_currentListItem != -1)
    {
		DatasetInfo* pDatasetInfo = ((DatasetInfo*)m_mainFrame->m_pCurrentSceneObject);
        if( pDatasetInfo != NULL)
		{
			DatasetHelper* pDatasetHelper = m_mainFrame->m_pDatasetHelper;
			if(pDatasetHelper)
			{
				FibersGroup* pFibersGroup;
				pDatasetHelper->getFibersGroupDataset(pFibersGroup);
				if(pFibersGroup)
				{
					pFibersGroup->OnToggleNormalColoring();
				}
			}
		}
	}
	m_mainFrame->refreshAllGLWidgets();
}

void PropertiesWindow::OnClickApplyBtn( wxEvent& WXUNUSED(event) )
{
	if (m_mainFrame->m_pCurrentSceneObject != NULL && m_mainFrame->m_currentListItem != -1)
    {
		DatasetInfo* pDatasetInfo = ((DatasetInfo*)m_mainFrame->m_pCurrentSceneObject);
        if( pDatasetInfo != NULL)
		{
			DatasetHelper* pDatasetHelper = m_mainFrame->m_pDatasetHelper;
			if(pDatasetHelper)
			{
				FibersGroup* pFibersGroup;
				pDatasetHelper->getFibersGroupDataset(pFibersGroup);
				if(pFibersGroup)
				{
					pFibersGroup->OnClickApplyBtn();
				}
			}
		}
	}
	m_mainFrame->refreshAllGLWidgets();
}

void PropertiesWindow::OnClickCancelBtn( wxEvent& WXUNUSED(event) )
{
	if (m_mainFrame->m_pCurrentSceneObject != NULL && m_mainFrame->m_currentListItem != -1)
    {
		DatasetInfo* pDatasetInfo = ((DatasetInfo*)m_mainFrame->m_pCurrentSceneObject);
        if( pDatasetInfo != NULL)
		{
			DatasetHelper* pDatasetHelper = m_mainFrame->m_pDatasetHelper;
			if(pDatasetHelper)
			{
				FibersGroup* pFibersGroup;
				pDatasetHelper->getFibersGroupDataset(pFibersGroup);
				if(pFibersGroup)
				{
					pFibersGroup->OnClickCancelBtn();
				}
			}
		}
	}
	m_mainFrame->refreshAllGLWidgets();
}

void PropertiesWindow::OnClickGenerateFiberVolumeBtn( wxEvent& WXUNUSED(event) )
{
	if (m_mainFrame->m_pCurrentSceneObject != NULL && m_mainFrame->m_currentListItem != -1)
    {
		DatasetInfo* pDatasetInfo = ((DatasetInfo*)m_mainFrame->m_pCurrentSceneObject);
        if( pDatasetInfo != NULL)
		{
			DatasetHelper* pDatasetHelper = m_mainFrame->m_pDatasetHelper;
			if(pDatasetHelper)
			{
				FibersGroup* pFibersGroup;
				pDatasetHelper->getFibersGroupDataset(pFibersGroup);
				if(pFibersGroup)
				{
					pFibersGroup->OnClickGenerateFiberVolumeBtn();
				}
			}
		}
	}
	m_mainFrame->refreshAllGLWidgets();
}

void PropertiesWindow::OnToggleShowFS( wxEvent& WXUNUSED(event) )
{
    Logger::getInstance()->print( _T( "Event triggered - PropertiesWindow::OnToggleShowFS" ), LOGLEVEL_DEBUG );

    if (m_mainFrame->m_pCurrentSceneObject != NULL && m_mainFrame->m_currentListItem != -1)
    {
        if( ! ((DatasetInfo*)m_mainFrame->m_pCurrentSceneObject)->toggleShowFS())
            m_mainFrame->m_pListCtrl->SetItem( m_mainFrame->m_currentListItem, 1, ((DatasetInfo*)m_mainFrame->m_pCurrentSceneObject)->getName().BeforeFirst( '.' ) + wxT( "*" ) );
        else
            m_mainFrame->m_pListCtrl->SetItem( m_mainFrame->m_currentListItem, 1, ((DatasetInfo*)m_mainFrame->m_pCurrentSceneObject)->getName().BeforeFirst( '.' ) );
        m_mainFrame->refreshAllGLWidgets();
    }
}

void PropertiesWindow::OnListItemShow( wxCommandEvent&  WXUNUSED(event) )
{
    Logger::getInstance()->print( _T( "Event triggered - PropertiesWindow::OnListItemShow" ), LOGLEVEL_DEBUG );

    if( m_mainFrame->m_pCurrentSceneObject == NULL && m_mainFrame->m_currentListItem != -1)
        return;

    DatasetInfo* l_info = (DatasetInfo*)m_mainFrame->m_pListCtrl->GetItemData( m_mainFrame->m_currentListItem );
    if( l_info->toggleShow() )
        m_mainFrame->m_pListCtrl->SetItem( m_mainFrame->m_currentListItem, 0, wxT( "" ), 0 );
    else
        m_mainFrame->m_pListCtrl->SetItem( m_mainFrame->m_currentListItem, 0, wxT( "" ), 1 );
	
	if(l_info->getType() == FIBERSGROUP)
	{
		FibersGroup* pFibersGroup = (FibersGroup*)l_info;
		if(pFibersGroup != NULL)
		{
			pFibersGroup->OnToggleVisibleBtn();
		}
	}

    m_mainFrame->refreshAllGLWidgets();
}


void PropertiesWindow::OnSliderIntensityThresholdMoved( wxCommandEvent& WXUNUSED(event) )
{
    if( m_mainFrame->m_pCurrentSceneObject != NULL && m_mainFrame->m_currentListItem != -1 )
    {
        DatasetInfo* l_current = (DatasetInfo*)m_mainFrame->m_pCurrentSceneObject;
        float l_threshold = (float)l_current->m_psliderThresholdIntensity->GetValue() / 100.0f;

        if( l_current->getUseTex() )
            m_mainFrame->m_pListCtrl->SetItem( m_mainFrame->m_currentListItem, 2, wxString::Format( wxT( "%.2f" ),   l_threshold * l_current->getOldMax() ) );
        else
            m_mainFrame->m_pListCtrl->SetItem( m_mainFrame->m_currentListItem, 2, wxString::Format( wxT( "(%.2f)" ), l_threshold * l_current->getOldMax() ) );

        l_current->setThreshold( l_threshold );
        if( l_current->getType() == SURFACE )
        {
            Surface* s = (Surface*)l_current;
            s->movePoints();
        }
        if( l_current->getType() == ISO_SURFACE && ! l_current->m_psliderThresholdIntensity->leftDown() )
        {
            CIsoSurface* s = (CIsoSurface*)l_current;
            s->GenerateWithThreshold();
        }
        if( l_current->getType() < RGB )
        {
            Anatomy* a = (Anatomy*)l_current;
            if( a->m_pRoi )
                a->m_pRoi->setThreshold( l_threshold );
        }
        // This slider will set the Brightness level. Currently only the glyphs uses this value.
        l_current->setBrightness( 1.0f - l_threshold );
        m_mainFrame->refreshAllGLWidgets();
    }
}

void PropertiesWindow::OnSliderOpacityThresholdMoved( wxCommandEvent& WXUNUSED(event) )
{
    if( m_mainFrame->m_pCurrentSceneObject != NULL && m_mainFrame->m_currentListItem != -1 )
    {
        DatasetInfo* l_current = (DatasetInfo*)m_mainFrame->m_pCurrentSceneObject;
        l_current->setAlpha( (float)l_current->m_psliderOpacity->GetValue() / 100.0f);
        m_mainFrame->refreshAllGLWidgets();
    }
}


void PropertiesWindow::OnEqualizeDataset( wxEvent& WXUNUSED(event) )
{
    Logger::getInstance()->print( _T( "Event triggered - PropertiesWindow::OnEqualizeDataset" ), LOGLEVEL_DEBUG );

    if( m_mainFrame->m_pCurrentSceneObject != NULL && m_mainFrame->m_currentListItem != -1 )
	{
        if( ((DatasetInfo*)m_mainFrame->m_pCurrentSceneObject)->getType() < MESH )
        {
            ((Anatomy*)m_mainFrame->m_pCurrentSceneObject)->toggleEqualization();
        }
	}
}

void PropertiesWindow::OnEqualizationSliderChange( wxCommandEvent& WXUNUSED(event) )
{
    Logger::getInstance()->print( _T( "Event triggered - PropertiesWindow::OnEqualizationSliderChange" ), LOGLEVEL_DEBUG );

    if( m_mainFrame->m_pCurrentSceneObject != NULL && m_mainFrame->m_currentListItem != -1 )
    {
        if( ((DatasetInfo*)m_mainFrame->m_pCurrentSceneObject)->getType() < MESH )
        {
            ((Anatomy*)m_mainFrame->m_pCurrentSceneObject)->equalizationSliderChange();
        }
    }
}

void PropertiesWindow::OnRename( wxCommandEvent& WXUNUSED(event) )
{
    if( m_mainFrame->m_pCurrentSceneObject != NULL && m_mainFrame->m_currentListItem != -1 )
    {
        wxTextEntryDialog dialog( this, _T( "Please enter a new name" ) );
        DatasetInfo* pInfo = (DatasetInfo*)m_mainFrame->m_pCurrentSceneObject;

        dialog.SetValue( pInfo->getName().BeforeFirst( '.' ) );

        wxString ext = pInfo->getName().AfterFirst( '.' );

        if( ( dialog.ShowModal() == wxID_OK ) && ( dialog.GetValue() != _T( "" ) ) )
        {
            pInfo->setName( dialog.GetValue() + wxT( "." ) + ext );

            //Change the name on the widget in the GUI
            long item = m_mainFrame->m_currentListItem;
            m_mainFrame->m_pListCtrl->SetItem(item, 1, pInfo->getName().BeforeFirst( '.' ) );

            DatasetInfo* info = ( (DatasetInfo*)m_mainFrame->m_pListCtrl->GetItemData( m_mainFrame->m_currentListItem ) );
            info->m_ptxtName->Clear();
            *info->m_ptxtName << pInfo->getName();
        }
    }

    m_mainFrame->refreshAllGLWidgets();
}

void PropertiesWindow::OnFlipX( wxCommandEvent& WXUNUSED(event) )
{
    DatasetInfo* pInfo = (DatasetInfo*)m_mainFrame->m_pListCtrl->GetItemData( m_mainFrame->m_currentListItem );
    pInfo->flipAxis(X_AXIS);
}

void PropertiesWindow::OnFlipY( wxCommandEvent& WXUNUSED(event) )
{
    DatasetInfo* pInfo = (DatasetInfo*)m_mainFrame->m_pListCtrl->GetItemData( m_mainFrame->m_currentListItem );
    pInfo->flipAxis(Y_AXIS);
}

void PropertiesWindow::OnFlipZ( wxCommandEvent& WXUNUSED(event) )
{
    DatasetInfo* pInfo = (DatasetInfo*)m_mainFrame->m_pListCtrl->GetItemData( m_mainFrame->m_currentListItem );
    pInfo->flipAxis(Z_AXIS);
}

void PropertiesWindow::OnDilateDataset( wxCommandEvent& WXUNUSED(event) )
{
    if( m_mainFrame->m_pCurrentSceneObject != NULL && m_mainFrame->m_currentListItem != -1 )
    {
        if( ((DatasetInfo*)m_mainFrame->m_pCurrentSceneObject)->getType() < MESH )
        {
            ((Anatomy*)m_mainFrame->m_pCurrentSceneObject)->dilate();
        }
    }
}

void PropertiesWindow::OnErodeDataset(wxCommandEvent& WXUNUSED(event))
{
    if( m_mainFrame->m_pCurrentSceneObject != NULL && m_mainFrame->m_currentListItem != -1 )
    {
        if( ((DatasetInfo*)m_mainFrame->m_pCurrentSceneObject)->getType() < MESH )
        {
            ((Anatomy*)m_mainFrame->m_pCurrentSceneObject)->erode();
        }
    }
}

void PropertiesWindow::OnMinimizeDataset( wxCommandEvent& WXUNUSED(event) )
{
    if( m_mainFrame->m_pCurrentSceneObject != NULL && m_mainFrame->m_currentListItem != -1 )
    {
        if( ((DatasetInfo*)m_mainFrame->m_pCurrentSceneObject)->getType() < MESH )
        {
            ((Anatomy*)m_mainFrame->m_pCurrentSceneObject)->minimize();
        }
    }
}

void PropertiesWindow::OnListItemCutOut( wxCommandEvent&  WXUNUSED(event) )
{
    m_mainFrame->m_pDatasetHelper->createCutDataset();
    m_mainFrame->refreshAllGLWidgets();
}

void PropertiesWindow::OnNewIsoSurface( wxCommandEvent& WXUNUSED(event) )
{
    m_mainFrame->m_pDatasetHelper->createIsoSurface();
}

void PropertiesWindow::OnNewOffsetSurface( wxCommandEvent& WXUNUSED(event ))
{
    m_mainFrame->m_pDatasetHelper->createDistanceMapAndIso();
}

void PropertiesWindow::OnNewDistanceMap (wxCommandEvent& WXUNUSED(event))
{
    m_mainFrame->m_pDatasetHelper->createDistanceMap();
}

void PropertiesWindow::OnNewVoiFromOverlay( wxCommandEvent& WXUNUSED(event) )
{
    SelectionObject *pSelectionObject  = NULL;
    Anatomy         *pAnatomy          = NULL;

    if(m_mainFrame->m_pCurrentSceneObject != NULL && m_mainFrame->m_currentListItem != -1)
    {
        if ( ((DatasetInfo*)m_mainFrame->m_pCurrentSceneObject)->getType() < RGB)
        {
            pAnatomy = (Anatomy*)m_mainFrame->m_pCurrentSceneObject;
            float trs = pAnatomy->getThreshold();
            if( trs == 0.0 )
                trs = 0.01f;

            pSelectionObject = new SelectionVOI( m_mainFrame->m_pDatasetHelper, pAnatomy, 
                                                  trs, THRESHOLD_GREATER_EQUAL );
        }
        else
        {
            return;
        }
    }
    else
    {
        return;
    }

    wxTreeItemId parentSelectionId = m_mainFrame->m_pTreeWidget->GetSelection();
    
    AddSelectionObjectToSelectionTree( pSelectionObject, parentSelectionId );

    m_mainFrame->refreshAllGLWidgets();
}

void PropertiesWindow::OnNewVoiFromClusters( wxCommandEvent& event )
{
    // First check is for the Shift key.
    bool shiftKeyDown = wxGetKeyState( WXK_SHIFT );
    
    // Get the anatomy
    Anatomy *pAnat( NULL );
    
    if( m_mainFrame->m_currentListItem != -1 )
    {
        DatasetInfo* pInfo = (DatasetInfo*)m_mainFrame->m_pListCtrl->GetItemData( m_mainFrame->m_currentListItem );
        if ( pInfo->getType() < RGB )
        {
            pAnat = (Anatomy*)pInfo;
            if( pAnat->getAnatType() != TYPE_FMRI_CLUSTERS_MAP )
            {
                return;
            }
        }
        else
        {
            return;
        }
    }
    else
    {
        return;
    }
            
    wxArrayString choices;
    set< float > clustersValues( pAnat->getValuesSet() );
    vector< float > clustersValuesArray( clustersValues.begin(), clustersValues.end() );
    
    for( set< float >::const_iterator valIt( clustersValues.begin() ); valIt != clustersValues.end(); ++valIt )
    {
        // Do not display the 0.0 value, since it does not represent a cluster.
        if( *valIt != 0.0f )
        {
            choices.Add( wxString::Format( wxT("%4f"), *valIt ) );
        }
    }

    wxMultiChoiceDialog clusterSelectionDiag( 0, wxT( "Please select the clusters for which you want to create the VOIs." ), wxT( "Clusters selection" ), choices );

    int returnValue = clusterSelectionDiag.ShowModal();
    
    if( returnValue == wxID_OK )
    {
        // Create the VOIs for each choice.
        wxArrayInt selectedIdx = clusterSelectionDiag.GetSelections();
        
        vector< SelectionObject * > selectionObjects;
        
        for( int idIdx( selectedIdx.GetCount() - 1 ); idIdx >= 0 ; --idIdx )
        {
            // Get the value, and create the VOI.
            // Need the + 1 since the 0.0 value is not displayed but still in the vector.
            float clusterValue = clustersValuesArray.at( selectedIdx.Item( idIdx ) + 1 );
            SelectionVOI *pSelVOI = new SelectionVOI( m_mainFrame->m_pDatasetHelper, pAnat, 
                                                      clusterValue, THRESHOLD_EQUAL );
            
            wxString convertedClusterVal;
            convertedClusterVal << clusterValue;
            wxString voiName( wxT( "[VOI] - [fMRI] ") + convertedClusterVal );
            
            pSelVOI->setName( voiName );
            
            selectionObjects.push_back( pSelVOI );
        }
        
        // If the Shift key was pressed when the button was clicked, the user wants to
        // place all selection VOIs in a hierarchical order. Else, in the default case,
        // we place the first one, then all the remaining ones are direct children of the
        // first one.
        AddSelectionObjectsToSelectionTree( selectionObjects, !shiftKeyDown );
    }
    
    m_mainFrame->refreshAllGLWidgets();
}

void PropertiesWindow::OnSegment(wxCommandEvent& WXUNUSED(event))
{
    m_mainFrame->m_pDatasetHelper->m_isSegmentActive = !m_mainFrame->m_pDatasetHelper->m_isSegmentActive;

    if(!m_mainFrame->m_pMainGL->object.empty())
        m_mainFrame->m_pMainGL->object.clear();
    if(!m_mainFrame->m_pMainGL->background.empty())
        m_mainFrame->m_pMainGL->background.clear();

    m_mainFrame->m_pDatasetHelper->m_isObjfilled = false;
    m_mainFrame->m_pDatasetHelper->m_isBckfilled = false;
    m_mainFrame->m_pDatasetHelper->m_isObjCreated = false;
    m_mainFrame->m_pDatasetHelper->m_isBckCreated = false;
}

void PropertiesWindow::OnFloodFill(wxCommandEvent& WXUNUSED(event))
{
    m_mainFrame->m_pDatasetHelper->m_SegmentMethod = 0;
    m_mainFrame->m_pDatasetHelper->m_isFloodfillActive = true;
    m_mainFrame->m_pDatasetHelper->m_isSelectBckActive = false;
    m_mainFrame->m_pDatasetHelper->m_isSelectObjActive = false;
    ((Anatomy*)m_mainFrame->m_pCurrentSceneObject)->toggleSegment();
}

void PropertiesWindow::OnSliderFloodMoved( wxCommandEvent& WXUNUSED(event) )
{
    float l_sliderValue = ((Anatomy*)m_mainFrame->m_pCurrentSceneObject)->m_pSliderFlood->GetValue() / 200.0f;
    ((Anatomy*)m_mainFrame->m_pCurrentSceneObject)->setFloodThreshold(l_sliderValue);
    ((Anatomy*)m_mainFrame->m_pCurrentSceneObject)->m_pTxtThresBox->SetValue(wxString::Format( wxT( "%.2f"), l_sliderValue));
}

void PropertiesWindow::OnSliderGraphSigmaMoved( wxCommandEvent& WXUNUSED(event) )
{
    ((Anatomy*)m_mainFrame->m_pCurrentSceneObject)->setGraphSigma(((Anatomy*)m_mainFrame->m_pCurrentSceneObject)->m_pSliderGraphSigma->GetValue());
    std::cout << (((Anatomy*)m_mainFrame->m_pCurrentSceneObject)->m_pSliderGraphSigma->GetValue()) << endl;
}

void PropertiesWindow::OnKmeans( wxCommandEvent& WXUNUSED(event) )
{
    m_mainFrame->m_pDatasetHelper->m_SegmentMethod = 2;
    m_mainFrame->m_pMainGL->segment();
}

void PropertiesWindow::OnSelectObj(wxCommandEvent& WXUNUSED(event))
{
    Logger::getInstance()->print( _T( "Event triggered - PropertiesWindow::OnSelectObj" ), LOGLEVEL_DEBUG );

    m_mainFrame->m_pDatasetHelper->m_SegmentMethod = 1;
    m_mainFrame->m_pDatasetHelper->m_isSelectBckActive = false;
    m_mainFrame->m_pDatasetHelper->m_isFloodfillActive = false;
    m_mainFrame->m_pDatasetHelper->m_isSelectObjActive = true;
}

void PropertiesWindow::OnSelectBck(wxCommandEvent& WXUNUSED(event))
{
    Logger::getInstance()->print( _T( "Event triggered - PropertiesWindow::OnSelectBck" ), LOGLEVEL_DEBUG );

    m_mainFrame->m_pDatasetHelper->m_SegmentMethod = 1;
    m_mainFrame->m_pDatasetHelper->m_isFloodfillActive = false;
    m_mainFrame->m_pDatasetHelper->m_isSelectBckActive = true;
    m_mainFrame->m_pDatasetHelper->m_isSelectObjActive = false;    
}

void PropertiesWindow::OnbtnGraphCut(wxCommandEvent& WXUNUSED(event))
{
    m_mainFrame->m_pDatasetHelper->m_SegmentMethod = 1;
    m_mainFrame->m_pDatasetHelper->m_isFloodfillActive = false;
    m_mainFrame->m_pMainGL->segment();    
}


void PropertiesWindow::OnClean( wxCommandEvent& WXUNUSED(event) )
{
    if(m_mainFrame->m_pCurrentSceneObject != NULL && m_mainFrame->m_currentListItem != -1 )
    {
        DatasetInfo* l_info = (DatasetInfo*)m_mainFrame->m_pCurrentSceneObject;
        if( l_info->getType() == MESH || l_info->getType() == ISO_SURFACE)
            l_info->clean();
    }
    m_mainFrame->refreshAllGLWidgets();
}

void PropertiesWindow::OnLoop( wxCommandEvent& WXUNUSED(event) )
{
    if( m_mainFrame->m_pCurrentSceneObject != NULL && m_mainFrame->m_currentListItem != -1 )
    {
        ((DatasetInfo*)m_mainFrame->m_pCurrentSceneObject)->smooth();
    }
    m_mainFrame->refreshAllGLWidgets();
}


void PropertiesWindow::OnToggleLIC( wxCommandEvent& WXUNUSED(event) )
{
    Logger::getInstance()->print( _T( "Event triggered - PropertiesWindow::OnToggleLIC" ), LOGLEVEL_DEBUG );

    if( m_mainFrame->m_pCurrentSceneObject != NULL && m_mainFrame->m_currentListItem != -1 && m_mainFrame->m_pDatasetHelper->m_vectorsLoaded )
    {
        ((DatasetInfo*) m_mainFrame->m_pCurrentSceneObject)->activateLIC();
    }
    m_mainFrame->refreshAllGLWidgets();
}

void PropertiesWindow::OnToggleDrawPointsMode( wxCommandEvent& event )
{
    Logger::getInstance()->print( _T( "Event triggered - PropertiesWindow::OnToggleDrawPointsMode" ), LOGLEVEL_DEBUG );

    m_mainFrame->onToggleDrawPointsMode(event);
}

void PropertiesWindow::OnMoveBoundaryPointsLeft( wxCommandEvent& event )
{
    m_mainFrame->onMoveBoundaryPointsLeft(event);
}

void PropertiesWindow::OnMoveBoundaryPointsRight(wxCommandEvent& event)
{
    m_mainFrame->onMoveBoundaryPointsRight(event);
}

void PropertiesWindow::OnFibersFilter( wxCommandEvent& event)
{
    Fibers* pTmpFib = NULL;
    m_mainFrame->m_pDatasetHelper->getSelectedFiberDataset(pTmpFib);    
    if(pTmpFib != NULL)
    {
        pTmpFib->updateFibersFilters();
    }
}

void PropertiesWindow::OnGenerateFiberVolume( wxCommandEvent& WXUNUSED(event) )
{
    Fibers* pTmpFib = NULL;
    m_mainFrame->m_pDatasetHelper->getSelectedFiberDataset(pTmpFib);
    if(pTmpFib != NULL)
    {
        pTmpFib->generateFiberVolume();
    }
}

void PropertiesWindow::OnListMenuThreshold( wxCommandEvent&  WXUNUSED(event) )
{
    Logger::getInstance()->print( _T( "Event triggered - PropertiesWindow::OnListMenuThreshold" ), LOGLEVEL_DEBUG );

    if( m_mainFrame->m_pCurrentSceneObject == NULL && m_mainFrame->m_currentListItem != -1)
        return;
    DatasetInfo* l_info = (DatasetInfo*)m_mainFrame->m_pCurrentSceneObject;
    if( l_info->getType() >= MESH )
    {
        if( ! l_info->toggleUseTex() )
            m_mainFrame->m_pListCtrl->SetItem( m_mainFrame->m_currentListItem, 2, wxT( "(" ) + wxString::Format( wxT( "%.2f" ), l_info->getThreshold() * l_info->getOldMax()) + wxT( ")" ) );
        else
            m_mainFrame->m_pListCtrl->SetItem( m_mainFrame->m_currentListItem, 2, wxString::Format( wxT( "%.2f" ), l_info->getThreshold() * l_info->getOldMax() ) );
    }
}

//////////////////////////////////////////////////////////////////////////
// This function will be called when the distance coloring option is
// selected when right-clicking on a fiber.
//////////////////////////////////////////////////////////////////////////
void PropertiesWindow::OnListMenuDistance( wxCommandEvent& WXUNUSED(event))
{
	Fibers* pFibers = NULL;
    m_mainFrame->m_pDatasetHelper->getSelectedFiberDataset(pFibers);
    if(pFibers != NULL)
	{
		Logger::getInstance()->print( _T( "Event triggered - PropertiesWindow::OnListMenuDistance" ), LOGLEVEL_DEBUG );

		if(pFibers->getColorationMode() != DISTANCE_COLOR)
		{
			pFibers->setColorationMode( DISTANCE_COLOR );
			pFibers->updateFibersColors();  
			pFibers->updateColorationMode();
		}
	}
}

//////////////////////////////////////////////////////////////////////////
// This function will be called when the minimum distance coloring option is
// selected when right-clicking on a fiber.
//////////////////////////////////////////////////////////////////////////
void PropertiesWindow::OnListMenuMinDistance( wxCommandEvent& WXUNUSED(event))
{
    Logger::getInstance()->print( _T( "Event triggered - PropertiesWindow::OnListMenuMinDistance" ), LOGLEVEL_DEBUG );

	Fibers* pFibers = NULL;
    m_mainFrame->m_pDatasetHelper->getSelectedFiberDataset(pFibers);
    if(pFibers != NULL)
	{

		if(pFibers->getColorationMode() != MINDISTANCE_COLOR)
		{
			pFibers->setColorationMode( MINDISTANCE_COLOR );
			pFibers->updateFibersColors();
			pFibers->updateColorationMode();
		}
    }
}


///////////////////////////////////////////////////////////////////////////
// This function will be triggered when the user click on the color with curvature button
// button that is located in the m_fibersInfoSizer.
///////////////////////////////////////////////////////////////////////////
void PropertiesWindow::OnColorWithCurvature( wxCommandEvent& WXUNUSED(event) )
{
	Fibers* pFibers = NULL;
    m_mainFrame->m_pDatasetHelper->getSelectedFiberDataset(pFibers);
    if(pFibers != NULL)
    {
		if(pFibers->getColorationMode() != CURVATURE_COLOR)
		{
			pFibers->setColorationMode( CURVATURE_COLOR );
			pFibers->updateFibersColors();
			pFibers->updateColorationMode();
		}
    }
}

///////////////////////////////////////////////////////////////////////////
// This function will be triggered when the user click on the display min/max cross section
// button that is located in the m_fibersInfoSizer.
///////////////////////////////////////////////////////////////////////////
void PropertiesWindow::OnColorWithTorsion( wxCommandEvent& WXUNUSED(event) )
{
	Fibers* pFibers = NULL;
    m_mainFrame->m_pDatasetHelper->getSelectedFiberDataset(pFibers);
    if(pFibers != NULL)
    {
		if(pFibers->getColorationMode() != TORSION_COLOR)
		{
			pFibers->setColorationMode( TORSION_COLOR );
			pFibers->updateFibersColors();
			pFibers->updateColorationMode();
		}
    }
}

void PropertiesWindow::OnNormalColoring( wxCommandEvent& WXUNUSED(event) )
{
	Fibers* pFibers = NULL;
    m_mainFrame->m_pDatasetHelper->getSelectedFiberDataset(pFibers);
    if(pFibers != NULL)
    {
		if(pFibers->getColorationMode() != NORMAL_COLOR)
		{
			pFibers->setColorationMode( NORMAL_COLOR );
			pFibers->updateFibersColors();
			pFibers->updateColorationMode();
		}
    }
}

///////////////////////////////////////////////////////////////////////////
// This function will be triggered when the user click on the normal coloring radio
// button located in the mean fiber coloring option
///////////////////////////////////////////////////////////////////////////
void PropertiesWindow::OnNormalMeanFiberColoring( wxCommandEvent& event )
{
   ( (SelectionObject*) m_mainFrame->m_pCurrentSceneObject )->setMeanFiberColorMode(NORMAL_COLOR); 
}

///////////////////////////////////////////////////////////////////////////
// This function will be triggered when the user click on the custom coloring radio
// button located in the mean fiber coloring option
///////////////////////////////////////////////////////////////////////////
void PropertiesWindow::OnCustomMeanFiberColoring( wxCommandEvent& event )
{
    ( (SelectionObject*) m_mainFrame->m_pCurrentSceneObject )->setMeanFiberColorMode(CUSTOM_COLOR);
}

///////////////////////////////////////////////////////////////////////////
// This function will be triggered when the user move the slider
// button located in the mean fiber coloring option
///////////////////////////////////////////////////////////////////////////
void PropertiesWindow::OnMeanFiberOpacityChange( wxCommandEvent& event )
{
    ( (SelectionObject*) m_mainFrame->m_pCurrentSceneObject )->updateMeanFiberOpacity();
}

//////////////////////////////////////////////////////////////////////////
// This function will be triggered when the user clicks on the "Set as distance
// anchor" option.
//////////////////////////////////////////////////////////////////////////
void PropertiesWindow::OnDistanceAnchorSet( wxCommandEvent& event )
{
    if (m_mainFrame->m_pDatasetHelper->m_lastSelectedObject!=NULL)
    {
        m_mainFrame->m_pDatasetHelper->m_lastSelectedObject->UseForDistanceColoring(!m_mainFrame->m_pDatasetHelper->m_lastSelectedObject->IsUsedForDistanceColoring());
        ColorFibers();
    }
}

///////////////////////////////////////////////////////////////////////////
// This function will call the updateFibersColors function on the currently loaded fiber set.
///////////////////////////////////////////////////////////////////////////
void PropertiesWindow::ColorFibers()
{   
    if (m_mainFrame->m_pCurrentSceneObject != NULL && m_mainFrame->m_currentListItem != -1)
    {
        ((Fibers*)m_mainFrame->m_pCurrentSceneObject)->updateFibersColors();  
    }  
}
void PropertiesWindow::OnOriginalShBasis( wxCommandEvent& WXUNUSED(event) )
{
    ODFs* l_dataset = new ODFs( m_mainFrame->m_pDatasetHelper );
    ((ODFs*)m_mainFrame->m_pCurrentSceneObject)->changeShBasis(l_dataset, m_mainFrame->m_pDatasetHelper, 0);

}

void PropertiesWindow::OnDescoteauxShBasis( wxCommandEvent& WXUNUSED(event) )
{
    ODFs* l_dataset = new ODFs( m_mainFrame->m_pDatasetHelper );
    ((ODFs*)m_mainFrame->m_pCurrentSceneObject)->changeShBasis(l_dataset, m_mainFrame->m_pDatasetHelper, 1);

}

void PropertiesWindow::OnTournierShBasis( wxCommandEvent& WXUNUSED(event) )
{
    ODFs* l_dataset = new ODFs( m_mainFrame->m_pDatasetHelper );
    ((ODFs*)m_mainFrame->m_pCurrentSceneObject)->changeShBasis(l_dataset, m_mainFrame->m_pDatasetHelper, 2);

}

void PropertiesWindow::OnPTKShBasis( wxCommandEvent& WXUNUSED(event) )
{
    ODFs* l_dataset = new ODFs( m_mainFrame->m_pDatasetHelper );
    ((ODFs*)m_mainFrame->m_pCurrentSceneObject)->changeShBasis(l_dataset, m_mainFrame->m_pDatasetHelper, 3);
}

///////////////////////////////////////////////////////////////////////////
// This function will be called when the slider for the min hue value in 
// the glyph options panel moved.
///////////////////////////////////////////////////////////////////////////
void PropertiesWindow::OnGlyphMinHueSliderMoved( wxCommandEvent& WXUNUSED(event) )
{
    updateGlyphColoration( MIN_HUE, ((Glyph*)m_mainFrame->m_pCurrentSceneObject)->m_psliderMinHueValue->GetValue() / 100.0f );
}

///////////////////////////////////////////////////////////////////////////
// This function will be called when the slider for the max hue value in 
// the glyph options panel moved.
///////////////////////////////////////////////////////////////////////////
void PropertiesWindow::OnGlyphMaxHueSliderMoved( wxCommandEvent& WXUNUSED(event) )
{
    updateGlyphColoration( MAX_HUE, ((Glyph*)m_mainFrame->m_pCurrentSceneObject)->m_psliderMaxHueValue->GetValue() / 100.0f );
}

///////////////////////////////////////////////////////////////////////////
// This function will be called when the slider for the saturation value in 
// the glyph options panel moved.
///////////////////////////////////////////////////////////////////////////
void PropertiesWindow::OnGlyphSaturationSliderMoved( wxCommandEvent& WXUNUSED(event) )
{
    updateGlyphColoration( SATURATION,((Glyph*)m_mainFrame->m_pCurrentSceneObject)->m_psliderSaturationValue->GetValue() / 100.0f );
}

///////////////////////////////////////////////////////////////////////////
// This function will be called when the slider for the luminance value in 
// the glyph options panel moved.
///////////////////////////////////////////////////////////////////////////
void PropertiesWindow::OnGlyphLuminanceSliderMoved( wxCommandEvent& WXUNUSED(event) )
{
    updateGlyphColoration( LUMINANCE, ((Glyph*)m_mainFrame->m_pCurrentSceneObject)->m_psliderLuminanceValue->GetValue() / 100.0f );
}

///////////////////////////////////////////////////////////////////////////
// This function will set the value of a glyph color modifier by the value
// on its corresponding slider.
//
// i_modifier       : The modifier indicating what GlyphColorModifier needs to be updated.
// i_value          : The value of the modifier to set.
///////////////////////////////////////////////////////////////////////////
void PropertiesWindow::updateGlyphColoration( GlyphColorModifier i_modifier, float i_value )
{
    if( m_mainFrame->m_pCurrentSceneObject != NULL && m_mainFrame->m_currentListItem != -1 )
    {            
        DatasetInfo* l_info = (DatasetInfo*)m_mainFrame->m_pCurrentSceneObject;    
        if( l_info->getType() == TENSORS || l_info->getType() == ODFS )
            ( (Glyph*)l_info )->setColor( i_modifier, i_value );
    }
}

///////////////////////////////////////////////////////////////////////////
// This function will set the LOD of a glyph by the value of its slider.
///////////////////////////////////////////////////////////////////////////
void PropertiesWindow::OnGlyphLODSliderMoved( wxCommandEvent& WXUNUSED(event) )
{
    if( m_mainFrame->m_pCurrentSceneObject != NULL && m_mainFrame->m_currentListItem != -1 )
    {            
        DatasetInfo* l_info = (DatasetInfo*)m_mainFrame->m_pCurrentSceneObject;
        if( l_info->getType() == TENSORS || l_info->getType() == ODFS )
        {
            ( (Glyph*)l_info )->setLOD( (LODChoices)((Glyph*)m_mainFrame->m_pCurrentSceneObject)->m_psliderLODValue->GetValue() );
        }
    }
}

///////////////////////////////////////////////////////////////////////////
// This function will set the light attenuation of a glyph by the value of its slider.
///////////////////////////////////////////////////////////////////////////
void PropertiesWindow::OnGlyphLightAttenuationSliderMoved( wxCommandEvent& WXUNUSED(event) )
{
    if( m_mainFrame->m_pCurrentSceneObject != NULL && m_mainFrame->m_currentListItem != -1 )
    {            
        DatasetInfo* l_info = (DatasetInfo*)m_mainFrame->m_pCurrentSceneObject;
        if( l_info->getType() == TENSORS || l_info->getType() == ODFS )
            ( (Glyph*)l_info )->setLighAttenuation( ((Glyph*)m_mainFrame->m_pCurrentSceneObject)->m_psliderLightAttenuation->GetValue() / 100.0f );
    }
}

///////////////////////////////////////////////////////////////////////////
// This function will be called when the light x position slider moved.
///////////////////////////////////////////////////////////////////////////
void PropertiesWindow::OnGlyphLightXDirectionSliderMoved( wxCommandEvent& WXUNUSED(event) )
{
    if( m_mainFrame->m_pCurrentSceneObject != NULL && m_mainFrame->m_currentListItem != -1 )
    {
        OnGlyphLightPositionChanged( X_AXIS, ((Glyph*)m_mainFrame->m_pCurrentSceneObject)->m_psliderLightXPosition->GetValue() / 100.0f  );
    }
}

///////////////////////////////////////////////////////////////////////////
// This function will be called when the light y position slider moved.
///////////////////////////////////////////////////////////////////////////
void PropertiesWindow::OnGlyphLightYDirectionSliderMoved( wxCommandEvent& WXUNUSED(event) )
{
    if( m_mainFrame->m_pCurrentSceneObject != NULL && m_mainFrame->m_currentListItem != -1 )
    {
        OnGlyphLightPositionChanged( Y_AXIS, ((Glyph*)m_mainFrame->m_pCurrentSceneObject)->m_psliderLightYPosition->GetValue() / 100.0f  );
    }
}

///////////////////////////////////////////////////////////////////////////
// This function will be called when the light z position slider moved.
///////////////////////////////////////////////////////////////////////////
void PropertiesWindow::OnGlyphLightZDirectionSliderMoved( wxCommandEvent& WXUNUSED(event) )
{
    if( m_mainFrame->m_pCurrentSceneObject != NULL && m_mainFrame->m_currentListItem != -1 )
    {
        OnGlyphLightPositionChanged( Z_AXIS, ((Glyph*)m_mainFrame->m_pCurrentSceneObject)->m_psliderLightZPosition->GetValue() / 100.0f  );
    }
}

///////////////////////////////////////////////////////////////////////////
// This function will set the light position for the proper axis.
//
// i_axisType       : The axis that we want to set the lght position for.
// i_position       : The value of the position.
///////////////////////////////////////////////////////////////////////////
void PropertiesWindow::OnGlyphLightPositionChanged( AxisType i_axisType, float i_position )
{
    if( m_mainFrame->m_pCurrentSceneObject != NULL && m_mainFrame->m_currentListItem != -1 )
    {
        ((Glyph*)m_mainFrame->m_pCurrentSceneObject)->setLightPosition( i_axisType, i_position);
    }
}

///////////////////////////////////////////////////////////////////////////
// This function will set the display value of a glyph by the value of its slider.
///////////////////////////////////////////////////////////////////////////
void PropertiesWindow::OnGlyphDisplaySliderMoved( wxCommandEvent& WXUNUSED(event) )
{
    if( m_mainFrame->m_pCurrentSceneObject != NULL && m_mainFrame->m_currentListItem != -1 )
    {
        ((Glyph*)m_mainFrame->m_pCurrentSceneObject)->setDisplayFactor( ((Glyph*)m_mainFrame->m_pCurrentSceneObject)->m_psliderDisplayValue->GetValue());
    }
}

///////////////////////////////////////////////////////////////////////////
// This function will set the scaling factor of a glyph by the value of its slider.
///////////////////////////////////////////////////////////////////////////
void PropertiesWindow::OnGlyphScalingFactorSliderMoved( wxCommandEvent& WXUNUSED(event) )
{
    if( m_mainFrame->m_pCurrentSceneObject != NULL && m_mainFrame->m_currentListItem != -1 )
    {
        ((Glyph*)m_mainFrame->m_pCurrentSceneObject)->setScalingFactor( ((Glyph*)m_mainFrame->m_pCurrentSceneObject)->m_psliderScalingFactor->GetValue()/10.0f);
    }
}

///////////////////////////////////////////////////////////////////////////
// This function will be called when the x flip check box in the 
// glyph options is checked/unchecked.
///////////////////////////////////////////////////////////////////////////
void PropertiesWindow::OnGlyphXAxisFlipChecked( wxCommandEvent& event )
{
    OnGlyphFlip( X_AXIS, event.IsChecked() );
}

///////////////////////////////////////////////////////////////////////////
// This function will be called when the y flip check box in the 
// glyph options is checked/unchecked.
///////////////////////////////////////////////////////////////////////////
void PropertiesWindow::OnGlyphYAxisFlipChecked( wxCommandEvent& event )
{
    OnGlyphFlip( Y_AXIS, event.IsChecked() );
}

///////////////////////////////////////////////////////////////////////////
// This function will be called when the z flip check box in the 
// glyph options is checked/unchecked.
///////////////////////////////////////////////////////////////////////////
void PropertiesWindow::OnGlyphZAxisFlipChecked( wxCommandEvent& event )
{
    OnGlyphFlip( Z_AXIS, event.IsChecked() );
}

///////////////////////////////////////////////////////////////////////////
// This function will simply find the currently displayed glyph and call 
// the flipAxis function with the proper parameter.
//
// i_axisType               : Determines on what axis we want to do the flip.
// i_isChecked              : Determines if the item is checked or not.
///////////////////////////////////////////////////////////////////////////
void PropertiesWindow::OnGlyphFlip( AxisType i_axisType, bool i_isChecked )
{
    if( m_mainFrame->m_pCurrentSceneObject != NULL && m_mainFrame->m_currentListItem != -1 )
    {
        ((Glyph*)m_mainFrame->m_pCurrentSceneObject)->flipAxis( i_axisType, i_isChecked );
    }
}

///////////////////////////////////////////////////////////////////////////
// This function will be called when the map on sphere radio button in the 
// glyph options is selected.
///////////////////////////////////////////////////////////////////////////
void PropertiesWindow::OnGlyphMapOnSphereSelected( wxCommandEvent& event )
{
    if( m_mainFrame->m_pCurrentSceneObject != NULL && m_mainFrame->m_currentListItem != -1 )
    {
        ((Glyph*)m_mainFrame->m_pCurrentSceneObject)->setDisplayShape( SPHERE );
    }
}

///////////////////////////////////////////////////////////////////////////
// This function will be called when the normal display radio button in the 
// glyph options is selected.
///////////////////////////////////////////////////////////////////////////
void PropertiesWindow::OnGlyphNormalSelected( wxCommandEvent& event )
{
    if( m_mainFrame->m_pCurrentSceneObject != NULL && m_mainFrame->m_currentListItem != -1 )
    {
        ((Glyph*)m_mainFrame->m_pCurrentSceneObject)->setDisplayShape( NORMAL );
    }
}

///////////////////////////////////////////////////////////////////////////
// This function will be called when the axes display radio button in the 
// glyph options is selected.
///////////////////////////////////////////////////////////////////////////
void PropertiesWindow::OnGlyphAxesSelected( wxCommandEvent& event )
{
    if( m_mainFrame->m_pCurrentSceneObject != NULL && m_mainFrame->m_currentListItem != -1 )
    {
        ((Glyph*)m_mainFrame->m_pCurrentSceneObject)->setDisplayShape( AXES );
    }
}

///////////////////////////////////////////////////////////////////////////
// This function will be called when the main axis display radio button in the 
// glyph options is selected.
///////////////////////////////////////////////////////////////////////////
void PropertiesWindow::OnGlyphMainAxisSelected( wxCommandEvent& event )
{
    if( m_mainFrame->m_pCurrentSceneObject != NULL && m_mainFrame->m_currentListItem != -1 )
    {
        if(((DatasetInfo*)m_mainFrame->m_pCurrentSceneObject)->getType() == ODFS && !((ODFs*)m_mainFrame->m_pCurrentSceneObject)->m_isMaximasSet)
        {
            ((ODFs*)m_mainFrame->m_pCurrentSceneObject)->extractMaximas();
        }
        ((Glyph*)m_mainFrame->m_pCurrentSceneObject)->setDisplayShape( AXIS );
        ((Glyph*)m_mainFrame->m_pCurrentSceneObject)->updatePropertiesSizer();
    }
}

///////////////////////////////////////////////////////////////////////////
// This function will be called when the color with position check box in the 
// glyph options is checked/unchecked.
///////////////////////////////////////////////////////////////////////////
void PropertiesWindow::OnGlyphColorWithPosition( wxCommandEvent& event )
{
    if( m_mainFrame->m_pCurrentSceneObject != NULL && m_mainFrame->m_currentListItem != -1 )
    {
        ((Glyph*)m_mainFrame->m_pCurrentSceneObject)->setColorWithPosition( event.IsChecked() );
    }
}

void PropertiesWindow::OnNormalizeTensors( wxCommandEvent& event )
{
    if( m_mainFrame->m_pCurrentSceneObject != NULL && m_mainFrame->m_currentListItem != -1 )
    {        
        ((Tensors*)m_mainFrame->m_pCurrentSceneObject)->normalize();        
    }
}

///////////////////////////////////////////////////////////////////////////
// This function will be triggered when the user click on the display fibers 
// info after a right click in the tree on a selectio object.
///////////////////////////////////////////////////////////////////////////
void PropertiesWindow::OnDisplayFibersInfo( wxCommandEvent& WXUNUSED(event) )
{
    Logger::getInstance()->print( _T( "Event triggered - PropertiesWindow::OnDisplayFibersInfo" ), LOGLEVEL_DEBUG );

// TODO remove when the bug with the wxChoice in Windows is fixed.
#ifndef __WXMSW__
    ((SelectionObject*)m_mainFrame->m_pCurrentSceneObject)->UpdateMeanValueTypeBox();
#endif
    ((SelectionObject*)m_mainFrame->m_pCurrentSceneObject)->SetFiberInfoGridValues();
    m_mainFrame->refreshAllGLWidgets();
}

///////////////////////////////////////////////////////////////////////////
// This function will be triggered when the user click on the display mean fiber 
// button that is located in the m_fibersInfoSizer.
///////////////////////////////////////////////////////////////////////////
void PropertiesWindow::OnDisplayMeanFiber( wxCommandEvent& WXUNUSED(event) )
{
    Logger::getInstance()->print( _T( "Event triggered - PropertiesWindow::OnDisplayMeanFiber" ), LOGLEVEL_DEBUG );

    ( (SelectionObject*)m_mainFrame->m_pCurrentSceneObject )->computeMeanFiber();
    m_mainFrame->refreshAllGLWidgets();
}

///////////////////////////////////////////////////////////////////////////
// This function will be triggered when the user click on the display convex hull
// button that is located in the m_fibersInfoSizer.
///////////////////////////////////////////////////////////////////////////
void PropertiesWindow::OnDisplayConvexHull( wxCommandEvent& WXUNUSED(event) )
{
    ( (SelectionObject*)m_mainFrame->m_pCurrentSceneObject )->computeConvexHull();
    m_mainFrame->refreshAllGLWidgets();
}

///////////////////////////////////////////////////////////////////////////
// This function will be triggered when the user click on the color button
// beside the display convex hull button that is located in the m_fibersInfoSizer.
///////////////////////////////////////////////////////////////////////////
void PropertiesWindow::OnConvexHullColorChange( wxCommandEvent& WXUNUSED(event) )
{
    if( ! m_mainFrame->m_pDatasetHelper->m_theScene )
        return;

    wxColourData l_colorData;

    for( int i = 0; i < 10; ++i )
    {
        wxColour l_color(i * 28, i * 28, i * 28);
        l_colorData.SetCustomColour(i, l_color);
    }

    int i = 10;
    wxColour l_color ( 255, 0, 0 );
    l_colorData.SetCustomColour( i++, l_color );
    wxColour l_color1( 0, 255, 0 );
    l_colorData.SetCustomColour( i++, l_color1 );
    wxColour l_color2( 0, 0, 255 );
    l_colorData.SetCustomColour( i++, l_color2 );
    wxColour l_color3( 255, 255, 0 );
    l_colorData.SetCustomColour( i++, l_color3 );
    wxColour l_color4( 255, 0, 255 );
    l_colorData.SetCustomColour( i++, l_color4 );
    wxColour l_color5( 0, 255, 255 );
    l_colorData.SetCustomColour( i++, l_color5 );

    wxColourDialog dialog( this, &l_colorData );
    wxColour l_col;
    if( dialog.ShowModal() == wxID_OK )
    {
        wxColourData l_retData = dialog.GetColourData();
        l_col = l_retData.GetColour();
    }
    else
    {
        return;
    }

    ((SelectionObject*)m_mainFrame->m_pCurrentSceneObject)->setConvexHullColor( l_col );
  
    m_mainFrame->refreshAllGLWidgets();
}

void PropertiesWindow::OnConvexHullOpacityChange( wxCommandEvent& WXUNUSED(event) )
{
    ( (SelectionObject*) m_mainFrame->m_pCurrentSceneObject )->updateConvexHullOpacity();
}

///////////////////////////////////////////////////////////////////////////
// This function will be triggered when the user click on the color palette
// button that is located aside of the Show mean fiber button
///////////////////////////////////////////////////////////////////////////
void PropertiesWindow::OnMeanFiberColorChange( wxCommandEvent& WXUNUSED(event) )
{
    if( ! m_mainFrame->m_pDatasetHelper->m_theScene )
        return;

    wxColourData l_colorData;

    for( int i = 0; i < 10; ++i )
    {
        wxColour l_color(i * 28, i * 28, i * 28);
        l_colorData.SetCustomColour(i, l_color);
    }

    int i = 10;
    wxColour l_color ( 255, 0, 0 );
    l_colorData.SetCustomColour( i++, l_color );
    wxColour l_color1( 0, 255, 0 );
    l_colorData.SetCustomColour( i++, l_color1 );
    wxColour l_color2( 0, 0, 255 );
    l_colorData.SetCustomColour( i++, l_color2 );
    wxColour l_color3( 255, 255, 0 );
    l_colorData.SetCustomColour( i++, l_color3 );
    wxColour l_color4( 255, 0, 255 );
    l_colorData.SetCustomColour( i++, l_color4 );
    wxColour l_color5( 0, 255, 255 );
    l_colorData.SetCustomColour( i++, l_color5 );

    wxColourDialog dialog( this, &l_colorData );
    wxColour l_col;
    if( dialog.ShowModal() == wxID_OK )
    {
        wxColourData l_retData = dialog.GetColourData();
        l_col = l_retData.GetColour();
    }
    else
    {
        return;
    }

    ((SelectionObject*)m_mainFrame->m_pCurrentSceneObject)->setMeanFiberColor( l_col);
    
    m_mainFrame->refreshAllGLWidgets();
}

///////////////////////////////////////////////////////////////////////////
// This function will be triggered when the user click on the display cross sections
// button that is located in the m_fibersInfoSizer.
///////////////////////////////////////////////////////////////////////////
void PropertiesWindow::OnDisplayCrossSections( wxCommandEvent& WXUNUSED(event) )
{
    Logger::getInstance()->print( _T( "Event triggered - PropertiesWindow::OnDisplayCrossSections" ), LOGLEVEL_DEBUG );

    ((SelectionObject*)m_mainFrame->m_pCurrentSceneObject)->m_displayCrossSections = (CrossSectionsDisplay)( ( (int)((SelectionObject*)m_mainFrame->m_pCurrentSceneObject)->m_displayCrossSections ) + 1 );
    if( ((SelectionObject*)m_mainFrame->m_pCurrentSceneObject)->m_displayCrossSections == CS_NB_OF_CHOICES )
	{
        ((SelectionObject*)m_mainFrame->m_pCurrentSceneObject)->m_displayCrossSections = CS_NOTHING;
	}
}

///////////////////////////////////////////////////////////////////////////
// This function will be triggered when the user click on the display dispersion tube
// button that is located in the m_fibersInfoSizer.
///////////////////////////////////////////////////////////////////////////
void PropertiesWindow::OnDisplayDispersionTube( wxCommandEvent& WXUNUSED(event) )
{
    Logger::getInstance()->print( _T( "Event triggered - PropertiesWindow::OnDisplayDispersionTube" ), LOGLEVEL_DEBUG );

    ((SelectionObject*)m_mainFrame->m_pCurrentSceneObject)->m_displayDispersionCone = (DispersionConeDisplay)( ( (int)((SelectionObject*)m_mainFrame->m_pCurrentSceneObject)->m_displayDispersionCone ) + 1 );
    if( ((SelectionObject*)m_mainFrame->m_pCurrentSceneObject)->m_displayDispersionCone == DC_NB_OF_CHOICES )
	{
        ((SelectionObject*)m_mainFrame->m_pCurrentSceneObject)->m_displayDispersionCone = DC_NOTHING;
	}
}

///////////////////////////////////////////////////////////////////////////
// This function will be called when the rename option is clicked on the right
// click menu of a SelectionObject item in the tree.
///////////////////////////////////////////////////////////////////////////
void PropertiesWindow::OnRenameBox( wxCommandEvent& WXUNUSED(event) )
{
    wxTreeItemId l_treeBoxId = m_mainFrame->m_pTreeWidget->GetSelection();
    
    TreeObjectType objectType = m_mainFrame->treeSelectedNew( l_treeBoxId );
    
    if( objectType == TYPE_SELECTION_OBJECT )
    {
        CustomTreeItem *pTreeItem = (CustomTreeItem*)m_mainFrame->m_pTreeWidget->GetItemData( l_treeBoxId );
        SelectionObject *pSelObject = m_mainFrame->m_pDatasetHelper->m_pSelectionTree->getObject( pTreeItem->getId() );

        wxTextEntryDialog dialog(this, _T( "Please enter a new name" ) );
        dialog.SetValue( pSelObject->getName() );

        if( ( dialog.ShowModal() == wxID_OK ) && ( dialog.GetValue() != _T( "" ) ) )
		{
            pSelObject->setName( dialog.GetValue() );
		}

        m_mainFrame->m_pTreeWidget->SetItemText( l_treeBoxId, pSelObject->getName() );
    }
    
    m_mainFrame->refreshAllGLWidgets();
}

void PropertiesWindow::OnToggleAndNot( wxCommandEvent& WXUNUSED(event) )
{
    Logger::getInstance()->print( _T( "Event triggered - PropertiesWindow::OnToggleAndNot" ), LOGLEVEL_DEBUG );

    if( !m_mainFrame->m_pDatasetHelper->m_theScene)
        return;

    // Get what selection object is selected.
    wxTreeItemId selectionObjectTreeId = m_mainFrame->m_pTreeWidget->GetSelection();

    if( m_mainFrame->treeSelectedNew( selectionObjectTreeId ) == TYPE_SELECTION_OBJECT )
    {
        CustomTreeItem *pTreeItem = (CustomTreeItem*)m_mainFrame->m_pTreeWidget->GetItemData( selectionObjectTreeId );
        SelectionObject *pSelObject = m_mainFrame->m_pDatasetHelper->m_pSelectionTree->getObject( pTreeItem->getId() );
        
        pSelObject->toggleIsNOT();
        
        if( pSelObject->getIsNOT() )
        {
            m_mainFrame->m_pTreeWidget->SetItemBackgroundColour( selectionObjectTreeId, *wxRED   );
        }
        else
        {
            m_mainFrame->m_pTreeWidget->SetItemBackgroundColour( selectionObjectTreeId, *wxGREEN );
        }
        
        wxTreeItemId parentId = m_mainFrame->m_pTreeWidget->GetItemParent( selectionObjectTreeId );
        
        if( m_mainFrame->treeSelectedNew( parentId ) == TYPE_SELECTION_OBJECT )
        {
            CustomTreeItem *pParentTreeItem = (CustomTreeItem*)m_mainFrame->m_pTreeWidget->GetItemData( parentId );
            SelectionObject *pParentSelObject = m_mainFrame->m_pDatasetHelper->m_pSelectionTree->getObject( pParentTreeItem->getId() );
            
            pParentSelObject->setIsDirty( true );
        }
    }
    m_mainFrame->refreshAllGLWidgets();
}


void PropertiesWindow::OnColorRoi( wxCommandEvent& WXUNUSED(event) )
{
    // TODO SELECTION TREE
    if( ! m_mainFrame->m_pDatasetHelper->m_theScene )
        return;

    // Get the currently selected object.
    wxTreeItemId l_selectionObjectTreeId = m_mainFrame->m_pTreeWidget->GetSelection();
    SelectionObject* l_selectionObject = (SelectionObject*)( m_mainFrame->m_pTreeWidget->GetItemData( l_selectionObjectTreeId ) );

    wxColourData l_colorData;

    for( int i = 0; i < 10; ++i )
    {
        wxColour l_color( i * 28, i * 28, i * 28 );
        l_colorData.SetCustomColour( i, l_color );
    }

    int i = 10;
    wxColour color ( 255, 0,   0   );
    wxColour color1( 0,   255, 0   );
    wxColour color2( 0,   0,   255 );
    wxColour color3( 255, 255, 0   );
    wxColour color4( 255, 0,   255 );
    wxColour color5( 0,   255, 255 );

    l_colorData.SetCustomColour( i++, color  );
    l_colorData.SetCustomColour( i++, color1 );
    l_colorData.SetCustomColour( i++, color2 );
    l_colorData.SetCustomColour( i++, color3 );
    l_colorData.SetCustomColour( i++, color4 );
    l_colorData.SetCustomColour( i++, color5 );
#ifdef __WXMAC__
    wxColourDialog dialog( this);
#else
    wxColourDialog dialog( this, &l_colorData );
#endif
    wxColour l_color;

    if( dialog.ShowModal() == wxID_OK )
    {
        wxColourData retData = dialog.GetColourData();
        l_color = retData.GetColour();
    }
    else
        return;

    l_selectionObject->setColor( l_color );

    m_mainFrame->refreshAllGLWidgets();
}

void PropertiesWindow::OnVoiFlipNormals( wxCommandEvent& WXUNUSED(event) )
{
    wxTreeItemId l_selectionObjectTreeId = m_mainFrame->m_pTreeWidget->GetSelection();
    SelectionObject* l_selectionObject = (SelectionObject*)( m_mainFrame->m_pTreeWidget->GetItemData( l_selectionObjectTreeId ) );

    if(l_selectionObject->getSelectionType() == CISO_SURFACE_TYPE)
    {
        l_selectionObject->FlipNormals();        
    }
}

void PropertiesWindow::OnDeleteTreeItem( wxTreeEvent&    event )
{
    m_mainFrame->onDeleteTreeItem(event);
    m_mainFrame->m_pMainGL->m_pRealTimeFibers->clearFibersRTT();
    m_mainFrame->m_pMainGL->m_pRealTimeFibers->clearColorsRTT();
    m_mainFrame->m_pDatasetHelper->m_isRTTDirty = false;
    m_mainFrame->m_pDatasetHelper->m_isRTTReady = false;
    m_mainFrame->m_pTrackingWindow->m_pBtnStart->Enable(false);

}

void PropertiesWindow::OnActivateTreeItem ( wxTreeEvent&    event )
{
    Logger::getInstance()->print( _T( "Event triggered - PropertiesWindow::OnActivateTreeItem" ), LOGLEVEL_DEBUG );

    m_mainFrame->onActivateTreeItem(event);
}

void PropertiesWindow::OnToggleShowSelectionObject( wxCommandEvent& WXUNUSED(event) )
{
    Logger::getInstance()->print( _T( "Event triggered - PropertiesWindow::OnToggleShowSelectionObject" ), LOGLEVEL_DEBUG );

    if( !m_mainFrame->m_pDatasetHelper->m_theScene)
        return;

    // Get the selected selection object.
    wxTreeItemId l_selectionObjectTreeId = m_mainFrame->m_pTreeWidget->GetSelection();

    if( m_mainFrame->treeSelected( l_selectionObjectTreeId ) == MASTER_OBJECT )
    {
        SelectionObject* l_selecitonObject = (SelectionObject*)( m_mainFrame->m_pTreeWidget->GetItemData( l_selectionObjectTreeId ) );
        l_selecitonObject->toggleIsVisible();
        m_mainFrame->m_pTreeWidget->SetItemImage( l_selectionObjectTreeId, l_selecitonObject->getIcon() );
        l_selecitonObject->setIsDirty( true );

        int l_childSelectionObjects = m_mainFrame->m_pTreeWidget->GetChildrenCount( l_selectionObjectTreeId );
        wxTreeItemIdValue childcookie = 0;
        for( int i = 0; i < l_childSelectionObjects; ++i )
        {
            wxTreeItemId l_childId = m_mainFrame->m_pTreeWidget->GetNextChild( l_selectionObjectTreeId, childcookie );
            if( l_childId.IsOk() )
            {
                SelectionObject* childBox = ( (SelectionObject*)( m_mainFrame->m_pTreeWidget->GetItemData( l_childId ) ) );
                childBox->setIsVisible( l_selecitonObject->getIsVisible() );
                m_mainFrame->m_pTreeWidget->SetItemImage( l_childId, childBox->getIcon() );
                childBox->setIsDirty( true );
            }
        }
    }
    else if( m_mainFrame->treeSelected( l_selectionObjectTreeId ) == CHILD_OBJECT )
    {
        SelectionObject *l_selectionObject = (SelectionObject*)( m_mainFrame->m_pTreeWidget->GetItemData( l_selectionObjectTreeId ) );
        l_selectionObject->toggleIsVisible();
        m_mainFrame->m_pTreeWidget->SetItemImage( l_selectionObjectTreeId, l_selectionObject->getIcon() );
        l_selectionObject->setIsDirty( true );
    }

    m_mainFrame->m_pDatasetHelper->m_selBoxChanged = true;
    m_mainFrame->refreshAllGLWidgets();
}

void PropertiesWindow::OnAssignColor( wxCommandEvent& WXUNUSED(event) )
{
    if( ! m_mainFrame->m_pDatasetHelper->m_theScene )
        return;

    wxColourData l_colorData;

    for( int i = 0; i < 10; ++i )
    {
        wxColour l_color(i * 28, i * 28, i * 28);
        l_colorData.SetCustomColour(i, l_color);
    }

    int i = 10;
    wxColour l_color ( 255, 0, 0 );
    l_colorData.SetCustomColour( i++, l_color );
    wxColour l_color1( 0, 255, 0 );
    l_colorData.SetCustomColour( i++, l_color1 );
    wxColour l_color2( 0, 0, 255 );
    l_colorData.SetCustomColour( i++, l_color2 );
    wxColour l_color3( 255, 255, 0 );
    l_colorData.SetCustomColour( i++, l_color3 );
    wxColour l_color4( 255, 0, 255 );
    l_colorData.SetCustomColour( i++, l_color4 );
    wxColour l_color5( 0, 255, 255 );
    l_colorData.SetCustomColour( i++, l_color5 );

    wxColourDialog dialog( this, &l_colorData );
    wxColour l_col;
    if( dialog.ShowModal() == wxID_OK )
    {
        wxColourData l_retData = dialog.GetColourData();
        l_col = l_retData.GetColour();
    }
    else
    {
        return;
    }

    if( m_mainFrame->m_currentListItem != -1 )
    {
        DatasetInfo *l_info = (DatasetInfo*)m_mainFrame->m_pCurrentSceneObject;
        if( l_info->getType() == MESH || l_info->getType() == ISO_SURFACE || l_info->getType() == SURFACE || l_info->getType() == VECTORS)
        {
            l_info->setColor( l_col );
            l_info->setuseTex( false );
            m_mainFrame->m_pListCtrl->SetItem( m_mainFrame->m_currentListItem, 2, wxT( "(") + wxString::Format( wxT( "%.2f" ), l_info->getThreshold() ) + wxT( ")" ) );           
        }
    }
    else if ( m_mainFrame->m_pDatasetHelper->m_lastSelectedObject != NULL )
    {
        SelectionObject *l_selObj = (SelectionObject*)m_mainFrame->m_pCurrentSceneObject;
        if (!l_selObj->getIsMaster())
        {
            wxTreeItemId l_parentId = m_mainFrame->m_pTreeWidget->GetItemParent( m_mainFrame->m_pDatasetHelper->m_lastSelectedObject->GetId());
            l_selObj = (SelectionObject*)m_mainFrame->m_pTreeWidget->GetItemData(l_parentId);
        }
        l_selObj->setFiberColor( l_col);
        l_selObj->setIsDirty( true );
        m_mainFrame->m_pDatasetHelper->m_selBoxChanged = true;
    }    
    m_mainFrame->refreshAllGLWidgets();
}

void PropertiesWindow::OnCreateFibersDensityTexture( wxCommandEvent& WXUNUSED(event) )
{
    Fibers* l_fibers = NULL;

    if( ! m_mainFrame->m_pDatasetHelper->getSelectedFiberDataset(l_fibers) )
        return ;

    int l_x,l_y,l_z;

    Anatomy* l_newAnatomy = new Anatomy( m_mainFrame->m_pDatasetHelper );
    l_newAnatomy->setZero( m_mainFrame->m_pDatasetHelper->m_columns, m_mainFrame->m_pDatasetHelper->m_rows, m_mainFrame->m_pDatasetHelper->m_frames );
    l_newAnatomy->setDataType( 16 );
    l_newAnatomy->setType( OVERLAY );
    float l_max = 0.0f;
    wxTreeItemId l_treeObjectId = m_mainFrame->m_pTreeWidget->GetSelection();
    if( m_mainFrame->treeSelected( l_treeObjectId ) == MASTER_OBJECT )
    {
        SelectionObject* l_object = (SelectionObject*)( m_mainFrame->m_pTreeWidget->GetItemData( l_treeObjectId ) );

        std::vector<float>* l_dataset = l_newAnatomy->getFloatDataset();

        for( int l = 0; l < l_fibers->getLineCount(); ++l )
        {
            if( l_object->m_inBranch[l] )
            {
                unsigned int pc = l_fibers->getStartIndexForLine(l)*3;

                for( int j = 0; j < l_fibers->getPointsPerLine(l) ; ++j )
                {
                    l_x = (int)( l_fibers->getPointValue(pc) / m_mainFrame->m_pDatasetHelper->m_xVoxel );
                    ++pc;
                    l_y = (int)( l_fibers->getPointValue(pc) / m_mainFrame->m_pDatasetHelper->m_yVoxel );
                    ++pc;
                    l_z = (int)( l_fibers->getPointValue(pc) / m_mainFrame->m_pDatasetHelper->m_zVoxel );
                    ++pc;

                    int index =( l_x + l_y * m_mainFrame->m_pDatasetHelper->m_columns + l_z * m_mainFrame->m_pDatasetHelper->m_columns * m_mainFrame->m_pDatasetHelper->m_rows );
                    l_dataset->at(index) += 1.0;
                    l_max = wxMax( l_max,l_dataset->at(index) );
                }
            }
        }
        for( int i = 0 ; i < m_mainFrame->m_pDatasetHelper->m_columns * m_mainFrame->m_pDatasetHelper->m_rows * m_mainFrame->m_pDatasetHelper->m_frames ; ++i )
        {
            l_dataset->at(i) /= l_max;
        }
    }

    l_newAnatomy->setName( wxT(" (fiber_density)" ) );
    l_newAnatomy->setOldMax( l_max );
    m_mainFrame->m_pListCtrl->InsertItem( 0, wxT( "" ), 0 );
    m_mainFrame->m_pListCtrl->SetItem( 0, 1, l_newAnatomy->getName() );
    m_mainFrame->m_pListCtrl->SetItem( 0, 2, wxT( "0.00" ) );
    m_mainFrame->m_pListCtrl->SetItem( 0, 3, wxT( "" ), 1 );
    m_mainFrame->m_pListCtrl->SetItemData( 0, (long) l_newAnatomy );
    m_mainFrame->m_pListCtrl->SetItemState( 0, wxLIST_STATE_SELECTED, wxLIST_STATE_SELECTED );
    m_mainFrame->refreshAllGLWidgets();
}

void PropertiesWindow::OnCreateFibersColorTexture( wxCommandEvent& WXUNUSED(event) )
{
    Fibers* l_fibers = NULL;

    if( !m_mainFrame->m_pDatasetHelper->getSelectedFiberDataset(l_fibers) )
        return ;

    int l_x,l_y,l_z;
    Anatomy* l_newAnatomy = new Anatomy( m_mainFrame->m_pDatasetHelper );
    l_newAnatomy->setRGBZero( m_mainFrame->m_pDatasetHelper->m_columns, m_mainFrame->m_pDatasetHelper->m_rows, m_mainFrame->m_pDatasetHelper->m_frames );

    wxTreeItemId l_treeObjectId = m_mainFrame->m_pTreeWidget->GetSelection();
    if(m_mainFrame-> treeSelected( l_treeObjectId ) == MASTER_OBJECT )
    {
        SelectionObject* l_object = (SelectionObject*)( m_mainFrame->m_pTreeWidget->GetItemData( l_treeObjectId ) );
        wxColour l_color = l_object->getFiberColor();

        std::vector<float>* l_dataset = l_newAnatomy->getFloatDataset();

        for( int l = 0; l < l_fibers->getLineCount(); ++l )
        {
            if( l_object->m_inBranch[l] )
            {
                unsigned int pc = l_fibers->getStartIndexForLine( l ) * 3;

                for( int j = 0; j < l_fibers->getPointsPerLine( l ) ; ++j )
                {
                    l_x = (int)( l_fibers->getPointValue( pc ) / m_mainFrame->m_pDatasetHelper->m_xVoxel );
                    ++pc;
                    l_y = (int)( l_fibers->getPointValue( pc ) / m_mainFrame->m_pDatasetHelper->m_yVoxel );
                    ++pc;
                    l_z = (int)( l_fibers->getPointValue( pc ) / m_mainFrame->m_pDatasetHelper->m_zVoxel );
                    ++pc;

                    int index = ( l_x + l_y * m_mainFrame->m_pDatasetHelper->m_columns + l_z *m_mainFrame-> m_pDatasetHelper->m_columns * m_mainFrame->m_pDatasetHelper->m_rows ) * 3;
                    l_dataset->at( index )     = l_color.Red()   / 255.0f;
                    l_dataset->at( index + 1 ) = l_color.Green() / 255.0f;
                    l_dataset->at( index + 2 ) = l_color.Blue()  / 255.0f;
                }
            }
        }
    }

    l_newAnatomy->setName( wxT( " (fiber_colors)" ) );
    m_mainFrame->m_pListCtrl->InsertItem( 0, wxT( "" ), 0 );
    m_mainFrame->m_pListCtrl->SetItem( 0, 1, l_newAnatomy->getName() );
    m_mainFrame->m_pListCtrl->SetItem( 0, 2, wxT( "0.00" ) );
    m_mainFrame->m_pListCtrl->SetItem( 0, 3, wxT( "" ), 1 );
    m_mainFrame->m_pListCtrl->SetItemData( 0, (long)l_newAnatomy );
    m_mainFrame->m_pListCtrl->SetItemState( 0, wxLIST_STATE_SELECTED, wxLIST_STATE_SELECTED );
    m_mainFrame->refreshAllGLWidgets();
}

void PropertiesWindow::OnMeanComboBoxSelectionChange( wxCommandEvent& event)
{
    Logger::getInstance()->print( _T( "Event triggered - PropertiesWindow::OnMeanComboBoxSelectionChange" ), LOGLEVEL_DEBUG );

    ((SelectionObject*)m_mainFrame->m_pCurrentSceneObject)->SetFiberInfoGridValues();
    m_mainFrame->refreshAllGLWidgets();
}

void PropertiesWindow::OnBoxPositionX( wxCommandEvent &event )
{    
    double posX = 0;
    Vector currPos;
    ((SelectionObject*)m_mainFrame->m_pCurrentSceneObject)->m_ctrlBoxX->GetValue().ToDouble(&posX);  
    currPos = ((SelectionObject*)m_mainFrame->m_pCurrentSceneObject)->getCenter();
    ((SelectionObject*)m_mainFrame->m_pCurrentSceneObject)->setCenter(posX,currPos.y,currPos.z);
}

void PropertiesWindow::OnBoxPositionY( wxCommandEvent &event )
{    
    double posY = 0;
    Vector currPos;
    ((SelectionObject*)m_mainFrame->m_pCurrentSceneObject)->m_ctrlBoxY->GetValue().ToDouble(&posY);  
    currPos = ((SelectionObject*)m_mainFrame->m_pCurrentSceneObject)->getCenter();
    ((SelectionObject*)m_mainFrame->m_pCurrentSceneObject)->setCenter(currPos.x,posY,currPos.z);
}

void PropertiesWindow::OnBoxPositionZ( wxCommandEvent &event )
{    
    double posZ = 0;
    Vector currPos;
    ((SelectionObject*)m_mainFrame->m_pCurrentSceneObject)->m_ctrlBoxZ->GetValue().ToDouble(&posZ);  
    currPos = ((SelectionObject*)m_mainFrame->m_pCurrentSceneObject)->getCenter();
    ((SelectionObject*)m_mainFrame->m_pCurrentSceneObject)->setCenter(currPos.x,currPos.y,posZ);
}

void PropertiesWindow::OnBoxSizeX( wxCommandEvent &event )
{    
    double sizeX = 0;
    Vector currSize;
    ((SelectionObject*)m_mainFrame->m_pCurrentSceneObject)->m_ctrlBoxSizeX->GetValue().ToDouble(&sizeX);  
    currSize = ((SelectionObject*)m_mainFrame->m_pCurrentSceneObject)->getSize();
    currSize.x = sizeX/m_mainFrame->m_pDatasetHelper->m_xVoxel;
    ((SelectionObject*)m_mainFrame->m_pCurrentSceneObject)->setSize(currSize);
}

void PropertiesWindow::OnBoxSizeY( wxCommandEvent &event )
{    
    double sizeY = 0;
    Vector currSize;
    ((SelectionObject*)m_mainFrame->m_pCurrentSceneObject)->m_ctrlBoxSizeY->GetValue().ToDouble(&sizeY);  
    currSize = ((SelectionObject*)m_mainFrame->m_pCurrentSceneObject)->getSize();
    currSize.y = sizeY/m_mainFrame->m_pDatasetHelper->m_yVoxel;
    ((SelectionObject*)m_mainFrame->m_pCurrentSceneObject)->setSize(currSize);
}

void PropertiesWindow::OnBoxSizeZ( wxCommandEvent &event )
{    
    double sizeZ = 0;
    Vector currSize;
    ((SelectionObject*)m_mainFrame->m_pCurrentSceneObject)->m_ctrlBoxSizeZ->GetValue().ToDouble(&sizeZ);  
    currSize = ((SelectionObject*)m_mainFrame->m_pCurrentSceneObject)->getSize();
    currSize.z = sizeZ/m_mainFrame->m_pDatasetHelper->m_zVoxel;
    ((SelectionObject*)m_mainFrame->m_pCurrentSceneObject)->setSize(currSize);
}

void PropertiesWindow::OnSliderAxisMoved( wxCommandEvent& WXUNUSED(event) )
{
    float l_sliderValue = ((ODFs*)m_mainFrame->m_pCurrentSceneObject)->m_pSliderFlood->GetValue() / 10.0f;
    ((ODFs*)m_mainFrame->m_pCurrentSceneObject)->m_axisThreshold = l_sliderValue;
    ((ODFs*)m_mainFrame->m_pCurrentSceneObject)->m_pTxtThresBox->SetValue(wxString::Format( wxT( "%.1f"), l_sliderValue));

    std::cout << ((ODFs*)m_mainFrame->m_pCurrentSceneObject)->m_axisThreshold << std::endl;
}

void PropertiesWindow::OnRecalcMainDir( wxCommandEvent& WXUNUSED(event) )
{   
    ((ODFs*)m_mainFrame->m_pCurrentSceneObject)->extractMaximas();

    /*for( int z = 0; z < m_mainFrame->m_datasetHelper->m_frames; z++ )
        for( int y = 0; y < m_mainFrame->m_datasetHelper->m_rows; y++ )
            for( int x = 0; x < m_mainFrame->m_datasetHelper->m_columns; x++ )
            {
                int  currentIdx = ((Glyph*)m_mainFrame->m_currentSceneObject)->getGlyphIndex( z, y, x );

                if(((ODFs*)m_mainFrame->m_currentSceneObject)->getCoeffs().at(currentIdx)[0] != 0)
                    ((ODFs*)m_mainFrame->m_currentSceneObject)->mainDirections[currentIdx] = ((ODFs*)m_mainFrame->m_currentSceneObject)->getODFmaxNotNorm(((ODFs*)m_mainFrame->m_currentSceneObject)->getCoeffs().at(currentIdx),
                    ((ODFs*)m_mainFrame->m_currentSceneObject)->getShMatrix()[NB_OF_LOD - 1], ((ODFs*)m_mainFrame->m_currentSceneObject)->getPhiTheta()[NB_OF_LOD - 1],((ODFs*)m_mainFrame->m_currentSceneObject)->m_axisThreshold,((ODFs*)m_mainFrame->m_currentSceneObject)->angle,((ODFs*)m_mainFrame->m_currentSceneObject)->m_nbors);
            }*/
}

void PropertiesWindow::OnToggleCrossingFibers( wxEvent& WXUNUSED(event) )
{
    if( m_mainFrame->m_pCurrentSceneObject != NULL && m_mainFrame->m_currentListItem != -1 )
    {
        if( ((DatasetInfo*)m_mainFrame->m_pCurrentSceneObject)->getType() == FIBERS )
        {
            ((Fibers*)m_mainFrame->m_pCurrentSceneObject)->toggleCrossingFibers();
        }
    }
}

void PropertiesWindow::OnCrossingFibersThicknessChange( wxCommandEvent& WXUNUSED(event) )
{
    if( m_mainFrame->m_pCurrentSceneObject != NULL && m_mainFrame->m_currentListItem != -1 )
    {
        if( ((DatasetInfo*)m_mainFrame->m_pCurrentSceneObject)->getType() == FIBERS )
        {
            ((Fibers*)m_mainFrame->m_pCurrentSceneObject)->updateCrossingFibersThickness();
        }
    }
}

void PropertiesWindow::AddSelectionObjectToSelectionTree( SelectionObject *pSelObj,
                                                          const wxTreeItemId &parentTreeId )
{
    wxTreeItemId newSelectionObjectId;
    
    if( m_mainFrame->m_pDatasetHelper->m_pSelectionTree->isEmpty() )
    {
        pSelObj->setIsMaster( true );
        int rootId = m_mainFrame->m_pDatasetHelper->m_pSelectionTree->setRoot( pSelObj );
        
        CustomTreeItem *pTreeItem = new CustomTreeItem( rootId );
        newSelectionObjectId = m_mainFrame->m_pTreeWidget->AppendItem( m_mainFrame->m_tSelectionObjectsId, pSelObj->getName(), 0, -1, pTreeItem );
        
        m_mainFrame->m_pTreeWidget->SetItemBackgroundColour( newSelectionObjectId, *wxCYAN );
    }
    else
    {
        CustomTreeItem *pItem = (CustomTreeItem*) m_mainFrame->m_pTreeWidget->GetItemData( parentTreeId );
        
        int parentId = pItem->getId();
        
        int childId = m_mainFrame->m_pDatasetHelper->m_pSelectionTree->addChildrenObject( parentId,  pSelObj );
        
        CustomTreeItem *pTreeItem = new CustomTreeItem( childId );
        newSelectionObjectId = m_mainFrame->m_pTreeWidget->AppendItem( parentTreeId, pSelObj->getName(), 0, -1, pTreeItem );
        
        m_mainFrame->m_pTreeWidget->SetItemBackgroundColour( newSelectionObjectId, *wxGREEN );
    }
    
    pSelObj->setTreeId( newSelectionObjectId );
    m_mainFrame->m_pTreeWidget->EnsureVisible( newSelectionObjectId );
    m_mainFrame->m_pTreeWidget->SetItemImage( newSelectionObjectId, pSelObj->getIcon() );
    m_mainFrame->m_pTreeWidget->SelectItem( newSelectionObjectId, true );
    
    m_mainFrame->m_pDatasetHelper->m_selBoxChanged = true;
}

void PropertiesWindow::AddSelectionObjectsToSelectionTree( const vector< SelectionObject* > &selObjects, bool addAsChildOfFirst /* = false */ )
{
    if( selObjects.empty() )
    {
        return;
    }
    
    wxTreeItemId baseSelectionId = m_mainFrame->m_pTreeWidget->GetSelection();
    
    AddSelectionObjectToSelectionTree( selObjects[0], baseSelectionId );
    
    wxTreeItemId parentSelectionId = m_mainFrame->m_pTreeWidget->GetSelection();
    
    for( unsigned int selObjIdx( 1 ); selObjIdx < selObjects.size(); ++selObjIdx )
    {
        AddSelectionObjectToSelectionTree( selObjects[ selObjIdx ], parentSelectionId );
        
        if( !addAsChildOfFirst )
        {
            parentSelectionId = m_mainFrame->m_pTreeWidget->GetSelection();
        }
    }
    
    return;
}

