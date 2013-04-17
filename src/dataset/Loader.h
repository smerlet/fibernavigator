#ifndef LOADER_H_
#define LOADER_H_

#include "DatasetInfo.h"
#include "DatasetManager.h"
#include "../Logger.h"
#include "../gui/ListCtrl.h"
#include "../gui/MainFrame.h"
#include "../gui/SceneManager.h"

#include <wx/file.h>
#include <wx/string.h>

class Loader
{
private:
    MainFrame *m_pMainFrame;
    ListCtrl *m_pListCtrl;
    unsigned int m_error;
    
    // Set this to true to force the loader to load
    // anats looking as r4 ODFs as a file with 5 peaks per voxel instead.
    bool m_loadAsPeaks;
    
    // Set this to true if loading EAPs.
    // Needs to have a special case since the number of coefficients
    // can differ wildly.
    bool m_loadEAPs;
public:
    Loader( MainFrame *pMainFrame, ListCtrl *pListCtrl, const bool loadAsPeaks = false,
            const bool loadAsEAPs = false ) 
    :   m_pMainFrame( pMainFrame ),
        m_pListCtrl( pListCtrl ),
        m_error( 0 ),
        m_loadAsPeaks( loadAsPeaks ),
        m_loadEAPs( loadAsEAPs )
    {
    }

    unsigned int getNbErrors() const { return m_error; }

    void operator()( const wxString &filename )
    {
        // check if i_fileName is valid
        if( !wxFile::Exists( filename ) )
        {
            Logger::getInstance()->print( wxString::Format( wxT( "File \"%s\" doesn't exist!" ), filename.c_str() ), LOGLEVEL_ERROR );
            ++m_error;
        }
        else
        {
            // If the file is in compressed formed, we check what kind of file it is
            wxString extension = filename.AfterLast( '.' );

            #ifdef __WXMSW__
            char separator = '\\';
            #else
            char separator = '/';
            #endif

            wxString name = filename.AfterLast( separator );

            if( wxT( "gz" ) == extension )
            {
                extension = filename.BeforeLast( '.' ).AfterLast( '.' );
            }

            if( wxT( "scn" ) == extension )
            {
                if( !SceneManager::getInstance()->load( filename ) )
                {
                    ++m_error;
                }

                m_pMainFrame->GetStatusBar()->SetStatusText( wxT( "Ready" ), 1 );
                m_pMainFrame->GetStatusBar()->SetStatusText( wxString::Format( wxT( "%s loaded" ), name.c_str() ), 2 );
            }
            else
            {
                DatasetIndex result;
                
                if( m_loadEAPs )
                {
                    result = DatasetManager::getInstance()->loadEAPs( filename, extension );
                }
                else
                {
                    DatasetManager::getInstance()->forceLoadingAsMaximas( m_loadAsPeaks );
                
                    result = DatasetManager::getInstance()->load( filename, extension );
                
                    DatasetManager::getInstance()->forceLoadingAsMaximas( false );
                }
                
                if( result.isOk() )
                {
                    DatasetInfo *pDataset = DatasetManager::getInstance()->getDataset( result );

                    switch( pDataset->getType() )
                    {
                        case HEAD_BYTE:
                        case HEAD_SHORT:
                        case OVERLAY:
                        case RGB:
                        {
                            if( 1 == DatasetManager::getInstance()->getAnatomyCount() )
                            {
                                m_pMainFrame->updateSliders();
                            }
                            break;
                        }
                        case FIBERS:
                        {
                            if( !DatasetManager::getInstance()->isFibersGroupLoaded() )
                            {
                                DatasetIndex result = DatasetManager::getInstance()->createFibersGroup();
                                m_pListCtrl->InsertItem( result );
                            }
                            break;
                        }
                        default:
                            break;
                    }

                    m_pListCtrl->InsertItem( result );

                    m_pMainFrame->GetStatusBar()->SetStatusText( wxT( "Ready" ), 1 );
                    m_pMainFrame->GetStatusBar()->SetStatusText( wxString::Format( wxT( "%s loaded" ), name.c_str() ), 2 );
                }
                else
                {
                    ++m_error;
                }
            }
        }
    }
};

#endif // LOADER_H_
