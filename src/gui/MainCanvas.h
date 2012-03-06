#ifndef MAINCANVAS_H_
#define MAINCANVAS_H_

#include "wx/wxprec.h"

#ifndef WX_PRECOMP
#include "wx/wx.h"
#endif

#include "../gfx/TheScene.h"
#include "wx/glcanvas.h"

#include "BoundingBox.h"
#include "../dataset/DatasetHelper.h"
#include "wx/math.h"
#include "../misc/Fantom/FArray.h"
#include <vector>
#include <deque>
#include <time.h>
#include "stdlib.h"
#include "../dataset/RTTFibers.h"

#define SWAP(a,b,t) t tempr=(a);(a)=(b);(b)=tempr
#define SQR(a) ((a)*(a))

class DatasetHelper;
class RTTFibers;

class MainCanvas: public wxGLCanvas
{
public:

    MainCanvas(DatasetHelper*, int,  wxWindow*, wxWindowID id = wxID_ANY,
        const wxPoint& pos = wxDefaultPosition,
        const wxSize& size = wxDefaultSize,
        long style = 0, const wxString& name = _T("GLCanvas"),
#ifdef CTX
        int* gl_attrib = NULL, wxGLContext* ctx = NULL);
#else
        int* gl_attrib = NULL,
        wxGLCanvas* shared=( wxGLCanvas* )NULL);
#endif
   ~MainCanvas();
    
    void OnPaint(wxPaintEvent& event);
    void OnSize(wxSizeEvent& event);
    void OnShow(wxShowEvent& event);
    void changeOrthoSize();
    void OnEraseBackground(wxEraseEvent& event);
    void OnMouseEvent(wxMouseEvent& event);

	wxPoint getMousePos() {return m_clicked;};
    float getDelta() {return m_delta;};
    int getPicked() {return m_hr.picked;};
    void updateView();
    void init();
    void render();
    void invalidate();
    void renderTestRay();
    void renderAxes();
	void renderRulerDisplay();
	void drawOnAnatomy();
    void segment();
    void KMeans(float i_means[2],float i_stddev[2],float i_apriori[2],std::vector<float>*,std::vector<float>*);
    void floodFill(std::vector<float>*, std::vector<float>*, Vector, float);
    void graphCut(std::vector<float>*, std::vector<float>*, float);
    float getElement(int,int,int,std::vector<float>*);

	void pushAnatomyHistory();
	void popAnatomyHistory();

    hitResult pick(wxPoint, bool i_isRulerOrDrawer);
    float getAxisParallelMovement(int, int, int, int, Vector);
    Vector getEventCenter();
    void setRotation();
    void OnChar(wxKeyEvent& event);

    std::vector<Vector> object;
    std::vector<Vector> background;
    bool m_isSlizesLocked;
    bool m_isSceneLocked;
    bool m_init;
    RTTFibers *m_pRealTimeFibers; 




    DECLARE_EVENT_TABLE()

private:
     int         m_view;           // identifier of gl widget
     wxPoint     m_clicked;        // last clicked point
     wxPoint     m_lastPos;        // last position to calculate right mouse dragging     
     float       m_delta;
     hitResult   m_hr;
     Vector      m_hitPts;
	 bool        m_isRulerHit;
	 bool        m_isDrawerHit;
    
     DatasetHelper *m_pDatasetHelper; 
     

     GLdouble m_pos1X, m_pos1Y, m_pos1Z;    // point of origin for picking
     GLdouble m_pos2X, m_pos2Y, m_pos2Z;    // distant point for picking


     //Matrix4fT m_transform;        // transformation matrix of current view
     Matrix3fT m_lastRot;
     Matrix3fT m_thisRot ;

     ArcBallT   *m_pArcBall;
     Point2fT   m_mousePt;

     int m_orthoSizeNormal;
     int m_orthoModX;
     int m_orthoModY;

     GLint m_viewport[4];                   // view context for picking
     GLdouble m_projection[16];             // view context for picking
     GLdouble m_modelview[16];              // view context for picking
     

};

#endif /*MAINCANVAS_H_*/
