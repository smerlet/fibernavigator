#ifndef THESCENE_H_
#define THESCENE_H_

#include "wx/wxprec.h"

#ifndef WX_PRECOMP
#include "wx/wx.h"
#endif

#include "theDataset.h"

#include "GLSL/GLSLShaderProgram.h"
#include "ArcBall.h"
#include "wx/glcanvas.h"

enum {
	axial,
	coronal,
	sagittal,
	mainView
};

class TheScene {

public:
	bool m_showXSlize;
	bool m_showYSlize;
	bool m_showZSlize;
	bool m_showMesh;
	
	TheScene();
	~TheScene();
	
	void initGL(int);
	
	void assignTextures();
	void addTexture();
	void swapTextures(int, int);
	void releaseTextures();
	GLuint makeCallList(DatasetInfo*);
	
	void initShaders();
	
	void setDataset(TheDataset*);
	void setDataListCtrl(wxListCtrl* value) {m_listctrl = value;};
	void setMainGLContext(wxGLContext* context) {m_mainGLContext = context;};
	wxGLContext* getMainGLContext() {return m_mainGLContext;};
	
	void renderScene(int, int);
	void renderNavView(int);
	void makeLights();
	void renderMesh();
	void colorMap(float);
	TheDataset* getDataset() {return m_dataset;};
	
	void updateView(float, float, float);
	
	
	bool m_texAssigned;
	
		
		
private:
	int m_countTextures;
	GLuint *m_texNames;
	FGLSLShaderProgram *m_textureShader; 
	FGLSLShaderProgram *m_meshShader;
	
	TheDataset* m_dataset;
	wxListCtrl* m_listctrl;
	wxGLContext* m_mainGLContext;
	
	float m_xSize;
	float m_ySize;
	float m_zSize;
	
	float m_xOffset0;
	float m_yOffset0;
	float m_xOffset1;
	float m_yOffset1;
	float m_xOffset2;
	float m_yOffset2;
	float m_xSlize;
	float m_ySlize;
	float m_zSlize;

	float m_ratio0;
	float m_ratio1;
	float m_ratio2;

	int m_quadrant;
	
	void bindTextures();
	void setTextureShaderVars();
	void setMeshShaderVars();
	void renderXSlize();
	void renderYSlize();
	void renderZSlize();
};

#endif /*THESCENE_H_*/
