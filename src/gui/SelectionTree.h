/*
 *  The SelectionTree and SelectionTreeNode classes declaration.
 *
 */

#ifndef SELECTIONTREE_H_
#define SELECTIONTREE_H_

#include "SelectionObject.h"

#include <vector>
using std::vector;

class Fibers;
class Octree;

class SelectionTree
{
public:
    typedef vector< SelectionObject* > SelectionObjectVector;
    
public:
    SelectionTree();
    ~SelectionTree();
    
    int setRoot( SelectionObject *pRootSelObject );
    
    int addChildrenObject( const int parentId, SelectionObject *pSelObject );
    
    bool removeObject( const int nodeId );
    
    SelectionObject* getObject( const int itemId ) const;
    
    SelectionObject* getParentObject( SelectionObject *pSelObj ) const;
    
    SelectionObjectVector getAllObjects() const;
    
    SelectionObjectVector getChildrenObjects( const int itemId ) const;
    
    int getActiveChildrenObjectsCount( SelectionObject *pSelObject ) const;
    
    bool isEmpty() const
    {
        return m_pRootNode == NULL;
    }
    
    bool containsId( const int itemId ) const;
    
    bool isRootObject( SelectionObject *pSelObj ) const
    {
        if( m_pRootNode != NULL )
        {
            return m_pRootNode->getSelectionObject() == pSelObj;
        }
        
        return false;
    }
    
    bool isRootItem( const int itemId ) const
    {
        if( m_pRootNode != NULL )
        {
            return m_pRootNode->getId() == itemId;
        }
        return false;
    }
    
    // Methods related to fiber selection.
    vector< bool > getSelectedFibers( const Fibers* const pFibers );
    
private:
    class SelectionTreeNode
    {
    public:
        SelectionTreeNode( const int id, SelectionObject *pSelObject );
        ~SelectionTreeNode();
        
        void setSelectionObject( SelectionObject *pSelObject );
        SelectionObject* getSelectionObject() const;
        
        SelectionObjectVector getAllSelectionObjects() const;
        SelectionObjectVector getAllChildrenSelectionObjects() const;
        
        int getActiveDirectChildrenCount() const;

        void addChildren( SelectionTreeNode *pNode );
        bool removeChildren( const int nodeId );
        void removeAllChildren();
        
        SelectionTreeNode * const findNode( const int nodeId );
        SelectionTreeNode * const findParentNode( const int searchedChildNodeId );
        
        SelectionTreeNode * const findNode( SelectionObject *pSelObj );
        
        void updateInObjectRecur( const int fibersCount, Octree *pCurOctree, const vector< int > &reverseIdx );
        void updateInBranchRecur( const int fibersCount );
        
        int getId() const;
        
    private:
        SelectionTreeNode();    // Disable default constructor.
        
    private:
        int m_nodeId;
        SelectionObject *m_pSelObject;
        vector< SelectionTreeNode* > m_children;
    };
    
    // Structure used with the find_if algorithm.
    struct SelectionTreeNodeFinder
    {
    public:
        SelectionTreeNodeFinder( const int searchedNodeId )
        : m_searchedId( searchedNodeId )
        {}
        
        bool operator() ( const SelectionTreeNode *treeNode )
        {
            return treeNode->getId() == m_searchedId;
        }
        
    private:
        int m_searchedId;
    };
    
private:
    SelectionTreeNode *m_pRootNode;
    
    int m_nextNodeId;
};



#endif // SELECTIONTREE_H_