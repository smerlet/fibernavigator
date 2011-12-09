/*
 *  The SelectionTree and SelectionTreeNode classes implementations.
 *
 */

#include "SelectionTree.h"
#include "SelectionObject.h"

#include "../dataset/Fibers.h"
#include "../dataset/Octree.h"

#include <algorithm>

// Anonymous namespace
namespace {
    void combineBoolVectors( vector< bool > &dest, const vector< bool > &source, const bool negateSrc, const bool useAnd )
    {
        if( dest.size() != source.size() )
        {
            return;
        }

        for( unsigned int elemIdx( 0 ); elemIdx < dest.size(); ++elemIdx )
        {
            if( !negateSrc && !useAnd )
            {
                dest[elemIdx] = dest[elemIdx] | source[elemIdx];
            }
            else if( !negateSrc && useAnd )
            {
                dest[elemIdx] = dest[elemIdx] & source[elemIdx];
            }
            else if( negateSrc && useAnd )
            {
                dest[elemIdx] = dest[elemIdx] & !source[elemIdx];
            }
            else // Negate and use Or
            {
                dest[elemIdx] = dest[elemIdx] | !source[elemIdx];
            }
        }
    }
};

/////
// SelectionTreeNode methods
/////

SelectionTree::SelectionTreeNode::SelectionTreeNode( const int id, SelectionObject *pSelObject )
    : m_nodeId( id ),
      m_pSelObject( pSelObject )
{}

void SelectionTree::SelectionTreeNode::setSelectionObject( SelectionObject *pSelObject )
{
    m_pSelObject = pSelObject;
}

SelectionObject* SelectionTree::SelectionTreeNode::getSelectionObject() const
{
    return m_pSelObject;
}

SelectionObjectVector SelectionTree::SelectionTreeNode::getAllSelectionObjects() const
{
    SelectionObjectVector objs;
    
    objs.push_back( m_pSelObject );
    
    SelectionObjectVector childObjs = getAllChildrenSelectionObjects();
    objs.insert( objs.end(), childObjs.begin(), childObjs.end() );
    
    return objs;
}

SelectionObjectVector SelectionTree::SelectionTreeNode::getAllChildrenSelectionObjects() const
{
    SelectionObjectVector objs;
    
    if( !m_children.empty() )
    {
        for( unsigned int childIdx( 0 ); childIdx < m_children.size(); ++childIdx )
        {
            SelectionObjectVector childObjs = m_children[childIdx]->getAllSelectionObjects();
            
            objs.insert( objs.end(), childObjs.begin(), childObjs.end() );
        }
    }
    
    return objs;
}

int SelectionTree::SelectionTreeNode::getActiveDirectChildrenCount() const
{
    int count( 0 );
    
    for( unsigned int childIdx( 0 ); childIdx < m_children.size(); ++childIdx )
    {
        if( m_children[ childIdx ]->m_pSelObject->getIsActive() )
        {
            ++count;
        }
    }
    
    return count;
}

void SelectionTree::SelectionTreeNode::addChildren( SelectionTreeNode *pNode )
{
    m_children.push_back( pNode );
}

bool SelectionTree::SelectionTreeNode::removeChildren( const int nodeId )
{
    vector< SelectionTreeNode* >::iterator foundPos = std::find_if( m_children.begin(), 
                                                                    m_children.end(),
                                                                   SelectionTreeNodeFinder( nodeId ) );
    
    if( foundPos != m_children.end() )
    {
        delete *foundPos;
        *foundPos = NULL;
        
        m_children.erase( foundPos );
        
        // TODO update selection
        
        return true;
    }
    
    return false;
}

void SelectionTree::SelectionTreeNode::removeAllChildren()
{
    for( vector< SelectionTreeNode* >::iterator nodeIt( m_children.begin() );
         nodeIt != m_children.end(); 
         ++nodeIt )
    {
        delete *nodeIt;
        *nodeIt = NULL;
    }
    
    m_children.clear();
}

SelectionTree::SelectionTreeNode * const
    SelectionTree::SelectionTreeNode::findNode( const int nodeId )
{
    if( getId() == nodeId )
    {
        return this;
    }
    
    for( unsigned int nodeIdx(0); nodeIdx < m_children.size(); ++nodeIdx )
    {
        SelectionTreeNode * const pReturnedNode = m_children[nodeIdx]->findNode( nodeId );
        
        if( pReturnedNode != NULL )
        {
            return pReturnedNode;
        }
    }
    
    return NULL;
}

SelectionTree::SelectionTreeNode * const
    SelectionTree::SelectionTreeNode::findParentNode( const int searchedChildNodeId )
{
    for( unsigned int childNodeIdx( 0 ); childNodeIdx < m_children.size(); ++childNodeIdx )
    {
        if( m_children[childNodeIdx]->getId() == searchedChildNodeId )
        {
            return this;
        }
        else
        {
            SelectionTreeNode *pFoundNode = m_children[childNodeIdx]->findParentNode( searchedChildNodeId );
            
            if( pFoundNode != NULL )
            {
                return pFoundNode;
            }
        }
    }
    
    return NULL;
}

SelectionTree::SelectionTreeNode * const
    SelectionTree::SelectionTreeNode::findNode( SelectionObject *pSelObj )
{
    if( m_pSelObject == pSelObj )
    {
        return this;
    }
    
    for( unsigned int nodeIdx(0); nodeIdx < m_children.size(); ++nodeIdx )
    {
        SelectionTreeNode * const pReturnedNode = m_children[nodeIdx]->findNode( pSelObj );
        
        if( pReturnedNode != NULL )
        {
            return pReturnedNode;
        }
    }
    
    return NULL;
}

void SelectionTree::SelectionTreeNode::updateInObjectRecur( const int fibersCount, Octree *pCurOctree, const vector< int > &reverseIdx )
{
    // TODO check if dirty, therefore if it is needed.
    
    vector< int > pointsInsideObject = pCurOctree->getPointsInside( m_pSelObject );
    
    //m_pSelObject->m_inBox.resize( fibersCount, false );
    m_pSelObject->m_inBox.assign( fibersCount, false );
    
    for( unsigned int ptIdx( 0 ); ptIdx < pointsInsideObject.size(); ++ptIdx )
    {
        m_pSelObject->m_inBox[ reverseIdx[ pointsInsideObject[ ptIdx ] ] ] = true;
    }
    
    // Call this recursively for all children.
    for( unsigned int childIdx( 0 ); childIdx < m_children.size(); ++childIdx )
    {
        m_children[ childIdx ]->updateInObjectRecur( fibersCount, pCurOctree, reverseIdx );
    }
}

void SelectionTree::SelectionTreeNode::updateInBranchRecur( const int fibersCount )
{
    // TODO find a way to check if dirty and propagate
    
    vector< bool > childInBranch( fibersCount, false );
    bool atLeastOneActiveChildren( false );

    // Call update for all children.
    for( unsigned int childIdx( 0 ); childIdx < m_children.size(); ++childIdx )
    {
        m_children[ childIdx ]->updateInBranchRecur( fibersCount );
        
        if( m_children[ childIdx ]->m_pSelObject->getIsActive() )
        {
            // Get the inBranch and the state, and combine.
            // First false probably not needed.
            combineBoolVectors( childInBranch, m_children[ childIdx ]->m_pSelObject->m_inBranch, false, false );
            atLeastOneActiveChildren = true;
        }
    }
    
    m_pSelObject->m_inBranch.assign( fibersCount, false );
    
    if( m_children.empty() || !atLeastOneActiveChildren )
    {
        // This node does not have any active child.
        // Use the elements in the current box to update the m_inBranch.
        bool negate( m_pSelObject->getIsNOT() );
        
        for( unsigned int elemIdx( 0 ); elemIdx < m_pSelObject->m_inBox.size(); ++elemIdx )
        {
            m_pSelObject->m_inBranch[ elemIdx ] = ( !negate ? m_pSelObject->m_inBox[ elemIdx ] :
                                                   !m_pSelObject->m_inBox[ elemIdx ] );
        }
        
        return;
    }

    
    // Combine the child state with the current.
    bool negate( m_pSelObject->getIsNOT() );
    for( unsigned int elemIdx( 0 ); elemIdx < m_pSelObject->m_inBox.size(); ++elemIdx )
    {
        m_pSelObject->m_inBranch[ elemIdx ] = ( !negate ? m_pSelObject->m_inBox[ elemIdx ] :
                                               !m_pSelObject->m_inBox[ elemIdx ] ) & childInBranch[ elemIdx ];
    }
    
}

int SelectionTree::SelectionTreeNode::getId() const
{
    return m_nodeId;
}

SelectionTree::SelectionTreeNode::~SelectionTreeNode()
{
    // TODO delete selection object. Ownership should be transferred to the selection tree.
    for( vector< SelectionTreeNode* >::iterator nodeIt( m_children.begin() );
        nodeIt != m_children.end(); 
        ++nodeIt )
    {
        delete *nodeIt;
        *nodeIt = NULL;
    }
}

/////
// SelectionTree methods
/////

SelectionTree::SelectionTree()
    : m_pRootNode( NULL ),
      m_nextNodeId( 0 )
{}

int SelectionTree::setRoot( SelectionObject *pRootSelObject )
{
    if( m_pRootNode == NULL )
    {
        // Create new node, add it
        m_pRootNode = new SelectionTreeNode( m_nextNodeId, pRootSelObject );
        ++m_nextNodeId;
    }
    else
    {
        // Set the selection object of the node.
        m_pRootNode->setSelectionObject( pRootSelObject );
    }
    
    // Update the selection because it has changed.
    // TODO
    
    return m_pRootNode->getId();
}

int SelectionTree::addChildrenObject( const int parentId, SelectionObject *pSelObject )
{
    // Find the node with id, if possible.
    SelectionTreeNode * const pParentNode = m_pRootNode->findNode( parentId );
    
    // If found, create the children node with the next id and set its selection object.
    if( pParentNode != NULL )
    {
        // Add the node to the children of the current.
        SelectionTreeNode *pChildrenNode = new SelectionTreeNode( m_nextNodeId, pSelObject );
        pParentNode->addChildren( pChildrenNode );
        
        // TODO update selection
        
        // Increment the nextId.
        ++m_nextNodeId;
        return pChildrenNode->getId();
    }

    return -1;
}

bool SelectionTree::removeObject( const int nodeId )
{
    // Find the node with the id
    if( m_pRootNode == NULL )
    {
        return false;
    }
    
    if( m_pRootNode->getId() == nodeId )
    {
        // We want to remove the current node. Propagate to children.
        m_pRootNode->removeAllChildren();

        delete m_pRootNode;
        m_pRootNode = NULL;
        
        return true;
    }
    
    SelectionTreeNode * const pParentNode = m_pRootNode->findParentNode( nodeId );
    // if found, remove it
    
    if( pParentNode != NULL )
    {
        pParentNode->removeChildren( nodeId );
        
        // TODO Update selection
        
        return true;
    }
    
    return false;
}

SelectionObject* SelectionTree::getObject( const int itemId ) const
{
    if( m_pRootNode != NULL )
    {
        SelectionTreeNode *pNode = m_pRootNode->findNode( itemId );
        
        if( pNode != NULL )
        {
            return pNode->getSelectionObject();
        }
    }
    
    return NULL;
}

SelectionObject* SelectionTree::getParentObject( SelectionObject *pSelObj ) const
{
    if( m_pRootNode == NULL )
    {
        return NULL;
    }
    
    SelectionTreeNode * const pTreeNode = m_pRootNode->findNode( pSelObj );
    
    if( pTreeNode != NULL )
    {
        SelectionTreeNode * const pParentNode = m_pRootNode->findParentNode( pTreeNode->getId() );
        
        if( pParentNode != NULL )
        {
            return pParentNode->getSelectionObject();
        }
    }
    
    return NULL;
}

SelectionObjectVector SelectionTree::getAllObjects() const
{
    SelectionObjectVector selObj;
    
    if( m_pRootNode != NULL )
    {
        selObj = m_pRootNode->getAllSelectionObjects();
    }

    return selObj;
}

SelectionObjectVector SelectionTree::getChildrenObjects( const int itemId ) const
{
    SelectionObjectVector selObjs;
    
    if( m_pRootNode != NULL )
    {
        SelectionTreeNode *pNode = m_pRootNode->findNode( itemId );
        
        if( pNode != NULL )
        {
            selObjs = pNode->getAllChildrenSelectionObjects();
        }
    }
    
    return selObjs;
}

int SelectionTree::getActiveChildrenObjectsCount( SelectionObject *pSelObj ) const
{
    int activeChildrenCount( 0 );
    
    if( m_pRootNode != NULL )
    {
        SelectionTreeNode * const pTreeNode = m_pRootNode->findNode( pSelObj );
     
        if( pTreeNode != NULL )
        {
            activeChildrenCount = pTreeNode->getActiveDirectChildrenCount();
        }
    }
    
    return activeChildrenCount;
}

bool SelectionTree::containsId( const int itemId ) const
{
    if( m_pRootNode == NULL )
    {
        return false;
    }
    
    SelectionTreeNode *pFoundNode = m_pRootNode->findNode( itemId );
    
    return pFoundNode != NULL;
}

vector< bool > SelectionTree::getSelectedFibers( const Fibers* const pFibers )
{
    if( pFibers == NULL )
    {
        // TODO determine what we do
    }
    
    if( m_pRootNode == NULL )
    {
        // TODO determine what we do
    }
    
    const int fibersCount( pFibers->getFibersCount() );
    const vector< int > reverseIndex( pFibers->getReverseIdx() );
    
    Octree *pCurOctree( pFibers->getOctree() );
    
    // Update all selection objects to make sure that each of them knows which 
    // fibers is in it.
    m_pRootNode->updateInObjectRecur( fibersCount, pCurOctree, reverseIndex);
    
    // Update all selection objects to make sure they take into account their state
    // and the selected fibers of its children.
    m_pRootNode->updateInBranchRecur( fibersCount );
    
    // TODO remove just to silence xcode temporarily
    return m_pRootNode->getSelectionObject()->m_inBranch;
}

SelectionTree::~SelectionTree()
{
    delete m_pRootNode;
    m_pRootNode = NULL;
}
