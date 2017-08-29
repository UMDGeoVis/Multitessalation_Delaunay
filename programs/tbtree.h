/*****************************************************************************
Delaunay Triangulator and MT constructor, version 1.0, 1999.
Copyright (C) 1999 DISI - University of Genova, Italy.
Group of Geometric Modeling and Computer Graphics DISI.
DISI - University of Genova, Via Dodecaneso 35, 16146 Genova - ITALY.

This program is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License
as published by the Free Software Foundation; either version 2
of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
*****************************************************************************/

// ----------------------------------------------------------------------------
//
//  file   : tbtree.h
//  author : Christian Melchiorre
//
//  Definition of Balanced Binary Trees of elements of type T.
//  For a complete description of such data type, see the text book:
//  N. Wirth, "Algorithms + Data Structures = Programs" (section 4.4.6
//  and following).
//


#ifndef _BTREE_H
#define _BTREE_H

#include "defs.h"
#include "error.h"
#include "ttriang.h"

int compare( double a , double b );

template <class T> class TBTreeNode;
template <class T> class TBTree;
template <class T> class TBTreeIterator;
// template <class T> class TBTreeGC;


// 
// Constants for balancement of the nodes of the binary tree
//

const char BALANCED         = 0; //  0 in Wirth
const char UNBALANCED_LEFT  = 1; // -1 in Wirth
const char UNBALANCED_RIGHT = 2; // +1 in Wirth


// ----------------------------------------------------------------------------
//
//  class TBTreeNode
//
//  Node of a balanced binary tree of elements of type T
//


template <class T> class TBTreeNode
{
   private:
    
       // 
       // Default constructor is private: the only way to create instances
       // of TBTreeNode is from within the friend class TBTree.
       //
       
       TBTreeNode()
          { error("private constructor TBTreeNode::TBTreeNode() called"); };
	  
       TBTreeNode( T iobject, TBTreeNode<T> *ileft, TBTreeNode<T> *iright, char ibal)
          : left(ileft), right(iright), bal(ibal), object(iobject)
	    {};
	
       TBTreeNode<T>  *left, *right;
       char bal; // UNBALANCED_LEFT, BALANCED, UNBALANCED_RIGHT

   public:

       T object;

       friend class TBTree<T>;
       friend class TBTreeIterator<T>;
       	     
};



// ----------------------------------------------------------------------------
//
//  class TBTree
//
//  Balanced binary tree of elements of type T.
//  Remark: for each type T, in order to use class TBTree<T> it is 
//  is necessary to provide a function:
//
//    int compare ( T t1, T t2 )
//
//  which returns -1 if t1<t2, 0 if t1=t2, and +1 if t1>t2
//


template <class T> class TBTree
{
   private:

      // Tree root.
      TBTreeNode<T> *root;

      // Auxiliary functions for insertion and deletion
      
      boolean InsertRec( TBTreeNode<T>* &, T );
      boolean BalanceAfterInsertL( TBTreeNode<T>* & );
      boolean BalanceAfterInsertR( TBTreeNode<T>* & );
      
      boolean RemoveRec( TBTreeNode<T>* &, T, T& );
      boolean DeleteRec( TBTreeNode<T>* &, TBTreeNode<T>* & );
      boolean BalanceAfterRemoveL( TBTreeNode<T>* & );
      boolean BalanceAfterRemoveR( TBTreeNode<T>* & );
      
      boolean RemoveMinRec( TBTreeNode<T>* &, T & );
      boolean RemoveMaxRec( TBTreeNode<T>* &, T & );
      
            
   public:

      // Constructor
      TBTree() : root(NULL) {};

      // Test if tree is empty
      boolean IsEmpty() { return( root == NULL ); };
      
      // Test if an element is in the tree
      boolean IsIn( T );

      
      // Insertion and deletion of an element.
      // Function Remove() returns the removed element since the removed
      // element might be just an element having a key equal to the key
      // of the parameter.
      void Insert( T );
      T Remove( T );
      
      // Deletion of the minimum / maximum element.
      T RemoveMin();
      T RemoveMax();   

      // Return the minimum / maximum element without deleting it.
      T GetMin();
      T GetMax();  
     
      void ClearTree();

      friend class TBTreeIterator<T>;
      
      friend int compare( T, T );

      #ifdef DEBUG
         void PrintRec( TBTreeNode<T> *Node );
         void Print();
      #endif

};


// ----------------------------------------------------------------------------
//
//  class TBTreeIterator
//
//  Iteratore for balanced binary trees of elements of type T
//


template <class T> class TBTreeIterator
{
   private:
   
      TBTree<T> *tree;
       
      TBTreeNode<T> *current;   // current node
      
   public:
   
      TBTreeIterator ( TBTree<T> *itree )
         : tree(itree), current(itree->root)
	  {};
	  
      boolean CanGoLeft()  { return( current != NULL && current->left != NULL ); };
      boolean CanGoRight() { return( current != NULL && current->right != NULL ); };
      boolean CanGoParent() { return( current != NULL && current != tree->root ); };
      
      void Restart();
      void GoLeft();
      void GoRight();
      void GoParent();

      TBTreeNode<T> *Current() { return( current ); };
      
};



#ifdef _GC_ON

// ---------------------------------------------------------------------------------
//
//   class TBTreeGC<T>
// 
//   Static class used as Garbage Collector for all instances of TBtree<T>
//

//   ...TO BE IMPLEMENTED ...IF NEEDED (SEE TLIST.H/TDDOUBLELIST.H)

#endif // _GC_ON



// ----------------------------------------------------------------------------
//
//  file   : tbtree.cpp
//
//  Implementation of Balanced Binary Trees of elements of type T.
//  For a complete description of such data type, see the text book:
//  N. Wirth, "Algorithms + Data Structures = Programs" (section 4.4.6
//  and following).
//



// ----------------------------------------------------------------------------
//
//  Methods of class TBTree<T>
//




// ----------------------------------------------------------------------------
//
//  boolean TBTree<T>::IsIn( T object )
// 
//  Return TRUE iff object belongs to the set of the objects contained
//  in the nodes of this tree.
//


template <class T>  inline
boolean TBTree<T>::IsIn( T object )
{
   TBTreeNode<T> *Node = this->root;
   
   while( Node != NULL )
   {
      switch( compare( object, Node->object ) )
      {
         case  0: return( TRUE );
	 case -1: Node = Node->left; break;
	 case  1: Node = Node->right; break;
	 default: error( "TBTree<T>::Remove(), invalid compare<T> return value" );
      }
   }
   
   return( FALSE );

}



// ----------------------------------------------------------------------------
//
//  void TBTree<T>::Insert( T )
//
//  Insert an element in the balanced binary tree. It uses the private 
//  auxiliary function TBTree<T>::InsertRec().
//

template <class T>  inline
void TBTree<T>::Insert( T object )
{
   // cerr << "Ins: " << object;
   InsertRec( this->root, object );
   // cerr << " OK " << endl;
}



// ----------------------------------------------------------------------------
//
//  boolean TBTree<T>::InsertRec( TBTreeNode<T>* &, T )
//
//  Private auxiliary function for inserting an element T in the tree.
//  Called by function TBTree<T>::Insert(). The returned value corresponds 
//  to parameter h of Wirth's book, and is TRUE iff the heigth of the 
//  subtree, in which the element has been inserted, has become larger.
//

template <class T>  inline
boolean TBTree<T>::InsertRec( TBTreeNode<T>* &Node, T object )
{

   if ( Node == NULL )
   {
   
      Node = new TBTreeNode<T>( object, NULL, NULL, 0 );
      // Remark: Node is passed by reference!
      
      check( (Node == NULL), "TBTree<T>::Insert(), insufficient memory" );
      
      return( TRUE );
      
   }
   else // Node != NULL
   {
   
      switch( compare( object, Node->object ) )
      {
         
           case -1: // object < Node->object     
           {     
               if ( this->InsertRec( Node->left, object ) )
	           return( BalanceAfterInsertL( Node ) );
	       else
	           return( FALSE );
	   }
	   
	   case 1: // object > Node->object
	   {
	       if ( this->InsertRec( Node->right,  object ) )
	           return( BalanceAfterInsertR( Node ) );
	       else 
	           return( FALSE );
	       
	   }
	   
	   case 0: // object == Node->object
	   {
	      return( FALSE );
	   }
	   	   
	   default:
	   {
	      error( "TBTree<T>::Insert(), wrong compare<T>() return value" );
	      return( FALSE ); // ...to avoid warning from the compiler
	   }
       
       } // end ...switch( compare )
       
   } // end ...else (Node != NULL)

}


// ----------------------------------------------------------------------------
//
//  boolean TBTree<T>::BalanceAfterInsertL( TBTreeNode<T>* & Node )
//
//  Private auxiliary function for rebalancing the tree after the
//  insertion of an element.
//

template <class T>  inline
boolean TBTree<T>::BalanceAfterInsertL( TBTreeNode<T>* &Node )
{	       

    switch( Node->bal )
    {
        case UNBALANCED_RIGHT:	
        {	 
            Node->bal = BALANCED;
            return(FALSE);
        }
	 
        case BALANCED:
	{
	    Node->bal = UNBALANCED_LEFT;
            return( TRUE );
        }
	 
	case UNBALANCED_LEFT:
	{
	    TBTreeNode<T> *lNode = Node->left;
			    
	    if ( lNode->bal == UNBALANCED_LEFT )
	    {
	         // single LL rotation
			       
	         Node->left = lNode->right;
	         lNode->right = Node;
	         Node->bal = BALANCED;
	         Node = lNode;
	    }
	    else
	    {
	         // double LR rotation
			       
	         check( (lNode == NULL), "TBTree<T>::BalanceAfterInsertL(), <1> this should not happen" );
			       
	         TBTreeNode<T> *lrNode = lNode->right;
			       
	         lNode->right = lrNode->left;
	         lrNode->left = lNode;
	         Node->left = lrNode->right;
	         lrNode->right = Node;
			       
	         if ( lrNode->bal == UNBALANCED_LEFT )
	             Node->bal = UNBALANCED_RIGHT;
	         else
	             Node->bal = BALANCED;
				  
	         if ( lrNode->bal == UNBALANCED_RIGHT )
	             lNode->bal = UNBALANCED_LEFT;
	         else
	             lNode->bal = BALANCED;	       
			  
	         Node = lrNode;
			    			       
	    }
			    			    
	    Node->bal = BALANCED;
			    
	    return( FALSE );
        }
	 
	default:
	{
	   error( "TBTree<T>::BalanceAfterInsertL(), <2> this should not happen" ); 
	   return( FALSE ); // ...per evitare warning
	}
      
   } // end ...switch( Node->bal )

}


// ----------------------------------------------------------------------------
//
//  boolean TBTree<T>::BalanceAfterInsertR( TBTreeNode<T>* & Node )
//
//  Private auxiliary function for rebalancing the tree after the 
//  insertion of an element.
//


template <class T>  inline
boolean TBTree<T>::BalanceAfterInsertR( TBTreeNode<T>* &Node )
{

   switch( Node->bal )
   {
		   
       case UNBALANCED_LEFT:
       {
           Node->bal = BALANCED;
   	   return( FALSE );
       }
		       
       case BALANCED:
       {
           Node->bal = UNBALANCED_RIGHT;
           return( TRUE );
       }
		       
       case UNBALANCED_RIGHT:
       {
           TBTreeNode<T> *rNode = Node->right;
			  
   	   if ( rNode->bal == UNBALANCED_RIGHT )
	   {
	       // single RR rotation
	     
	       Node->right = rNode->left;
	       rNode->left = Node;
	       Node->bal = BALANCED;
	       Node = rNode;
	   }
	   else
	   {
	       // double RL rotation
			     
	       check( (rNode == NULL), "TBTree<T>::BalanceAfterInsertR(), <1> this should not happen" );
			     
	       TBTreeNode<T> *rlNode = rNode->left;
			     
	       rNode->left = rlNode->right;
	       rlNode->right = rNode;
	       Node->right = rlNode->left;
	       rlNode->left = Node;
	     
	       if ( rlNode->bal == UNBALANCED_RIGHT )
	           Node->bal = UNBALANCED_LEFT;
	       else
	           Node->bal = BALANCED;
				
	       if ( rlNode->bal == UNBALANCED_LEFT )
	           rNode->bal = UNBALANCED_RIGHT;
	       else
	           rNode->bal = BALANCED;		     
			     
	       Node = rlNode; 
			     
	   }
			  
	   Node->bal = BALANCED;
			  
	   return( FALSE );
			  
       }
		       
       default:
       {
           error( "TBTree<T>::BalanceAfterInsertR(), <2> this should not happen" );
	   return( FALSE ); // to avoid warning
       }
		       
   } // end ...switch( Node->bal )	          
	       
}



// ----------------------------------------------------------------------------
//
//  T TBTree<T>::Remove( T )
//
//  Remove the node containing T from the tree. It uses the private 
//  auxiliary function TBTree<T>::RemoveRec(). Return the removed object.
//  In fact, such object might be different from the one passed as a 
//  parameter, but simply have an equal key according to function
//  compare( T, T )
//


template <class T>  inline
T TBTree<T>::Remove( T object )
{
   // cerr << "Rem : " << object;

   T robj;

   this->RemoveRec( this->root, object, robj );

   return( robj );

   // cerr << " OK " << endl;
}



// ----------------------------------------------------------------------------
//
//  boolean TBTree<T>::RemoveRec( TBTreeNode<T>* &, T, T& )
//
//  Private auxiliary function for deleting an element from the tree.
//  Called by function TBTree<T>::Remove(). The returned vale corresponds
//  to parameter h in Wirth's book, and is TRUE iff the heigth of the
//  subtree, from which the element has been deleted, has become smaller.
//  In the parameter passed by reference (robj), the removed object is
//  returned. In fact, it might be a different object from the given one,
//  and simply have an equal key according to function compare.
//

template <class T>  inline
boolean TBTree<T>::RemoveRec( TBTreeNode<T>* &Node, T object, T& robj )
{
    
    if ( Node == NULL )
    {
        //
	// object non e' nell'albero
	//
	
	return( FALSE );
    }
    else // Node != NULL
    {
        
	switch( compare( object, Node->object ) )
	{
	    case -1: // object < Node->object
	    {
	       if ( this->RemoveRec( Node->left, object, robj ) )
	           return( BalanceAfterRemoveL( Node ) );
	       else 
	           return( FALSE );
	    }
	    
	    case  1: // object > Node->object
	    {
	       if ( this->RemoveRec( Node->right, object, robj ) )
	           return( BalanceAfterRemoveR( Node ) );
	       else
	           return( FALSE );    
	    }
	    
	    case 0: // object == Node->object
	    {
               robj = Node->object;

	       TBTreeNode<T> *dNode = Node;
	       
	       if ( dNode->right == NULL )
	       {
	          Node = dNode->left;
			  delete( dNode ); dNode = NULL;
		  return(TRUE);
	       }
	       else if ( dNode->left == NULL )
	       {
	          Node = dNode->right;
		  delete( dNode ); dNode = NULL;
		  return( TRUE );
	       }
	       else // dNode->left != NULL && dNode->right != NULL
	       {
	          if ( DeleteRec( Node, Node->left ) )
		     return( BalanceAfterRemoveL( Node ) );
		  else
		     return( FALSE );		     
	       }
	       
	    }
	    
	    default:
	    {
	       error( "TBTree<T>::Remove(), invalid compare<T> return value" );
	       return( FALSE ); // to avoid warning
	    }
	
	} // end ...switch( compare );
    
    } // end ...else ( Node != NULL )
        
}


// ----------------------------------------------------------------------------
//
//  boolean TBTree<T>::BalanceAfterRemoveL( TBTreeNode<T>* & )
//
//  Private auxiliary function for rebalancing the tree after the
//  deletion of an element (see procedure balance1 in Wirth'd s book).
//


template <class T>  inline
boolean TBTree<T>::BalanceAfterRemoveL( TBTreeNode<T>* & Node )
{
   
   switch( Node->bal )
   {
      case UNBALANCED_LEFT:
      {
         Node->bal = BALANCED;
	 return( TRUE );
      }
      
      case BALANCED:
      {
         Node->bal = UNBALANCED_RIGHT;
	 return( FALSE );
      }
      
      case UNBALANCED_RIGHT:
      {
         TBTreeNode<T> *rNode = Node->right;
	 char rBal = rNode->bal;
	 
	 if ( rBal == BALANCED || rBal == UNBALANCED_RIGHT )
	 {
	    // single RR rotation
	    
	    Node->right = rNode->left;
	    rNode->left = Node;
	    
	    if ( rBal == BALANCED )
	    {
	       Node->bal = UNBALANCED_RIGHT;
	       rNode->bal = UNBALANCED_LEFT;
	       Node = rNode;
	    
	       return( FALSE );
	    }
	    else // rBal == UNBALANCED_RIGHT
	    {
	       Node->bal = BALANCED;
	       rNode->bal = BALANCED;
	       Node = rNode;
	    
	       return( TRUE );
	    }
	    	    
	 }
	 else // rBal == UNBALANCED_LEFT
	 {
	    // double RL rotation
	    
	    check( (rNode == NULL), "TBTree<T>::BalanceAfterRemoveL(), <1> this should not happen" );
	    
	    TBTreeNode<T> *rlNode = rNode->left;
	    int rlBal = rlNode->bal;
	    
	    rNode->left = rlNode->right;
	    rlNode->right = rNode;
	    Node->right = rlNode->left;
	    rlNode->left = Node;
	    
	    if ( rlBal == UNBALANCED_RIGHT )
	        Node->bal = UNBALANCED_LEFT;
            else
	        Node->bal = BALANCED;	    
	    
	    if ( rlBal == UNBALANCED_LEFT )
	        rNode->bal = UNBALANCED_RIGHT;
	    else
	        rNode->bal = BALANCED;
	    
	    Node = rlNode;
	    rlNode->bal = BALANCED;
	    
	    return( TRUE );
	 }
	 
      }
      
      default:
      {
         error( "TBTree<T>::BalanceAfterRemoveL(), <2> this should not happen" );
	 return( FALSE ); // to avoid warning
      } 
   
   } // end ...switch( Node->bal )
   
}


// ----------------------------------------------------------------------------
//
//  boolean TBTree<T>::BalanceAfterRemoveR( TBTreeNode<T>* & Node )
//
//  Private auxiliary function for rebalancing the tree after the
//  deletion of an element (see procedure balance2 in Wirth'd s book).
//


template <class T>  inline
boolean TBTree<T>::BalanceAfterRemoveR( TBTreeNode<T>* &Node )
{

   switch( Node->bal )
   {
   
      case UNBALANCED_RIGHT:
      {
         Node->bal = BALANCED;
	 return( TRUE );
      }
      
      case BALANCED:
      {
         Node->bal = UNBALANCED_LEFT;
	 return( FALSE );
      }
      
      case UNBALANCED_LEFT:
      {
         TBTreeNode<T> *lNode = Node->left;
	 int lBal = lNode->bal;
	 
	 if ( lBal == BALANCED || lBal == UNBALANCED_LEFT )
	 {
	    // single LL rotation
	    
	    Node->left = lNode->right;
	    lNode->right = Node;
	    
	    if ( lBal == BALANCED )
	    {
	       Node->bal = UNBALANCED_LEFT;
	       lNode->bal = UNBALANCED_RIGHT;
	       Node = lNode;
	       return( FALSE );
	    }
	    else // lBal == UNBALANCED_LEFT
	    {
	       Node->bal = BALANCED;
	       lNode->bal = BALANCED;
	       Node = lNode;
	       return( TRUE );
	    }
	    
	 }
	 else // lBal == UNBALANCED_RIGHT
	 {
	    // double LR rotation
	    
            check( (lNode == NULL), "TBTree<T>::BalanceAfterRemoveR(), <1> this should not happen" );

	    TBTreeNode<T> *lrNode = lNode->right;
	    int lrBal = lrNode->bal;
	    
	    lNode->right = lrNode->left;
	    lrNode->left = lNode;
	    Node->left = lrNode->right;
	    lrNode->right = Node;
	    
	    if ( lrBal == UNBALANCED_LEFT )
	        Node->bal = UNBALANCED_RIGHT;
            else
	        Node->bal = BALANCED;
		
	    if ( lrBal == UNBALANCED_RIGHT )
	        lNode->bal = UNBALANCED_LEFT;
	    else
	        lNode->bal = BALANCED;
		
	    Node = lrNode;
	    lrNode->bal = BALANCED;
	    
	    return( TRUE );
	    
	 }
	 
      }
      
      default:
      {
         error( "TBTree<T>::BalanceAfterRemoveR(), <2> this should not happen" );
	 return( FALSE ); // to avoid warning
      }
   
   
   } // end ...switch( Node->bal )

}


// ----------------------------------------------------------------------------
//
//  boolean TBTree<T>::DeleteRec( TBTreeNode<T>* &dNode, TBTreeNode<T>* &psNode,
//                             TBTreeNode<T>* &sNode )
//
//  Auxiliary function called by TBTree<T>::RemoveRec( T ).
//  Perform the removal of node dNode, found by TBTree<T>::RemoveRec()
//  by replacing it with the maximum of its left descendants.
//  The search for such replacing node is done through recursion of
//  TBTree<T>::DeleteRec(). The returned value has the same meaning as
//  in the previous functions.
//


template <class T>  inline
boolean TBTree<T>::DeleteRec( TBTreeNode<T>* &dNode, TBTreeNode<T>* &sNode )
{
 
    if ( sNode->right != NULL )
    {
        if ( DeleteRec( dNode, sNode->right ) )
	   return( BalanceAfterRemoveR( sNode ) );
	else   
	   return( FALSE );
    }
    else
    {
        dNode->object = sNode->object;
        TBTreeNode<T> *OldsNode = sNode;
	sNode = sNode->left;
	delete( OldsNode ); OldsNode = NULL;
	return( TRUE );
    }

}


// ----------------------------------------------------------------------------
//
//  T TBTree<T>::RemoveMin()
//
//  Find the node containing the minimum element in the tree (it is
//  the leftmost node among the descendants of the root), and delete it.
//  Return the element contained in such node.
//


template <class T>  inline
T TBTree<T>::RemoveMin()
{
    check( (this->root == NULL), "TBTree<T>::RemoveMin(), called on an empty tree");
    T robj;    
    RemoveMinRec( this->root, robj );    
    return( robj );
}


// ----------------------------------------------------------------------------
//
//  boolean TBTree<T>::RemoveMinRec( TBTreeNode<T>* &Node, T& robj )
//
//  Private auxiliary function used by TBTree<T>::RemoveMin() to perform
//  the recursive search and deletion of the mimimum node. The code is
//  very similar to TBTree<T>::RemoveRec(), but the descent is always
//  continued in the left branch (no comparison with compare<T>() is done).
//  In variable robj (passed by reference) the value T, contained in the 
//  deleted node, is returned.
//


template <class T>  inline
boolean TBTree<T>::RemoveMinRec( TBTreeNode<T>* &Node, T& robj )
{

    check( (Node == NULL), "TBTree<T>::RemoveMinRec(), <1> this should not happen" );
    
    if ( Node->left != NULL )
    {
        if ( this->RemoveMinRec( Node->left, robj ) )
	    return( BalanceAfterRemoveL( Node ) );
	else
	    return( FALSE );
    }
    else // mimimun has been found, Node->left == NULL;
    {
	TBTreeNode<T> *dNode = Node;
	robj = dNode->object;	
	Node = dNode->right;
	delete( dNode ); dNode = NULL;
	return( TRUE );	
    }

}


// ----------------------------------------------------------------------------
//
//  T TBTree<T>::RemoveMax()
//
//  Find the node containing the maximum element in the tree (it is
//  the rightmost node among the descendants of the root), and delete it.
//  Return the element contained in such node. The code is very similar to
//  TBTree<T>::RemoveRec(), but the descent is always continued in the
//  right branch (no comparison with compare<T>() is done).
//

template <class T>  inline
T TBTree<T>::RemoveMax()
{
   check( (this->root == NULL), "TBTree<T>::RemoveMax(), called on an empty tree");
   T robj;
   RemoveMaxRec( this->root, robj );
   return( robj );
}


// ----------------------------------------------------------------------------
//
//  T TBTree<T>::RemoveMaxRec()
//
//  Private auxiliary function used by TBTree<T>::RemoveMax() to
//  perform the recursive search and deletion of the maximum node.
//  Similar to RemoveMinRec().
//


template <class T>  inline
boolean TBTree<T>::RemoveMaxRec( TBTreeNode<T>* &Node, T& robj )
{

    check( (Node == NULL), "TBTree<T>::RemoveMaxRec(), <1> this should not happen" );
    
    if ( Node->right != NULL )
    {
        if ( this->RemoveMaxRec( Node->right, robj ) )
	    return( BalanceAfterRemoveR( Node ) );
	else
	    return( FALSE );
    }
    else // trovato il massimo, Node->right == NULL;
    {
	TBTreeNode<T> *dNode = Node;
	robj = dNode->object;	
	Node = dNode->left;
	delete( dNode ); dNode = NULL;	
        return( TRUE );
    }

}


// ----------------------------------------------------------------------------
//
//  T TBTree<T>::GetMin()
//
//  Find the minimum element of the tree and return the object contained 
//  in the corresponding node, without deleting the node.
//


template <class T>  inline
T TBTree<T>::GetMin()
{
   check( (this->root == NULL), "TBTree<T>::GetMin(), called on an empty tree");
   TBTreeNode<T> *Node = this->root;
   while( Node->left != NULL ) Node = Node->left;
   return( Node->object );
}


// ----------------------------------------------------------------------------
//
//  T TBTree<T>::GetMax()
//
//  Find the maximum element of the tree and return the object contained 
//  in the corresponding node, without deleting the node.
//


template <class T>  inline
T TBTree<T>::GetMax()
{
   check( (this->root == NULL), "TBTree<T>::GetMax(), called on an empty tree");
   TBTreeNode<T> *Node = this->root;
   while( Node->right != NULL ) Node = Node->right;
   return( Node->object );
}


// ----------------------------------------------------------------------------
//
//  T TBTree<T>::ClearTree()
//
//  Empty the tree by deleting all its nodes.
//


template <class T>  inline
void TBTree<T>::ClearTree()
{
   while( !IsEmpty() )
      RemoveMin();
}



// ---------------------------------------------------------------------------------
//
//  Methods of class TBTreeIterator<T>
//



// ---------------------------------------------------------------------------------
//
//  void TBTreeIterator<T>::Restart()
//
//  Place iterator at the root of the tree.
//

template <class T>  inline
void TBTreeIterator<T>::Restart()
{
   check( (tree == NULL), "TBTreeIterator<T>::Restart(), iterator pointing to NULL tree" );
   
   current = tree->root;
}


// ---------------------------------------------------------------------------------
//
//  void TBTreeIterator<T>::GoLeft()
//
//  Move iterator to left child of current node.
//

template <class T>  inline
void TBTreeIterator<T>::GoLeft()
{
   check( (tree == NULL), "TBTreeIterator<T>::GoLeft(), iterator pointing to NULL tree" );
   check( (current == NULL), "TBTreeIterator<T>::GoLeft(), current undefined" );
   
   current = current->left;
}

   
// ---------------------------------------------------------------------------------
//
//  void TBTreeIterator<T>::GoRight()   
//
//  Move iterator to right child of current node.
//

template <class T>  inline
void TBTreeIterator<T>::GoRight()
{
   check( (tree == NULL), "TBTreeIterator<T>::GoRight(), iterator pointing to NULL tree" );
   check( (current == NULL), "TBTreeIterator<T>::GoRight(), current undefined" );
   
   current = current->right;

}

  
// ---------------------------------------------------------------------------------
//
//  void TBTreeIterator<T>::GoParent()   
//
//  Move back to parent of current node. 
//  Implemented by searching for the current node, starting from the root.
//  Another way would be maintaining, for each node, a pointer to the 
//  parent (but more memory needed).
//


template <class T>  inline
void TBTreeIterator<T>::GoParent()
{
   check( (tree == NULL), "TBTreeIterator<T>::GoParent(), iterator pointing to NULL tree" );
   check( (current == NULL), "TBTreeIterator<T>::GoParent(), current undefined" );
   
   TBTreeNode<T> Node = tree->root;
   TBTreeNode<T> Parent = NULL;
   
   while( Node != NULL ) 
   {
       switch( compare( Node->object, current->object ) )
       {
	    case -1: Parent = Node; Node = Node->left; break;
	    case  1: Parent = Node; Node = Node->right; break;
	    case  0: current = Parent; return; 
	    default: error( "TBTreeIterator<T>::GoParent(), invalid compare<T>() return value" );
       }    
	    	 
   } // end...while
   
   error( "TBTreeIterator<T>::GoParent(), current not found in tree" );
}


#ifdef DEBUG

template <class T>  inline
void TBTree<T>::PrintRec( TBTreeNode<T> *Node )
{
   DEBUG << Node << "[ ( " 
        << Node->object->x << ", " << Node->object->y << ", " << Node->object->z 
        << "=> " << Node->object->Error << " ), L:" ;
   if ( Node->left != NULL ) DEBUG << Node->left; else DEBUG << "NULL";
   DEBUG << ", R:";
   if ( Node->right != NULL) DEBUG << Node->right; else DEBUG << "NULL";
   DEBUG << " ]" << endl;

   if ( Node->left != NULL ) PrintRec( Node->left );
   if ( Node->right != NULL ) PrintRec( Node->right );
}

template <class T>  inline
void TBTree<T>::Print()
{
   DEBUG << endl << "TREE: " << endl;
   if (root != NULL ) PrintRec( root );
}

#endif


#endif // _BTREE_H
