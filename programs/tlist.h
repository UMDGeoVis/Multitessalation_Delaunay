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

// ---------------------------------------------------------------------------------
//
//   file   : tlist.h
//   author : Christian Melchiorre
//
//   Implementation of a class for a doubly-linekd list of elements of 
//   type T (template)
//



#ifndef _TLIST_H
#define _TLIST_H

#include "defs.h"
#include "error.h"


template <class T> class TListNode;
template <class T> class TList;
template <class T> class TListIterator;
template <class T> class TListGC;

// ---------------------------------------------------------------------------------
//
//   class TListNode<T>
//
//   Node of a linked list
//


template <class T> class TListNode
{
    private: // PROVATO A RIMETTERLO PRIVATO, VED SE FUNGE %PAOLA%
       //
       // Default constructor is private: the only way to create instances
       // of TListNode<T> is from inside friend class TList<T>
       //

       TListNode()
	  { error("private constructor TListNode<T>::TListNode() called"); };

       TListNode( T iobject, TListNode<T> *inext )
	  : next(inext), object(iobject)
	    {};

       TListNode<T> *next;

    public:

	T object;
	
	TListNode<T>* AttachList() {return(next);}
	TListNode<T>* Succ() { return next;}
	int IsEmpty() { return object==NULL;}

	friend class TList<T>;
	friend class TListIterator<T>;
	friend class TListGC<T>;

};



// ---------------------------------------------------------------------------------
//
//   class TList<T>
//
//   Linked list
//

template <class T> class TList
{
   private:

      TListNode<T> *first;
      TListNode<T> *last;

   public:  // protected

	   TList() : first(NULL),last(NULL) {};

      void Empty() { first->next = NULL; }
      int IsEmpty() { return( first == NULL ); }
      int OnlyOne() { return( (first->next) == NULL ); }
	  
      TListNode<T>*  GetFirst() { return first; }; 
      void Assign(TListNode<T> *node) { first = node; };  // formerly (*node).next


//
//   void TList<T>::AddHead/AddTail( T object )
//  
//   Add object to the head / tail of the list.
//
    
     void AddHead( T object )
     {

   #ifdef _GC_ON
       TListNode<T> *NewNode = TListGC<T>::NewNode( object, this->first );
   #else
       TListNode<T> *NewNode = new TListNode<T>( object, this->first );
   #endif

       if (first == NULL)  last = NewNode;   // list was empty
       first = NewNode;   
     }
  


     void AddTail(  T object )
     {

   #ifdef _GC_ON
       TListNode<T> *NewNode = TListGC<T>::NewNode( object, this->last );
   #else
       TListNode<T> *NewNode = new TListNode<T>( object, this->last );
   #endif

      if (last != NULL)  last->next = NewNode;
      else  first = NewNode;   // list was empty
     
      last = NewNode;
      last->next=NULL;
     }




//
//   T TList<T>::GetHead()
//
//   Return the object of type T that is contained in the head of the list.
//   If list is empty raises an error.
//

      T GetHead() 
      {
        check( (first == NULL), "TList<T>::GetHead(), called on an empty list");
        return( first->object );
      }


// ---------------------------------------------------------------------------------
//
//   T TList<T>::RemoveHead()
//
//   Remove the element at the head of the list, return the object T
//   contained in it.
//
   
      T RemoveHead()
      {
        check( (first==NULL), "TList<T>::RemoveHead(), called on a empty list");
        
        TListNode<T> *OldNode = this->first;
        first = first->next;
        T robj = OldNode->object;

        #ifdef _GC_ON
           TListGC<T>::DeleteNode( OldNode );
        #else
           delete( OldNode );
        #endif

        return(robj);
   
      }


// ---------------------------------------------------------------------------------
//
//   T TList<T>::Pop()
//
//   Remove the element at the head of the list, return the object T
//   contained in it, the only difference w.r.t. RemoveHead() is that it
//   does not raise an error but returns NULL instead.
//

     T Pop()
     {
       if (first==NULL) return NULL;

       TListNode<T> *OldNode = this->first;
       first = first->next;
       T robj = OldNode->object;

       #ifdef _GC_ON
          TListGC<T>::DeleteNode( OldNode );
       #else
          delete( OldNode );
       #endif

       return(robj);   
     }  
 
      
// ---------------------------------------------------------------------------------
//
//   void TList<T>::AddAfter( TListNode<T>* Node, T object )
//
//   Insert object in the list after the node pointed by the first parameter,
//   this function will typically be called with the first argument being 
//   the Current() of an iterator on the same list:
//
//      TListIterator<T> I( &L );
//          .
//          .
//          .
//      L.AddAfter( I.Current(), t );
//
//   REMARK: no check is made to make sure that pointer Node points to a
//   node belonging to the list on which method AddAfter() is called.
//   This method must be used with great care. On the other side, such 
//   check would have been too expensive.
//        

      void AddAfter( TListNode<T> *Node, T object )
      {
        check( (Node == NULL), "TList<T>::AddAfter(), called with NULL parameter" );
        check( (first == NULL), "TList<T>::AddAfter(), called on an empty list" );
   
        #ifdef _GC_ON
           TListNode<T> *NewNode = TListGC<T>::NewNode( object, Node->next );
        #else
           TListNode<T> *NewNode = new TListNode<T>( object, Node->next );
        #endif
   
       Node->next = NewNode;
      };
      
    
      
// ---------------------------------------------------------------------------------
//
//   T TList<T>::RemoveAfter( TListNode<T> * )
//
//   Remove the node following a given node in a list and return the 
//   object contained in it.
//

      T RemoveAfter( TListNode<T> *Node )
      {

         check( (Node == NULL), "TList<T>::RemoveAfter(), called with NULL parameter" );
         check( (Node->next == NULL),"TList<T>::RemoveAfter(), called on the last element of the list" );
         check( (first == NULL), "TList<T>::RemoveAfter(), called on an empty list" );

         TListNode<T> *OldNode = Node->next;
         Node->next = Node->next->next;
         T robj = OldNode->object;

      #ifdef _GC_ON
         TListGC<T>::DeleteNode( OldNode );
      #else
         delete( OldNode );
      #endif

         return(robj);
      }


// ---------------------------------------------------------------------------------
//
//   void TList<T>::ClearList()
//
//   Empty this list by deleting all its nodes.
//

      void ClearList()
      {
         TListNode<T> *Node, *NextNode;
         Node = first;        

         while ( Node != NULL )
         {
            NextNode = Node->next;

      #ifdef _GC_ON
            TListGC<T>::DeleteNode( Node );
      #else
            delete( Node );
      #endif

            Node = NextNode;
         }
         first = NULL;
      }

      

      friend class TListIterator<T>;

};




// ---------------------------------------------------------------------------------
//
//   class TListIterator<T>
//
//   Iteratore for the linked list
//


template <class T> class TListIterator
{

   private:
   
       TList<T> *list;
       TListNode<T> *current;
       
   public:
   
       TListIterator( TList<T> *ilist )
          : list(ilist), current(ilist->first)
	    {};
       
 
       
// ---------------------------------------------------------------------------------
//
//   void TListIterator<T>::Restart()
//
//   Reposition iterator at the first node of the list.
//
     
       void Restart()
       {
         check( (list == NULL),"TListIterator<T>::Restart(), iterator pointing to NULL list" );
         current = list->first;
       }

   
// ---------------------------------------------------------------------------------
//
//   void TListIterator<T>::GoNext()
//
//   Moves iterator one node forward. If iterator is at the end of the list
//   (current = NULL), it does nothing.
//
      
       void GoNext()
       {
          check( (list == NULL),"TListIterator<T>::GoNext(), iterator pointing to NULL list" );
          check( (current == NULL),"TListIterator<T>::GoNext(), current undefined" );
          current = current->next;
       }
   

       // Check if iterator is at the beginning / end of the list.
       // Note that EndOfList is true when current == NULL, not when 
       // it points to the last node of the list.
       boolean EndOfList() { return( current == NULL ); };
       boolean StartOfList() { return( current == list->first ); };
       
       TListNode<T> *Current() { return(current); };
};




// ---------------------------------------------------------------------------------
//
//   class TListGC<T>
// 
//   Static class used as a Garbage Collector for all instances of TList<T>
//

#ifdef _GC_ON


template <class T> class TListGC
{
   private:
   
       static TListNode<T> *GCList;
       static int GCListSize;
       static int GCCapacity;

//
//   TListNode<T> *TListGC<T>::NewNode( 
//                   T, TListNode<T> *, TListNode<T> * );
//
//   Allocate a new node, first checking if a node is available in GCList.
//

       static TListNode<T> *NewNode( T object, TListNode<T> *next )   
       {
          TListNode<T> *retNode = NULL;
   
          if ( GCList != NULL )
          {
              retNode = GCList;
              GCList = GCList->next;
              GCListSize--;  
       
              retNode->object = object;
              retNode->next   = next;
          }
          else
          {
              retNode = new TListNode<T>( object, next );
     
            #ifdef DEBUG
              DEBUG << "\nGC(size=" << GCListSize << ") : allocato nuovo nodo";
            #endif // DEBUG
          }
   
         check( (retNode == NULL), "TListGC<T>::NewNode(), insufficient memory");
   
         return( retNode );

       }    


// ---------------------------------------------------------------------------------
//
//   void TListGC<T>::DeleteNode( TListNode<T> * );
//
//   Delete a new node, without releasing memory if there is still
//   space in GCList.
//
  
       static void DeleteNode( TListNode<T> *Node )
       {

           if ( GCListSize < GCCapacity )
           {
              Node->next = GCList;
              GCList = Node;
              GCListSize++;
           } 
           else 
	   {
            #ifdef DEBUG
              DEBUG << "\nGC: node has been deallocated";
            #endif // DEBUG
 
              delete( Node );
           }
       }  

  
       static void SetCapacity( int newCapacity ) { GCCapacity = newCapacity; };

       friend class TList<T>;
  	 
};


#endif // _GC_ON

#endif // _TLIST_H
