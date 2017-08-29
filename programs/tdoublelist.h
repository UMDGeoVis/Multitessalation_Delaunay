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
//   file   : tdoublelist.h
//   author : Christian Melchiorre
//
//   Implementation of class for doubly linked list of elements
//   type T (template)
//


#ifndef _TDOUBLELIST_H
#define _TDOUBLELIST_H

#include "defs.h"
#include "error.h"


//template <class T> class TDoubleListNode;
template <class T> class TDoubleList;
template <class T> class TDoubleListIterator;
template <class T> class TDoubleListGC;

// ---------------------------------------------------------------------------------
//
//   class TDobleListNode<T>
//
//   node of a doubly linked list
//


template <class T> class TDoubleListNode
{
    private:
    
       // 
       // Default constructor is private: the only way to create instances
       // of TDoubleListNode<T> is from inside friend class TDoubleList<T>
       //
       
       TDoubleListNode()
          { error("private constructor TDoubleListNode<T>::TDoubleListNode() called"); };
	  
       TDoubleListNode( T iobject, TDoubleListNode<T> *iprev, TDoubleListNode<T> *inext )
          : next(inext), prev(iprev), object(iobject)
	    {};
	
       TDoubleListNode<T> *next, *prev;
       
    public:
     
        T object;
	
	friend class TDoubleList<T>;
	friend class TDoubleListIterator<T>;
	friend class TDoubleListGC<T>;

};


// ---------------------------------------------------------------------------------
//
//   class TDobleList<T>
//
//   Doubly linked list
//

template <class T> class TDoubleList
{
   private:
   
      TDoubleListNode<T> *first, *last;
      
      int lenght;
      
   public:
   
      TDoubleList()
         : first(NULL), last(NULL), lenght(0)
	   {};
	   
      boolean IsEmpty() { return( first == NULL ); };
      
      int Lenght() { return(lenght); };


//
//   T TDoubleList<T>::GetHead()
//
//   Return the object of type T contained in the head of the list.
//   Error if list is empty.
//

      T GetHead()
      {
         check( (first == NULL), "TDoubleList<T>::GetHead(), called on an empty list");
         return( first->object );
      }


//
//   T TDoubleList<T>::GetLast()
//
//   Return the object of type T contained in the last node of the list.
//   Error if list is empty.
//

      T GetLast()
      {
        check( (first == NULL), "TDoubleList<T>::GetLast(), called on an empty list");
        return( last->object );
      }


//
//   void TDoubleList<T>::AddHead/AddTail( T object )
//
//   Add object to the head/tail of the list.
//
    
     void AddHead( T object )
     {

   #ifdef _GC_ON
       TDoubleListNode<T> *NewNode = TDoubleListGC<T>::NewNode( object, NULL, this->first );
   #else
       TDoubleListNode<T> *NewNode = new TDoubleListNode<T>( object, NULL, this->first );
   #endif

       lenght++;
  
       if (first != NULL)  first->prev = NewNode;
       else last = NewNode;   // ...la lista era vuota

       first = NewNode;   
     }


     void AddTail(  T object )
     {

   #ifdef _GC_ON
       TDoubleListNode<T> *NewNode = TDoubleListGC<T>::NewNode( object, this->last, NULL );
   #else
       TDoubleListNode<T> *NewNode = new TDoubleListNode<T>( object, this->last, NULL );
   #endif

       lenght++;

       if (last != NULL)  last->next = NewNode;
       else  first = NewNode;   // ...la lista era vuota
     
       last = NewNode;
     }


//
//   T TDoubleList<T>::RemoveHead()
//
//   Remove the first element of the list, return the object T contained
//   in it.
//
      
      T RemoveHead()
{
   check( (first==NULL), "TDoubleList<T>::RemoveHead(), called on a empty list");

   TDoubleListNode<T> *OldNode = this->first;


   if (first->next != NULL) first->next->prev = NULL;
   
   if ( first == last ) // ...Un solo elemento in lista, ora resta vuota
     first = last = NULL;
   else
     first = first->next;
     
   T robj = OldNode->object;

#ifdef _GC_ON
   TDoubleListGC<T>::DeleteNode( OldNode );
#else
   delete( OldNode );
#endif

   lenght--;

   return(robj);

}

      T RemoveLast()
      {
   check( (first==NULL), "TDoubleList<T>::RemoveLast(), called on a empty list");

   TDoubleListNode<T> *OldNode = this->last;

   if ( last->prev != NULL ) last->prev->next = NULL;

   if ( first == last )
      first = last = NULL;
   else
      last = last->prev;

   T robj = OldNode->object;

#ifdef _GC_ON
   TDoubleListGC<T>::DeleteNode( OldNode );
#else
   delete( OldNode );
#endif

   lenght--;

   return( robj );

}


//
//   void TDoubleList<T>::AddAfter( TDoubleListNode<T>* Node, T object )
//
//   Insert the object T in the list after the node pointed by the first
//   parameter. Usually this function is called with first parameter equal
//   to Current() of an iterator on the same list:
//
//      TDoubleListIterator<T> I( &L );
//          .
//          .
//          .
//      L.AddAfter( I.Current(), t );
//
//   ATTENTION: no check is made to make sure that pointer Node points to a
//   node belonging to the list on which method AddAfter() is called.
//   This function, as well as functions AddBefore(), RemoveAt(),
//   RemoveAfter(), RemoveBefore()... must be used with great care.
//   On the other hand, the check would be too expensive.
//

void AddAfter( TDoubleListNode<T> *Node, T object )
{
   check( (Node == NULL), "TDoubleList<T>::AddAfter(), called with NULL parameter" );
   check( (first == NULL), "TDoubleList<T>::AddAfter(), called on an empty list" );

#ifdef _GC_ON
   TDoubleListNode<T> *NewNode = TDoubleListGC<T>::NewNode( object, Node, Node->next );
#else
   TDoubleListNode<T> *NewNode = new TDoubleListNode<T>( object, Node, Node->next );
#endif

   lenght++;

   if ( Node->next != NULL )
      Node->next->prev = NewNode;
   else // ...Node was the last node of the list
      last = NewNode;

   Node->next = NewNode;
}



//
//   void TDoubleList<T>::AddBefore( TDoubleListNode<T>* Node, T object )
//
//   Insert object T in the list before the node pointed by the first
//   parameter.
//



void AddBefore( TDoubleListNode<T> *Node, T object )
{
   check( (Node == NULL), "TDoubleList<T>::AddBefore(), called with NULL parameter" );
   check( (first == NULL), "TDoubleList<T>::AddBefore(), called on an empty list" );

#ifdef _GC_ON
   TDoubleListNode<T> *NewNode = TDoubleListGC<T>::NewNode( object, Node->prev, Node );
#else
   TDoubleListNode<T> *NewNode = new TDoubleListNode<T>(object, Node->prev, Node );
#endif

   lenght++;

   if ( Node->prev != NULL )
      Node->prev->next = NewNode;
   else // ...Node was the first element of the list
      first = NewNode;

   Node->prev = NewNode;

}



//
//   T TDoubleList<T>::RemoveAt( TDoubleListNode<T> * )
//
//   Remove the node pointed by Node and return the object contained in it.
//


T RemoveAt( TDoubleListNode<T> *Node )
{

   check( (Node == NULL), "TDoubleList<T>::RemoveAt(), called with NULL parameter" );
   check( (first == NULL), "TDoubleList<T>::RemoveAt(), called on an empty list" );

   if ( first == Node ) first = first->next;
   if ( last == Node ) last = last->prev;

   if ( Node->prev != NULL ) Node->prev->next = Node->next;
   if ( Node->next != NULL ) Node->next->prev = Node->prev;

   T robj = Node->object;

#ifdef _GC_ON
   TDoubleListGC<T>::DeleteNode( Node );
#else
   delete( Node );
#endif

   lenght--;

   return(robj);

}


// ---------------------------------------------------------------------------------
//
//   T TDoubleList<T>::RemoveAfter( TDoubleListNode<T> * )
//
//   Remove the node following Node in the list and return the object
//   contained in it.
//


T RemoveAfter( TDoubleListNode<T> *Node )
{

   check( (Node == NULL), "TDoubleList<T>::RemoveAfter(), called with NULL parameter" );
   check( (Node->next == NULL),
      "TDoubleList<T>::RemoveAfter(), called on the last element of the list" );
   check( (first == NULL), "TDoubleList<T>::RemoveAfter(), called on an empty list" );

   return( this->RemoveAt( Node->next ) );
}



// ---------------------------------------------------------------------------------
//
//   T TDoubleList<T>::RemoveBefore( TDoubleListNode<T> * )
//
//   Remove the node preceding Node in the list and return the object
//   contained in it.
//


T RemoveBefore( TDoubleListNode<T> *Node )
{

   check( (Node == NULL), "TDoubleList<T>::RemoveBefore(), called with NULL parameter" );
   check( (Node->prev == NULL),
      "TDoubleList<T>::RemoveBefore(), called on the last element of the list" );
   check( (first == NULL), "TDoubleList<T>::RemoveBefore(), called on an empty list" );

   return( this->RemoveAt( Node->prev ) );
}

//
//  void Append(TDoubleList<T>  L)
//
//  Append the elements of L to this list. This function is destructive:
//  it empties list L and moves its elements to this list.
//

void Append(TDoubleList<T> & L)
{
  while(! L.IsEmpty())
  {
     AddTail(L.RemoveHead());
  }
}


//
//   void TDoubleList<T>::ClearList()
//
//   Empty this list from its nodes.
//

void ClearList()
{
   TDoubleListNode<T> *Node, *NextNode;

   Node = first;

   while ( Node != NULL )
   {
      NextNode = Node->next;

#ifdef _GC_ON
      TDoubleListGC<T>::DeleteNode( Node );
#else
      delete( Node );
#endif

      Node = NextNode;
   }

   first = last = NULL;
   lenght = 0;
}

      friend class TDoubleListIterator<T>;
};



// ---------------------------------------------------------------------------------
//
//   class TDobleListIterator<T>
//
//   Iteratore for doubly linked list.
//


template <class T> class TDoubleListIterator
{

   private:

       TDoubleList<T> *list;
       TDoubleListNode<T> *current;

   public:

       TDoubleListIterator( TDoubleList<T> *ilist )
	  : list(ilist), current(ilist->first){};

//
//   void TDoubleListIterator<T>::Restart()
//
//   Move the iterator back to the first node of the list.
//

void Restart()
{
   check( (list == NULL),
     "TDoubleListIterator<T>::Restart(), iterator pointing to NULL list" );
   current = list->first;
}


//
//   void TDoubleListIterator<T>::GoNext()
//
//   Move the iterator one node forward. If iterator is at the end of
//   the list (current = NULL) then do nothing.
//

void GoNext()
{
   check( (list == NULL),"TDoubleListIterator<T>::GoNext(), iterator pointing to NULL list" );
   check( (current == NULL),"TDoubleListIterator<T>::GoNext(), current undefined" );
   current = current->next;
}

//
//   void TDoubleListIterator<T>::GoPrev()
//
//   Move the iterator one node backward. If iterator is at the beginning
//   of the list (current = first) then do nothing.
//

void GoPrev()
{
   check( (list == NULL),"TDoubleListIterator<T>::GoPrev(), iterator pointing to NULL list" );
   check( (current == NULL),"TDoubleListIterator<T>::GoPrev(), current undefined" );
   if ( current != list->first ) current = current->prev;
}

//
//   void TDoubleListIterator<T>::GoLast()
//
//   Move iterator to the last node of the list. If list is empty, then
//   do nothing.
//

void GoLast()
{
   check( (list == NULL),"TDoubleListIterator<T>::GoLast(), iterator pointing to NULL list" );
   current = list->last;
}


       //
       // Note that EndOfList is true when current == NULL, not when it
       // points to the last node of the list.
       //

       boolean EndOfList() { return( current == NULL ); };
       boolean StartOfList() { return( current == list->first ); };

       TDoubleListNode<T> *Current() { return(current); };
};




// ---------------------------------------------------------------------------------
//
//   class TDoubleListGC<T>
//
//   Static class acting as Garbage Collector for all instances of 
//   TDoubleList<T>
//

#ifdef _GC_ON

template <class T> class TDoubleListGC
{
   private:

       static TDoubleListNode<T> *GCList=NULL;
       static int GCListSize=0;
       static int GCCapacity = GC_CAPACITY;

//
//   TDoubleListNode<T> *TDoubleListGC<T>::
//             NewNode( T, TDoubleListNode<T> *, TDoubleListNode<T> * );
//
//   Allocate a new node, first searching if a node is available in GCList.
//


static TDoubleListNode<T>* NewNode( T object, TDoubleListNode<T> *prev, TDoubleListNode<T> *next )
{
   TDoubleListNode<T> *retNode = NULL;

   if ( GCList != NULL )
   {
       retNode = GCList;
       GCList = GCList->next;
       GCListSize--;

       retNode->object = object;
       retNode->prev   = prev;
       retNode->next   = next;
   }
   else
   {
       retNode = new TDoubleListNode<T>( object, prev, next );

       #ifdef DEBUG
	 DEBUG << "\nGC(size=" << GCListSize << ") : new node allocated" << endl;
       #endif // DEBUG
   }

   check( (retNode == NULL), "TDoubleListGC<T>::NewNode(), insufficient memory");

   return( retNode );

}


//
//   void TDoubleListGC<T>::DeleteNode( TDoubleListNode<T> * );
//
//   Delete a node, without freeing its memory if there is still space 
//   left in GCList.
//

static void DeleteNode( TDoubleListNode<T> *Node )
{
    if ( GCListSize < GCCapacity )
    {
       Node->prev = NULL;
       Node->next = GCList;
       GCList = Node;
       GCListSize++;
    }
    else {
       #ifdef DEBUG
	 DEBUG << "\nGC: node memory freed";
       #endif // DEBUG

       delete( Node );
    }
}

       static void SetCapacity( int newCapacity ) { GCCapacity = newCapacity; };
       friend class TDoubleList<T>;
};

#endif // _GC_ON
#endif // _TDOUBLELIST_H
