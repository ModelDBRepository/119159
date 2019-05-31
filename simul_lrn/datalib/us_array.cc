/*********************************************************************
 file us_array.cc                                                       
 contains the implementation of class members of class               
 us_array                                                             
*********************************************************************/

#include "us_array.h"

/*********************************************************************
 implementation of class us_array                                       
*********************************************************************/


/*********************************************************************
 constructor
*********************************************************************/

template <class type>
us_array<type>::us_array(int d)
{
  dim= d;
  thelist= new label_dlist<void *>;
  this->sz= 0;
}

template <class type>
us_array<type>::~us_array()
{
  delete thelist;
}



/*********************************************************************
 public member function add                                       
 takes an array of labels and data and stores them at the appropriate
 place
*********************************************************************/

#define A 0.75
#define B 0.25
#define weight_avg(X,Y) (A*(*(X))+B*(*(Y)))

template <class type>
void us_array<type>::add(int *lb, const type &dat)
{
  label_dlist<void *> *lst= thelist; 
  label_dlist<void *> *lst2; 
  label_dlist_item<type> *it;
  
  for (int d= 0; d < dim-1; d++) {
    lst2= (label_dlist<void *> *) (lst->label_find_data(lb[d]));
    if (lst2 == NULL) {
      lst2= new label_dlist<void *>;
      lst->add(lb[d],(void *)lst2);
    }
    lst= lst2;
  }
  it= (label_dlist_item<type> *) (lst->label_find(lb[dim-1]));
  if (it == NULL) {
    lst->add(lb[dim-1], dat);
    this->sz++;
  }
  else {
    *(it->data)= weight_avg(it->data, dat); 
  }
}    


template <class type>
int us_array<type>::in(int *lb)
{
  int d= 0;
  label_dlist<void *> *lst= thelist;
  while ((lst != NULL) && (d < dim-1)) {
    lst= (label_dlist<void *> *) (lst->label_find_data(lb[d++]));
  }
  if (d < dim-1) {
    return 0;
  }
  else {
    return lst->label_in(lb[dim-1]);
  }
}


template <class type>
int us_array<type>::del(int *lb)
{
  int d= 0;
  label_dlist<void *> *lst[dim];
  lst[0]= thelist;
  while ((lst[d] != NULL) && (d < dim-1)) {
    lst[d+1]= (label_dlist<void *> *) (lst[d]->label_find_data(lb[d]));
    d++;
  }
  if (d < dim-1) {
    return 0;
  }
  else {
    lst[dim-1]->label_del(lb[dim-1]);
    for (int d= dim-1; d > 0; d--) {
      if (lst[d]->empty()) {
	lst[d-1]->del(lst[d]);
	delete lst[d];
      }
    }
    this->sz--;
    return 1;
  }
}

template <class type>
void us_array<type>::neigh_add(int *lb, int *lbreal, label_dlist<type> *lst, slset<tupel<int *, type> > *nset)
{
  assert(lst != NULL);
  int *id;
  tupel<int *, type> tp;
  label_dlist_item<type> *entry;
  if (!lst->find_left(lb[dim-1], entry) && (entry->label <= lb[dim-1])
      && (entry->next != NULL)) {
    // the label not found, it IS <= lb[dim-1] and there IS a right neighbour
    lbreal[dim-1]= entry->next->label;
    id= new int[dim]; for (int i= 0; i < dim; i++) id[i]= lbreal[i];
    tp.x= id;
    tp.y= entry->next->data;
    nset->add(tp);
  }
  lbreal[dim-1]= entry->label;
  id= new int[dim]; for (int i= 0; i < dim; i++) id[i]= lbreal[i];
  tp.x= id;
  tp.y= entry->data;
  nset->add(tp);
}


template <class type>
void us_array<type>::neigh(int *lb, int *lbreal, label_dlist<void *> *lst, int d, slset<tupel<int *,type> > *nset)
{
  assert(lst != NULL);
  label_dlist_item<void *> *item;
  if (d == dim-2) {
    if (!lst->find_left(lb[d], item) && (item->label <= lb[d])
	&& (item->next != NULL)) {
      // the label not found, lb IS greater/equal and there IS a right neighb.
      lbreal[d]= item->next->label;
      neigh_add(lb, lbreal, (label_dlist<type> *) item->next->data, nset);
    }
    lbreal[d]= item->label;
    neigh_add(lb, lbreal, (label_dlist<type> *) item->data, nset);
  }
  else {
    if (!lst->find_left(lb[d], item) && (item->label <= lb[d])
	&& (item->next != NULL)) {
      // the label not found, lb IS greater/equal and there IS a right neighb.
      lbreal[d]= item->next->label;
      neigh(lb, lbreal, (label_dlist<void *> *) item->next->data, d+1, nset);
    }
    lbreal[d]= item->label;
    neigh(lb, lbreal, (label_dlist<void *> *) item->data, d+1, nset);
  }
}

template <class type>
void us_array<type>::neighbors(int *lb, slset<tupel<int *,type> > *nset)
{
  int lbreal[dim];
  neigh(lb, lbreal, thelist, 0, nset);
}
