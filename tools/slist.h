typedef struct slist_s slist;

struct slist_s {
  char* ele;
  slist *nxt;
};




/* finds element "newele" in list "list" (NULL: not in list) */
slist* find_lele(slist **list, char* newele);



/* delete element newele from list */
void del_lele(slist **list, char* newele);



/*   insert element "newele" into sorted list "list".
 *   "list" is entry-point of sorted list.
 *   pp and p point to the element before and after newele
 */
void ins_lele(slist **list, char* newele);



/* delete and dellocate all elements from list */
void zap_lele(slist **list);
