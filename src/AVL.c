#include <stdlib.h>
#include <stdio.h>

#include "AVL.h"

typedef struct {   
  int incr;
  PtrAVL p_AVL;
} AVLPaire ;



void AVLafficherNoeud(PtrAVL p_a) {

  fprintf(stderr,"  %d : %d\n",p_a->val,p_a->equilibre);
  fprintf(stderr,"  / \\\n"); 
  fprintf(stderr,"%2d  %2d\n",(p_a->lSon==NULL)?-1:p_a->lSon->val,
      (p_a->rSon==NULL)?-1:p_a->rSon->val); 
  fprintf(stderr,"-----\n");

}

void AVLafficherInfixe(PtrAVL p_a) {  
  if (p_a != NULL) {
    AVLafficherInfixe(p_a->lSon);
    AVLafficherNoeud(p_a);
    AVLafficherInfixe(p_a->rSon);
  }
}
void AVLafficherPrefixe(PtrAVL p_a) {  
  if (p_a != NULL) {
    AVLafficherNoeud(p_a);
    AVLafficherPrefixe(p_a->lSon);
    AVLafficherPrefixe(p_a->rSon);
  }
}

void AVLafficher(PtrAVL p_a) {  
  fprintf(stderr,"Affichage\n");
  AVLafficherPrefixe(p_a);

}




PtrAVL   AVLnew   () {
  
  PtrAVL p_a;

  p_a = (PtrAVL) malloc(sizeof(struct AVL));
  p_a->lSon = NULL;
  p_a->rSon = NULL;
  p_a->data = NULL;
  p_a->equilibre = 0;
  //p_a->val = 0;
  
  return p_a;

}

PtrAVL   lRot     (PtrAVL p_a) {
/*#define TRACE*/ 
/*

       (A)                    (B)
       / \                    / \
      /   \                  /   \ 
     /\   (B)              (A)    \
    /Ag\  / \       -->    / \     \
1-  ---- /   \            /   \     \
        /\    \          /\   /\    /\
       /Bg\    \        /Ag\ /Bg\  /Bd\
2-     ----     \       ---- ----  ----
                /\
               /Bd\
3-             ----
*/

  PtrAVL p_b;

  int vA, vB, vAnew, vBnew;


  p_b = p_a->rSon;
  
  vA = p_a->equilibre; vB = p_b->equilibre;
  
  p_a->rSon = p_b->lSon;
  p_b->lSon = p_a;

  /* recalculer les champ equilibre */
#define MIN(x,y) (((x)<(y)) ? (x) : (y))
#define MAX(x,y) (((x)<(y)) ? (y) : (x))

  vAnew = vA - 1 - MAX(0,vB);
  vBnew = - 1 + MIN(0,vAnew) + vB;

#undef MIN
#undef MAX  

  p_a->equilibre = vAnew;
  p_b->equilibre = vBnew;

  return p_b; 

}

PtrAVL     rRot     (PtrAVL p_a) {
  /*
    
             (A)                  (B)
             / \                  / \
            /   \                /   \ 
          (B)   /\              /    (A)
          / \  /Ad\     -->    /    /  \
1-       /   \ ----           /    /    \
        /    /\              /\   /\    /\
       /    /Bd\            /Bg\ /Bd\  /Ad\
2-    /     ----            ---- ----  ----
     /\                
    /Bg\                             
3-  ----   
*/
  PtrAVL p_b;

  int vA, vB, vAnew, vBnew;

  p_b = p_a->lSon;
  
  vA = p_a->equilibre; vB = p_b->equilibre;
  
  p_a->lSon = p_b->rSon;
  p_b->rSon = p_a;

  /* recalculer les champ equilibre */

#define MIN(x,y) (((x)<(y)) ? (x) : (y))
#define MAX(x,y) (((x)<(y)) ? (y) : (x))

  vAnew = vA + 1 - MIN(0,vB);
  vBnew = 1 + MAX(0,vAnew) + vB;

#undef MIN
#undef MAX

  p_a->equilibre = vAnew;
  p_b->equilibre = vBnew;

  return p_b; /* nouvelle racine */

}


void     AVLfree  (PtrAVL p_a) {
  if (p_a != NULL) {
    AVLfree(p_a->lSon);
    p_a->lSon = NULL;
    
    AVLfree(p_a->rSon);
    p_a->rSon = NULL;
    
    free(p_a);
  }
}


AVLPaire    AVLequilibrerAjout (PtrAVL p_a) {

  /* renvoie -1 si la profondeur a diminuee              */
  /* renvoie  0 si la profondeur est inchangee           */
  /* a verifier si OK */
  
  AVLPaire p ;
  
  p.incr = -1;
  if (p_a->equilibre<-1) {
    /* branche gauche trop profonde */

    if (p_a->lSon->equilibre<0) {
      p_a = rRot(p_a); 
      }
    else {
      p_a->lSon = lRot(p_a->lSon);
      p_a = rRot(p_a);  
    }
  }
  else if (p_a->equilibre>1) { 
    /* branche droite trop profonde */
    if (p_a->rSon->equilibre>0) {
      p_a = lRot(p_a); 
    }
    else {  
      p_a->rSon = rRot(p_a->rSon); 
      p_a = lRot(p_a);     
    }
  }
  else {
    /* aucun equilibrage effectué */
    p.incr = 0;
  }
  
  p.p_AVL = p_a;

  return p;
}


AVLPaire    AVLequilibrerSupp (PtrAVL p_a) {

  /* renvoie -1 si la profondeur a diminuee              */
  /* renvoie  0 si la profondeur est inchangee           */
  /* a verifier si OK */
  
  AVLPaire p ;
 
  p.incr = -1;
  if (p_a->equilibre<-1) {
    /* branche gauche trop profonde */
    if (p_a->lSon->equilibre<0) {
      /* branche gauche trop profonde */
      p_a = rRot(p_a);
    
      }
    else {
      if (p_a->lSon->equilibre==0) {
    p_a = rRot(p_a);
    p.incr = 0;
    
      } 
      else {
    p_a->lSon = lRot(p_a->lSon);
    p_a = rRot(p_a);
      }
    }
  }
  else if (p_a->equilibre>1) { 
    /* branche droite trop profonde */ 
    if (p_a->rSon->equilibre>0) {
      /* branche droite trop profonde */
     
      p_a = lRot(p_a); 
    }
    else {   
      if (p_a->rSon->equilibre==0) {
    p_a = lRot(p_a);
    p.incr = 0;

      } 
      else { 
    p_a->rSon = rRot(p_a->rSon); 
    p_a = lRot(p_a);     
      }
    }
  }
  else {
    /* aucun equilibrage effectué */
    p.incr = 0;
  }
  
  p.p_AVL = p_a;
 

  return p;


}

AVLPaire    AVLpush1    (PtrAVL p_a, void* data, int val) {
  /* renvoie 1 si la profondeur a augmentee   */
  /* renvoie 0 si la profondeur est inchangee */


  int incr, r;
  AVLPaire p;
  
  /* ajout */

  if (p_a == NULL) {
    p_a = AVLnew();
    p_a->equilibre = 0;
    p_a->val = val;
    p_a->data = data;
    r = 1;
  }
  else {
    if (p_a->val > val) {
      /* ajout a gauche */
      p = AVLpush1(p_a->lSon,data,val);
      incr = -p.incr;
      p_a->lSon = p.p_AVL;
    }
    else {
      /* ajout a droite */
      p = AVLpush1(p_a->rSon,data,val);
      incr = p.incr;
      p_a->rSon = p.p_AVL;
    }
  
    /* equilibrage */
    p_a->equilibre = p_a->equilibre + incr;

    /* si (incr est différent de 0 ET p_a->equilibre != 0 )         */
    /* alors cela signifie qu'il y a eu un changement de profondeur */
    /* pour le sous arbre concidere                                 */
    
    if (incr != 0 && p_a->equilibre != 0) {
      
      p = AVLequilibrerAjout(p_a);
      
      p_a = p.p_AVL;
      r = 1 + p.incr; 
      /* si on a equilibre on a egalise la hauteur r = 0 */
      /* sinon, on a ajoute un degre de profondeur r = 1 */  
      
    }
    else {
      /* incr == 0 || incr !=0 && p_a->equilibre==0  */
      /* donc la hauteur n'a pas changee */
      r = 0;
    }
  }

  /* retour des resultats : decalage et arbre */
  p.incr = r;
  p.p_AVL = p_a;
  
  return p;
}

void          AVLpush     (PtrAVL *p_p_a, void* data, int val) {
  (*p_p_a) = ((AVLpush1((*p_p_a),data,val)).p_AVL);
}

PtrAVL        AVLgetMin   (PtrAVL p_a) {
  PtrAVL retour;

  if (p_a == NULL) retour = NULL;
  else if (p_a->lSon == NULL) retour = p_a;
  else {
    while (p_a->lSon != NULL)
      p_a = p_a->lSon;
    retour = p_a;
  }
  return retour;
}

PtrAVL        AVLgetMax   (PtrAVL p_a) {
  PtrAVL retour;
  if (p_a == NULL) retour = NULL;
  else if (p_a->rSon == NULL) retour = p_a;
  else {
    while (p_a->rSon != NULL)
      p_a = p_a->rSon;
    retour = p_a;
  }
  return retour;
}

AVLPaire        AVLremoveMin1   (PtrAVL p_a) {
  /* AVL non vide en entree */

  /* renvoie -1 si la profondeur a diminuee   */
  /* renvoie 0 si la profondeur est inchangee */

  AVLPaire p;
  int incr, r;
  PtrAVL tmp;
  
  if (p_a->lSon == NULL ) {
       
    /* ici p_a designe le noeud a enlever */ 
    r = -1;
    tmp = p_a;
    p_a = p_a->rSon;  
   
    free(tmp);
   
  } 
  else {
    /* il faut chercher plus a gauche */
    p = AVLremoveMin1(p_a->lSon); 
       
    incr = -p.incr;
    p_a->lSon = p.p_AVL;
           
    
    
    /* equilibrage */
    p_a->equilibre = p_a->equilibre + incr;
    
    if (incr != 0 && p_a->equilibre != 0) { 
   
      p = AVLequilibrerSupp(p_a);
   
      p_a = p.p_AVL;
      r = p.incr; 
      /* si on a equilibre on a reduit la hauteur r = -1 */
      /* sinon, on a rien change r = 0                   */  
    }
    else {
      /* incr == 0 || incr != 0 && p_a->equilibre==0  */
      /* donc la hauteur n'a pas changee */ 
      r = (incr!=0) ? -1 : 0;
    }
  }        
 
  p.incr = r;
  p.p_AVL = p_a; 
  return p;
}


AVLPaire        AVLremoveMax1   (PtrAVL p_a) {
    /* AVL non vide en entree */

  /* renvoie -1 si la profondeur a diminuee   */
  /* renvoie 0 si la profondeur est inchangee */

  AVLPaire p;
  int incr, r;
  PtrAVL tmp;
  if (p_a->rSon == NULL ) {
    /* ici p_a designe le noeud a enlever */ 
    r = -1;
    tmp = p_a;
    p_a = p_a->lSon;
    free(tmp);
  }
  else {
    /* il faut chercher plus a droite */
    p = AVLremoveMax1(p_a->rSon);
    incr = p.incr;
    p_a->rSon = p.p_AVL;
           
   
    /* equilibrage */
    p_a->equilibre = p_a->equilibre + incr;
    
    if (incr != 0 && p_a->equilibre != 0) {
      p = AVLequilibrerSupp(p_a);
      
      p_a = p.p_AVL;
      r = p.incr; 
      /* si on a equilibre on a reduit la hauteur r = -1 */
      /* sinon, on a rien change r = 0                   */  
    }
    else {
      /* incr == 0 || incr != 0 && p_a->equilibre==0  */
      /* donc la hauteur n'a pas changee */
      r = (incr!=0) ? -1 : 0;
    }
  
  }      
  
  p.incr = r;
  p.p_AVL = p_a;
  
  return p;
}

AVLPaire        AVLremove1   (PtrAVL p_a, int val) {
  /* renvoie -1 si la profondeur a diminuee   */
  /* renvoie 0 si la profondeur est inchangee */

  AVLPaire p;
  int incr, r;
  PtrAVL p_min, p_max;

  if (p_a == NULL) {
    p.incr = 0;
    p.p_AVL = NULL;
  }
  else {
    if (p_a->val > val ) {
      p = AVLremove1(p_a->lSon,val);
      incr = -p.incr;
      p_a->lSon = p.p_AVL;
    } 
    else if (p_a->val < val ) {
      p = AVLremove1(p_a->rSon,val);
      incr = p.incr;
      p_a->rSon = p.p_AVL;
    }
    else {
      /* ici p_a designe le noeud a enlever */
      if (p_a->lSon != NULL) {
    p_max = AVLgetMax(p_a->lSon);
    p_a->data = p_max->data;
    p_a->val = p_max->val;
    p_max = NULL;
    
    p = AVLremoveMax1(p_a->lSon);
    
    incr = -p.incr;
    p_a->lSon = p.p_AVL;
    
      }
      else if (p_a->rSon != NULL) {
    p_min = AVLgetMin(p_a->rSon);
    p_a->data = p_min->data;
    p_a->val = p_min->val;
    p_min = NULL;
    
    p = AVLremoveMin1(p_a->rSon);
    
    incr = p.incr;
    p_a->rSon = p.p_AVL;
      }
      else {
    free(p_a);
    incr = -1;
    p_a = NULL;
      }
    }
    
    /* retour et equilibrage */
    if (p_a == NULL) {
      r = incr;
    }
    else {
      /* equilibrage */
      p_a->equilibre = p_a->equilibre + incr;
      
      if (incr != 0 && p_a->equilibre != 0) {
    p = AVLequilibrerSupp(p_a);
    
    p_a = p.p_AVL;
    r = p.incr;
    /* si on a equilibre on a reduit la hauteur r = -1 */
    /* sinon, on a rien change r = 0                   */  
      }
      else {
    /* incr == 0 || incr != 0 && p_a->equilibre==0  */
    /* donc la hauteur n'a pas changee */
    r = (incr!=0) ? -1 : 0;
      }
    }

    p.incr = r;
    p.p_AVL = p_a;

  }
    return p;
}

void     AVLremove      (PtrAVL* p_p_a, int val) {
  (*p_p_a) = ((AVLremove1((*p_p_a),val)).p_AVL);
}

AVLPaire        AVLremoveData1   (int* p_effectue, PtrAVL p_a, int val, void* data, AVLcompareF) {
  /* renvoie -1 si la profondeur a diminuee   */
  /* renvoie 0 si la profondeur est inchangee */
  
  AVLPaire p;
  int incr, r;
  PtrAVL p_min, p_max;
 
 
  if (p_a == NULL) {
  
    p.incr = 0;
    p.p_AVL = NULL;
  }
  else { 
    
    if (p_a->val > val ) {
      p = AVLremoveData1(p_effectue, p_a->lSon, val, data, AVLcompareFunc);
      incr = -p.incr;
      p_a->lSon = p.p_AVL;
    } 
    else if (p_a->val < val ) { 
      p = AVLremoveData1(p_effectue, p_a->rSon, val, data, AVLcompareFunc);
      incr = p.incr;
      p_a->rSon = p.p_AVL;
    }
    else {
      /* ici p_a designe le noeud a enlever */
      /* seulement si le noeud est egal à data */
    
      if (AVLcompareFunc(p_a->data, data)) {

    (*p_effectue) = 1;
    //fprintf(stderr,"Il dit qu'on va se faire defoncer la gueule !\n");

    if (p_a->lSon != NULL) {
      p_max = AVLgetMax(p_a->lSon);
      p_a->data = p_max->data;
      p_a->val = p_max->val;
      p_max = NULL;
      
      p = AVLremoveMax1(p_a->lSon);
      
      incr = -p.incr;
      p_a->lSon = p.p_AVL;
      
    }
    else if (p_a->rSon != NULL) {
      p_min = AVLgetMin(p_a->rSon);
      p_a->data = p_min->data;
      p_a->val = p_min->val;
      p_min = NULL;
      
      p = AVLremoveMin1(p_a->rSon);
      
      incr = p.incr;
      p_a->rSon = p.p_AVL;
    }
    else {
      free(p_a);
      incr = -1;
      p_a = NULL;
    }
      }
      else {
    /* alors il faut chercher dans le fils droit */
    /* et le fils gauche */

    //fprintf(stderr,"La Main Gauche\n");
    p = AVLremoveData1(p_effectue, p_a->lSon, val, data, AVLcompareFunc);
    
    incr = -p.incr;
    p_a->lSon = p.p_AVL;

    if ((*p_effectue) == 0) {
      //fprintf(stderr,"La Main Droite\n");
      p = AVLremoveData1(p_effectue, p_a->rSon, val, data, AVLcompareFunc);
      
      incr = p.incr;
      p_a->rSon = p.p_AVL;
    }
      }
    }
    
    /* retour et equilibrage */
    if (p_a == NULL) {
      r = incr;
    }
    else {
      /* equilibrage */
      p_a->equilibre = p_a->equilibre + incr;
      
      if (incr != 0 && p_a->equilibre != 0) {
    p = AVLequilibrerSupp(p_a);
    
    p_a = p.p_AVL;
    r = p.incr;
    /* si on a equilibre on a reduit la hauteur r = -1 */
    /* sinon, on a rien change r = 0                   */  
      }
      else {
    /* incr == 0 || incr != 0 && p_a->equilibre==0  */
    /* donc la hauteur n'a pas changee */
    r = (incr!=0) ? -1 : 0;
      }
    }

    p.incr = r;
    p.p_AVL = p_a;

  }
    return p;
}

void     AVLremoveData   (PtrAVL* p_p_a, int val, void* data, int (*compFunc) (const void* ,const void*)) {
  int effectue = 0;

  fprintf(stderr,"appel AVLremoveData1\n");
  (*p_p_a) = ((AVLremoveData1(&effectue,(*p_p_a),val, data, compFunc)).p_AVL); 
  fprintf(stderr,"retour AVLremoveData1\n");
}

void*    AVLget      (PtrAVL p_a, int val) {
  
  void* retour;

  if (p_a == NULL) {
    retour = NULL;
  }
  else {
    if (p_a->val > val) {
      retour = AVLget(p_a->lSon,val);
    }
    else if (p_a->val < val){
      retour = AVLget(p_a->rSon,val);
    }
    else {
      retour = p_a->data;
    }
  }
  
  return retour;
}

void*    AVLpop      (PtrAVL* p_p_a, int val) {
  void* retour;

  retour = AVLget((*p_p_a),val);
  AVLremove(p_p_a, val);

  return retour;
}

void*    AVLpopMin  (PtrAVL *p_p_a) {
  void* retour;
  PtrAVL p;
  
  p = AVLgetMin((*p_p_a));
  retour = p->data;
  
  (*p_p_a) = (AVLremoveMin1(*p_p_a)).p_AVL;
  
 
  return retour;
}

void*    AVLpopMax  (PtrAVL *p_p_a) {
  void* retour;

  retour = AVLgetMax((*p_p_a))->data;
  (*p_p_a) = (AVLremoveMax1(*p_p_a)).p_AVL;

  return retour;
}


void     AVLevalPostfixe (PtrAVL p_a, evalF) {
  
  if (p_a != NULL) {
    AVLevalPostfixe(p_a->lSon, evalFunc);
    AVLevalPostfixe(p_a->rSon, evalFunc); 
    evalFunc(p_a);
  }
}





int compint (const void* i1,const void* i2) {
  fprintf(stderr,"test Mortel retourne : %d == %d ? = %d\n",((int)i1),((int)i2),(((int)i1) - ((int)i2))==0);
  return (((int)i1) - ((int)i2))==0; 
}

int main (int argc, char** argv) {
  
  PtrAVL p_a=NULL;
  int i,nb,valPrec;
  int val;
  int SEED = 1;
  int NBINAVL = 100000;
  int MAXVAL = 10;
  
  srandom(SEED);

  for (i=0;i<NBINAVL;i++) {
    val=(((int)random())%(MAXVAL-1))+1;
    AVLpush(&p_a,(void*) val,val);
  }
  AVLpush(&p_a,(void*) val,val);
  fprintf(stderr,"%6d\n",i+1);

  fprintf(stderr,"\n");

  
  for(i=0;i<5000;i++) {
    valPrec = (((int)random())%(MAXVAL-1))+1;
    val = (int) AVLpop(&p_a,valPrec);
    fprintf(stderr,"On essaye d'enlever %4d, on obtient %4d\n",valPrec, val);
  }
  
  for(i=0;i<500000;i++) {
    val = (int) AVLpop(&p_a,valPrec);
    if (val != NULL) 
        fprintf(stderr,"On essaye d'enlever %4d, on obtient %4d (%4d)\r",valPrec, val, i);
    else {
        fprintf(stderr,"\nOn essaye d'enlever %4d, on obtient %4d\n",valPrec, val);
        break;
    }
  }
  nb=0;
  valPrec = 0;
  while (p_a!=NULL) {
    val =  (int)AVLpopMin(&p_a);  

    fprintf(stderr,"%6d\r",nb);
//    fprintf(stderr,"%d\n",val);
    nb++;
    if (val<valPrec) {
      fprintf(stderr,"Boum !\n"); 
      exit(0);
    }
    valPrec = val;
  }
  fprintf(stderr,"\n");

  fprintf(stderr,"*********************\n");
  fprintf(stderr,"%6d %6d\n",NBINAVL,nb);
  

  return (NBINAVL==nb);


} 


