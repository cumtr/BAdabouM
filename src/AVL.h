#ifndef _H_AVLH 
#define _H_AVLH

#define AVLcompareF int (*AVLcompareFunc)(const void *,\
                                          const void *)
#define evalF void (*evalFunc)(PtrAVL p_a)

struct AVL;
typedef struct AVL *PtrAVL;

struct AVL {   
  int val;
  int equilibre;
  void* data;   
  PtrAVL lSon;
  PtrAVL rSon;
};


void     AVLafficher     (PtrAVL p_a);
void     freeAVL         (PtrAVL p_a);
void     AVLpush         (PtrAVL *p_p_a, void* data, int val);
void     AVLremove       (PtrAVL *p_p_a, int val);
void     AVLremoveData   (PtrAVL* p_p_a, int val, void* data, AVLcompareF);
void*    AVLget          (PtrAVL p_a, int val);
PtrAVL   AVLgetMin       (PtrAVL p_a);
void*    AVLpop          (PtrAVL *p_p_a, int val);
void*    AVLpopMin       (PtrAVL *p_p_a);
void*    AVLpopMax       (PtrAVL *p_p_a);
void     AVLevalPostfixe (PtrAVL p_a, evalF);


#endif
