#include "general.h"

#if defined(_WIN32) || defined(WIN32)
# include <Windows.h>
#else
# include <unistd.h>
#endif

/*Define Boolean type */
#if !defined(__cplusplus)
  typedef enum bool { false, true } bool;
#endif


#define SORTCLEFTSBY num_spheres
#define SIZE_TEMP_LINES 2000

typedef struct tCleftStructure tsCleft;
typedef tsCleft *tCleft;

typedef struct tSphereStructure tsSphere;
typedef tsSphere *tSphere;

typedef struct tAtomStructure tsAtom;
typedef tsAtom *tAtom;

typedef struct tResStructure tsRes;
typedef tsRes *tRes;

struct tCleftStructure{
  tSphere start,end;
  float   volume;
  //int   num_gpoints;
  int     label;
  int     id;
  int     num_spheres;
  float   center[3];
  float   effrad;
  int     toprint;
  tCleft  next,prev;
};

struct tSphereStructure{
  float    center[3];
  float    radius;
  int      inum;
  //int    num_gpoints;
  tCleft   cleft;
  tSphere  next,prev;
};

struct tAtomStructure{
  char    line[82];   // PDB line
  char    name[5];    // atom name
  float   radius;     // atomic radius
  float   coor[3];    // coordinates
  tSphere ofSphere;   // pointer to sphere to which atom belongs
  int     isprot;     // 1 if ATOM - 0 if HETATM
  int     inum;
  tAtom   next,prev;  // connects the list of all atoms
  tAtom   next_inres; // link to next in list of all atoms in the same residue
  tAtom   prev_inres; // link to prev in list of all atoms in the same residue
  tRes    ofres;
  int     out;        // 1 if atom is to be outputed, 0 otherwise
};

struct tResStructure{
  char   name[4];     // residue name
  int    num;         // residue number
  char   chain;       // residue chain
  char   alt;
  tAtom  atom_anchor; // first atom read found to be part of this residue
  int    cleft_anchor;// 1 if residue should mark cleft to be outputed
  int    hetero;      // 0 for aminoacids, 1 for hetero groups
  tRes   next,prev;   // link to other residues in the list
  tCleft ofcleft;
  int    out;         // 1 if residue contains atoms to be outputed;
};


tAtom   atoms   = NULL;
tSphere spheres = NULL;
tCleft  clefts  = NULL;
tRes    resids  = NULL;

int   num_spheres=0;
int   top_clefts=0;

char  pdb_base[100];
char  out_base[100];

int   output_het;
int   het_whole;
int   complete_residue;
int   calpha;
int   cbeta;
int   output_extra_atoms;
int   output_spheres;

char  anchor_nam[4];
char  anchor_nam_copy[4];
int   anchor_num;
char  anchor_chn;
int   anchor_flg;
char  anchor_alt;

float contact_threshold;
int   output_im_cleft;
int   output_omim_clefts;

char  pdb_file[100];
char  chn[10];
int   chn_counter;
double pi;

float  sphere_lwb;
float  sphere_upb;

void    get_clefts();
float   assign_radius(char atm[]);
tSphere RotateSphereList(int n);
tCleft  MakeNullCleft(void);
void    AddNewSphere(float center[], float radius);
void    read_pdb(char file[]);
void    merge_clefts();
void    print_cleft(int label);
void    swap_clefts(tCleft b, tCleft d);
void    sort_clefts();
void    assign_atoms_to_clefts();
void    output_atoms_in_cleft(tCleft c);
void    output_atoms_in_interaction_model_cleft(tCleft c, tRes residue);
void    output_spheres_in_cleft(tCleft c);

tRes    assign_residue(char rnam[],int rnum, char chn, char alt, tAtom a);
void    get_contacting_chains(char filename[]);
void    read_commandline(int argc, char *argv[]);
float   dist3d(float a[], float b[]);
void    print_cleft_surface(char filename[], tCleft c);
/***********************************************************************/
/*        1         2         3         4         5         6         7*/
/*234567890123456789012345678901234567890123456789012345678901234567890*/
/*        1         2         3         4         5         6         7*/
/***********************************************************************/
int main(int argc, char *argv[]){

  tCleft c;
  tRes res;
  int id=1;
  //tRes  tmpresids=NULL;
  //tAtom tmpatoms=NULL;

  pi=acos(-1.0);

  // read command line arguments
  read_commandline(argc, argv);
  

  if(anchor_flg == 1) get_contacting_chains(pdb_file);

  // load PDB atoms
  read_pdb(pdb_file);

  /*
    a=atoms;
    do{
    printf("%s\n",a->line);
    a=a->next;
    }while(a != atoms);
    PAUSE;
  */
  
  res=resids;
  do{
    if(strcmp(res->name,anchor_nam) == 0 &&
       res->num == anchor_num &&
       res->chain == anchor_chn) res->cleft_anchor=1;
    //printf("res->(name=<%s> num=<%d> chain=<%c> anchor=%d) [%s %d %c]\n",res->name,
    //  res->num,res->chain,res->cleft_anchor,anchor_nam,
    // anchor_num,anchor_chn);
    //Pause;
    //if(res->cleft_anchor == 1) PAUSE;
    res=res->next;
  }while(res != resids);
  //PAUSE;
  

  //printf("%d: %s\n",atoms->inum,atoms->line);
  //printf("%d: %s\n",atoms->prev->inum,atoms->prev->line);
  //printf("%d: %s\n",atoms->prev->next->inum,atoms->prev->next->line);
  //PAUSE;

  // calculate clefts
  get_clefts();

  assign_atoms_to_clefts();
  
  //printf("all right so far!\n");
  //PAUSE;
  
  // assign an incremental id to clefts
  c=clefts;
  do{
	  c->id = id++;
	  c=c->next;
  }while(c != clefts);
  
  if(top_clefts > 0){
    c=clefts;

    while(c->label <= top_clefts){
      //printf("cleft label=%d\n",c->label); 
      output_atoms_in_cleft(c);
      if(output_spheres == 1) output_spheres_in_cleft(c);
      c=c->next;
      
      if(c->label==1){
	printf("there are a total of %d Clefts\n",c->prev->label);
	break;
      }
    }
  }else{
    //printf("got here and anchor_flg = %d\n",anchor_flg);
    //PAUSE;
    if(anchor_flg == 1){
      res=resids;
      do{
	//printf("res->(name=<%s> num=<%d> chain=<%c>",res->name,
	//       res->num,res->chain);
	//printf("anchor=%d ofcleft=%d) [%s %d %c]\n",res->cleft_anchor,res->ofcleft->label,anchor_nam,
	//       anchor_num,anchor_chn);

	if(res->cleft_anchor == 1){
	  //printf("residue=%s%d%c",res->name,res->num,res->chain);PAUSE;
	  
	  if(res->ofcleft == NULL){
	    fprintf(stderr,"Residue %s%d%c is not apart a cleft.\n",res->name,res->num,res->chain);
	    exit(1);
	  }

	  //printf("ofcleft=%d\n",res->ofcleft);PAUSE;
	  //printf("ofcleft=%d\n",res->ofcleft->label);PAUSE;

	  //if(output_im_cleft==0 || output_omim_clefts==1){
	  //  output_atoms_in_cleft(res->ofcleft);
	  //  if(output_spheres == 1) output_spheres_in_cleft(res->ofcleft);
	  //}
	  //if(output_im_cleft==1 || output_omim_clefts==1){
	  //  output_atoms_in_interaction_model_cleft(res->ofcleft,res);	    
	  //}
	  // -b
	  if (output_omim_clefts==1){            
	    output_atoms_in_cleft(res->ofcleft);
	    if(output_spheres == 1) output_spheres_in_cleft(res->ofcleft);
	    output_atoms_in_interaction_model_cleft(res->ofcleft,res);
	  }
	  // -a
	  else if (output_im_cleft==0){
	    output_atoms_in_cleft(res->ofcleft);
	    if(output_spheres == 1) output_spheres_in_cleft(res->ofcleft);
	  }
	  // -i
	  else if (output_im_cleft==1){
	    output_atoms_in_interaction_model_cleft(res->ofcleft,res);	    
	  }
	}
	//printf("res=%3s%d%c\n",res->name,res->num,res->chain);
	res=res->next;
	//printf("res=%3s%d%c\n",res->name,res->num,res->chain);

      }while(res != resids);
    }else{
      c=clefts;
      printf("there are a total of %d Clefts\n",clefts->prev->label);
      do{
	//printf("cleft label=%d\n",c->label); 
	output_atoms_in_cleft(c);
	if(output_spheres == 1) output_spheres_in_cleft(c);
	c=c->next;
      }while(c != clefts);
    }
  }

  return 0;
}
/***********************************************************************/
/*        1         2         3         4         5         6         7*/
/*234567890123456789012345678901234567890123456789012345678901234567890*/
/*        1         2         3         4         5         6         7*/
/***********************************************************************/
void get_contacting_chains(char filename[]){
  FILE  *infile_ptr;        // pointer to input file
  char  buffer[81];         // a line from the INPUT file
  char  coor_char[10];      // string used to read the coordinates
  char  field[7];
  int   i,j;

  char  rnam[4];
  int   rnum;
  char  chain;
  float coor[3];
  tAtom a;
  tAtom lig_tmp=NULL;

  int  chn_c,chn_flag;

  infile_ptr = fopen(filename, "r");

  if(infile_ptr == NULL){
    fprintf(stderr, "ERROR: unable to open input file: %s\n",filename);
    exit(8);
  }
  while (fgets(buffer, sizeof(buffer),infile_ptr)){
    for(i=0;i<6;i++) field[i]=buffer[i];
    field[6]='\0';
    if(strcmp(field,"ATOM  ") == 0 || strcmp(field,"HETATM") == 0){
      for(i=0;i<3;i++) rnam[i]=buffer[i+17];
      rnam[3]='\0';
      for(i=0;i<4;i++){coor_char[i]=buffer[22+i];}
      coor_char[4]='\0';
      sscanf(coor_char,"%d",&rnum);
      chain=buffer[21];
      for (j=0;j<=2;j++){
	for(i=0;i<=7;i++){coor_char[i]=buffer[30+j*8+i];}
	coor_char[8]='\0';
	sscanf(coor_char,"%f",&coor[j]);
      }
      if(strcmp(rnam,anchor_nam) == 0 && 
	 rnum == anchor_num &&
	 chain == anchor_chn){
	
	NEW( a, tsAtom );
	ADD( lig_tmp, a );
	for(i=0;i<3;i++) a->coor[i] = coor[i];
      }
    }
  }
  rewind(infile_ptr);    

  chn_counter=-1;

  while (fgets(buffer, sizeof(buffer),infile_ptr)){
    for(i=0;i<6;i++) field[i]=buffer[i];
    field[6]='\0';
    
    if(strcmp(field,"ATOM  ") == 0){
      for(i=0;i<3;i++) rnam[i]=buffer[i+17];
      rnam[3]='\0';
      for(i=0;i<4;i++){coor_char[i]=buffer[22+i];}
      coor_char[4]='\0';
      sscanf(coor_char,"%d",&rnum);
      chain=buffer[21];
      for (j=0;j<=2;j++){
	for(i=0;i<=7;i++){coor_char[i]=buffer[30+j*8+i];}
	coor_char[8]='\0';
	sscanf(coor_char,"%f",&coor[j]);
      }
      if(strcmp(rnam,anchor_nam) != 0 || 
	 rnum != anchor_num ||
	 chain != anchor_chn){
	
	a=lig_tmp;
	do{
	  if(dist3d(a->coor,coor) < contact_threshold){
	    chn_flag=0;
	    if(chn_counter > -1){
	      chn_c=0;
	      while(chn_c  <= chn_counter && chn_flag == 0){
		if(chain == chn[chn_c]) chn_flag=1;
		chn_c++;
	      }
	    }
	    if(chn_flag == 0){
	      chn_counter++;
	      chn[chn_counter]=chain;
	    }
	  }
	  
	  a=a->next;
	}while(a != lig_tmp);
      }
    }
  }
  fclose(infile_ptr);    

  /* 
     for(chn_c=0;chn_c <= chn_counter;chn_c++){
     printf("chn=<%c> ",chn[chn_c]);
     }
     PAUSE;
  */

  return;
}
/***********************************************************************/
/*        1         2         3         4         5         6         7*/
/*234567890123456789012345678901234567890123456789012345678901234567890*/
/*        1         2         3         4         5         6         7*/
/***********************************************************************/
//**************** NEEDS TO BE REMOVED !!!!!! **************************
void print_atom_contacts(tRes residue){
  tRes r;
  tAtom a,b;
  FILE *outfile_ptr;
  char filename[100];
  char chn,alt;
  int flag;
  int  print_residue_flag;
  char temp_lines[SIZE_TEMP_LINES];

  
  if(residue->chain == ' '){
    chn='-';}
  else { 
    chn=residue->chain; }
 
  if(residue->alt == ' '){
    alt='-';}
  else { 
    alt=residue->alt; }

  //sprintf(filename,"%s_clf_%d_%s%d%c%c_env.pdb",out_base,residue->ofcleft->label,residue->name,residue->num,chn,alt);
  sprintf(filename,"%s_clf_%d_env.pdb",out_base,residue->ofcleft->label);
  //printf("filename for atoms in contact=%s\n",filename);
  //PAUSE;

  outfile_ptr = fopen(filename, "w");
  if(outfile_ptr == NULL){
    fprintf(stderr, "ERROR: unable to open output file %s for write\n",filename);
    exit(8);
  }

  if(output_extra_atoms == 1){
    if(complete_residue == 1){
      r=resids;
      do{
	if(r->hetero == 0){
	  temp_lines[0]='\0';
	  print_residue_flag=0;
	  a=r->atom_anchor;
	  do{
	    strcat(temp_lines,a->line);
	    strcat(temp_lines,"\n");
	    if(strlen(temp_lines) > SIZE_TEMP_LINES){
		    printf("IO ERROR: temp_lines size exeeded (%d)!\n",(int)strlen(temp_lines));
	      exit(8);
	    }
	    
	    b=residue->atom_anchor;
	    flag=0;
	    do{
	      if(flag == 0 && dist3d(a->coor,b->coor) < contact_threshold) flag=1;	      
	      b=b->next_inres;
	    }while(b != residue->atom_anchor->prev_inres);

	    if(flag == 1) print_residue_flag=1;

	    a=a->next_inres;
	  }while(a != r->atom_anchor);

	  if(print_residue_flag==1) {fprintf(outfile_ptr,"%s",temp_lines);}
	  
	}
      	r=r->next;
      }while(r != resids);
    }else{
      fprintf(stderr, "Wrong choice of parameters: no output!\n");
      exit(8);
    }
  }else{
    a=atoms;
    do{      
      if(a->isprot != 0){
	flag=0;
	b=residue->atom_anchor;
	do{
	  if(a->ofres != b->ofres && flag == 0){
	    if(dist3d(a->coor,b->coor) < contact_threshold){
	      fprintf(outfile_ptr,"%s\n",a->line);
	      flag=1;
	    }
	  }
	  b=b->next_inres;
	}while(b != residue->atom_anchor->prev_inres);
      }
      a=a->next;
    }while(a != atoms);
  }
  fclose(outfile_ptr);
  
  return;
}
/***********************************************************************/
/*        1         2         3         4         5         6         7*/
/*234567890123456789012345678901234567890123456789012345678901234567890*/
/*        1         2         3         4         5         6         7*/
/***********************************************************************/
void output_spheres_in_cleft(tCleft c){
  tSphere s;
  FILE *outfile_ptr;
  char filename[100];
  int j;
  //int pdb=1;
  
  //sprintf(filename,"%s_sph_%d.pdb",out_base,c->label);
  sprintf(filename,"%s_sph_%d.pdb",out_base,c->id);
  outfile_ptr = fopen(filename, "w");
  if(outfile_ptr == NULL){
    fprintf(stderr, "ERROR: unable to open output file %s for write\n",filename);
    exit(8);
  }
  // if(pdb == 0){
  //  fprintf(outfile_ptr,"radius  x-coor    y-coor   z-coor\n");
  //    s=c->start;
  //  do{
  //    fprintf(outfile_ptr,"%5.3f ",s->radius);
  //    for(i=0;i<3;i++) fprintf(outfile_ptr,"%8.3f ",s->center[i]);
  //    fprintf(outfile_ptr,"\n");
  //     s=s->next;
  //  }while(s != c->end->next);
  //}
  //else{
    s=c->start;

    fprintf(outfile_ptr,"#REMARK  PARENTFILE  %s\n", pdb_file);
    
    do{    
      fprintf(outfile_ptr,"ATOM  %5d  C   SPH Z   1    ",s->inum);
      for(j=0;j<3;j++) fprintf(outfile_ptr,"%8.3f",s->center[j]);
      fprintf(outfile_ptr,"  1.00  %3.2f ",s->radius);
      fprintf(outfile_ptr,"\n");
      s=s->next;
    }while(s != c->end->next);
    // }
  return;
}
/***********************************************************************/
/*        1         2         3         4         5         6         7*/
/*234567890123456789012345678901234567890123456789012345678901234567890*/
/*        1         2         3         4         5         6         7*/
/***********************************************************************/
void output_atoms_in_interaction_model_cleft(tCleft c, tRes residue){
  tAtom a,b;
  tRes r;
  //float d_thres=0.30f;
  FILE *outfile_ptr;
  char filename[100];
  char temp_lines[SIZE_TEMP_LINES];

  temp_lines[0]='\0';

  //sprintf(filename,"%s_clf_%d_IM.pdb",out_base,c->label);
  sprintf(filename,"%s_clf_%d_IM.pdb",out_base,c->id);
  //printf("filename for IM cleft: %s\n",filename);
  //PAUSE;
  outfile_ptr = fopen(filename, "w");
  if(outfile_ptr == NULL){
    fprintf(stderr, "ERROR: unable to open output file %s for write\n",filename);
    exit(8);
  }

  //printf("output_extra_atoms = %d\n",output_extra_atoms);
  //PAUSE;

  // reset output atoms and residues
  a=atoms;
  do{
    a->out=0;
    a->ofres->out=0;

    a=a->next;
  }while(a != atoms);

  a=atoms;
  do{    
    if(a->ofSphere != NULL && a->ofSphere->cleft == c && a->isprot == 1){
      
      b=residue->atom_anchor;
      do{
	if(dist3d(a->coor,b->coor) < contact_threshold){
	  a->out=1;
	  a->ofres->out=1;
	  //printf("%s\n",a->line);
	  //printf("%s %c %d %d\n",a->ofres->name,
	  //a->ofres->chain,a->ofres->num,a->ofres->out);
	}
	b=b->next_inres;
      }while(a->out == 0 && b != residue->atom_anchor->prev_inres);
      
    }
    a=a->next;
  }while(a != atoms);
  
  //PAUSE;

  if(output_extra_atoms == 1){
    r=resids;
    do{
      //printf("%s %c %d ",r->name,r->chain,r->num);
      if(r->out == 1){
	//printf("\n");
	a=r->atom_anchor;
	do{
	  //printf("\t%s :: ",a->line);
	  if(complete_residue == 1){
	    a->out=1;
	  }else{
	    if(calpha == 1 && strcmp(a->name," CA ") == 0) a->out=1;
	    if(cbeta  == 1 && strcmp(a->name," CB ") == 0) a->out=1;
	  }
	  //printf("%d\n",a->out);
	  //PAUSE;
	  a=a->next_inres;
	}while(a != r->atom_anchor);
      }//else{printf("skip\n");}
      
      r=r->next;
    }while(r != resids);
  }
  //PAUSE;

  fprintf(outfile_ptr,"#REMARK  PARENTFILE  %s\n", pdb_file);

  a=atoms;
  do{
    if(a->out == 1) fprintf(outfile_ptr,"%s\n",a->line);
    //if(a->out == 1) printf("%s\n",a->line);
    a=a->next;
  }while(a != atoms);
  fclose(outfile_ptr);

  return;
}
/***********************************************************************/
/*        1         2         3         4         5         6         7*/
/*234567890123456789012345678901234567890123456789012345678901234567890*/
/*        1         2         3         4         5         6         7*/
/***********************************************************************/
void output_atoms_in_cleft(tCleft c){
  tAtom a;
  tRes r;
  //float d_thres=0.30f;
  FILE *outfile_ptr;
  char filename[100];
  char temp_lines[SIZE_TEMP_LINES];

  temp_lines[0]='\0';

  //sprintf(filename,"%s_clf_%d.pdb",out_base,c->label);
  sprintf(filename,"%s_clf_%d.pdb",out_base,c->id);
  outfile_ptr = fopen(filename, "w");
  if(outfile_ptr == NULL){
    fprintf(stderr, "ERROR: unable to open output file %s for write\n",filename);
    exit(8);
  }

  // reset output atoms
  a=atoms;
  do{
    a->out=0;
    a->ofres->out=0;
    
    a=a->next;
  }while(a != atoms);
  
  
  a=atoms;
  do{
    if(a->ofSphere != NULL && a->ofSphere->cleft == c && a->isprot == 1){
      a->out=1;  
      if(output_extra_atoms == 1) { a->ofres->out=1; }
    }
    a=a->next;
  }while(a != atoms);


  r=resids;
  do{
    if(r->out == 1){
      a=r->atom_anchor;

      do{
        if((complete_residue == 1)||
           (calpha == 1 && strcmp(a->name," CA ") == 0)||
           (cbeta  == 1 && strcmp(a->name," CB ") == 0))  a->out=1;

	a=a->next_inres;
      }while(a != r->atom_anchor);
    }
    r=r->next;
  }while(r != resids);
  
  fprintf(outfile_ptr,"#REMARK  PARENTFILE  %s\n", pdb_file);

  a=atoms;
  do{
    if(a->out == 1) fprintf(outfile_ptr,"%s\n",a->line);
    a=a->next;
  }while(a != atoms);
  fclose(outfile_ptr);

  if(output_het == 1){
    sprintf(filename,"%s_clf_%d_het.pdb",out_base,c->label);
    outfile_ptr = fopen(filename, "w");
    if(outfile_ptr == NULL){
      fprintf(stderr, "ERROR: unable to open output file %s for write\n",filename);
      exit(8);
    }
    
    if(het_whole == 0){
      a=atoms;
      do{
	if(a->ofSphere != NULL && 
	   a->ofSphere->cleft == c && a->isprot == 0) fprintf(outfile_ptr,"%s\n",a->line);
	a=a->next;
      }while(a != atoms);
    }else{
      r=resids;
      do{
	if(r->hetero == 1 && r->ofcleft == c){
	  a=r->atom_anchor;
	  do{
	    fprintf(outfile_ptr,"%s\n",a->line);
	    a=a->next_inres;
	  }while(a != r->atom_anchor);
	}
	r=r->next;
      }while(r != resids);
    }
    fclose(outfile_ptr);
  }
  
  //printf("so far ok inside sub end\n");
  //PAUSE;

  
  return;
}
/***********************************************************************/
/*        1         2         3         4         5         6         7*/
/*234567890123456789012345678901234567890123456789012345678901234567890*/
/*        1         2         3         4         5         6         7*/
/***********************************************************************/
void assign_atoms_to_clefts(){
  tAtom a;
  tSphere s;
  tRes r;
  tCleft c;
  float d;
  float d_thres=0.30f;

  r=resids;
  do{
    r->ofcleft=NULL;
    r=r->next;
  }while(r != resids);
  
  a=atoms;
  do{
    a->ofSphere=NULL;
    a=a->next;
  }while(a != atoms);
  
  c=clefts;

  if(clefts->start == NULL){
    fprintf(stderr,"No clefts found.\n");
    exit(1);
  }

  do{
    a=atoms;
    do{
      if(a->ofSphere == NULL){
	s=c->start;
	do{
	  //printf("x=%8.3f\ty=%8.3f\tz=%8.3f\n",a->coor[0],a->coor[1],a->coor[2]);
	  //printf("x=%8.3f\ty=%8.3f\tz=%8.3f\n",s->center[0],s->center[1],s->center[2]);
	  d=dist3d(a->coor,s->center) - (s->radius + a->radius);
	  if(d < d_thres){
	    a->ofSphere=s;
	    a->ofres->ofcleft=c;
	    //printf("res->(name=<%s> num=<%d> chain=<%c> ofcleft=%d)\n",
	    //   a->ofres->name,a->ofres->num,a->ofres->chain,a->ofres->ofcleft->label);
	    //PAUSE;
	  }
	  s=s->next;
	}while(a->ofSphere == NULL && s != c->end->next);
      }
      a=a->next;
    }while(a != atoms);
    c=c->next;
  }while(c != clefts);

  
  //a=atoms;
  //do{
  //  printf("ares->(name=<%s> num=<%d> chain=<%c>",a->ofres->name,a->ofres->num,a->ofres->chain);
  //  printf(" ofcleft=%d)\n",a->ofres->ofcleft->label);
  //  a=a->next;
  //}while(a != atoms);

  /*
    r=resids;
    do{
    printf("res->(name=<%s> num=<%d> chain=<%c> ofcleft=",r->name,r->num,r->chain);
    if(r->ofcleft != NULL){
    printf("%d)\n",r->ofcleft->label);
    }else{
    printf("NULL)\n");
    }
    PAUSE;
    r=r->next;
    }while(r != resids);
  */

  return;
}
/***********************************************************************/
/*        1         2         3         4         5         6         7*/
/*234567890123456789012345678901234567890123456789012345678901234567890*/
/*        1         2         3         4         5         6         7*/
/***********************************************************************/
void get_clefts(){
  tAtom   a,b,c;
  int i;
  float maxrad,effrad;
  float center[3];
  //int   num_clefts=0;

  a=atoms;
  do{
    if(a->isprot == 1){
      b=a->next;
      do{
	if(b->isprot == 1){
	  maxrad=(dist3d(a->coor,b->coor)-(a->radius + b->radius))/2;
	  if(maxrad <= sphere_upb){
	    for(i=0;i<3;i++){center[i]=(a->coor[i] + b->coor[i])/2;}
	    c=b->next;
	    do{
	      if(c != a && c->isprot == 1){
		effrad=dist3d(center,c->coor) - c->radius;
		if(maxrad > effrad){maxrad=effrad;}
		if(maxrad < sphere_lwb){c=b->prev;}
	      }
	      c=c->next;
	    }while(c != b);
	    if(maxrad >= sphere_lwb){
	      //printf("adding new sphere\n");
	      AddNewSphere(center,maxrad);
	      num_spheres++;
	    }
	  }
	}
	b=b->next;
      }while(b != atoms);
    }
    a=a->next;
  }while(a != atoms->prev);
  
  //printf("Initial number of spheres: %d\n",num_spheres);
  
  if(clefts==NULL){
    fprintf(stderr,"No clefts found.");
    exit(1);
  }
  merge_clefts();
  
  //printf("clefts merged\n");
  
  /*
    s=spheres;
    do{
    for(i=0;i<3;i++){
    if(pdb_min[i] > s->center[i] - s->radius)	
    pdb_min[i] = s->center[i] - s->radius;
    if(pdb_max[i] < s->center[i] + s->radius)	
    pdb_max[i] = s->center[i] + s->radius;
    }    
    s=s->next;
    }while(s != spheres);
  */

#ifdef SORTCODE
  volume_of_clefts();
#endif

  sort_clefts();
  
  /*
    f=clefts;
    do{
    printf("cleft %3d: volume=%9.3f spheres=%4d\n",
    f->label,f->volume,f->num_spheres);
    f=f->next;
    }while(f != clefts);
  */

  return;
}
/***********************************************************************/
/*        1         2         3         4         5         6         7*/
/*234567890123456789012345678901234567890123456789012345678901234567890*/
/*        1         2         3         4         5         6         7*/
/***********************************************************************/
void sort_clefts(){
  tCleft f,m,x;
  
  // Shift start of circular list to top cleft
  m=clefts;
  f=clefts->next;
  do{
    if(f->SORTCLEFTSBY > m->SORTCLEFTSBY) m=f;
    f=f->next;
  }while(f != clefts);
  clefts=m;
  clefts->label=1;
  
  // sort rest by each time finding the maximum in remaining of list and
  // swapping the maximum with current point if maximum > current
  //z=clefts;
  //do{
  //z=z->next;
  //}while(z != clefts);
  
  f=clefts->next;
  do{
    x=f;
    m=f->next;
    do{
      if(m->SORTCLEFTSBY > x->SORTCLEFTSBY) x=m;
      m=m->next;
    }while(m != clefts);
    
    if(x != f){
      swap_clefts(f,x);
      f=x;
    }
       
    f->label=f->prev->label + 1;
    f=f->next;
  }while(f != clefts->prev);
  f->label=f->prev->label + 1;

  return;
}
/***********************************************************************/
/*        1         2         3         4         5         6         7*/
/*234567890123456789012345678901234567890123456789012345678901234567890*/
/*        1         2         3         4         5         6         7*/
/***********************************************************************/
void swap_clefts(tCleft b, tCleft d){
  tCleft dprev,dnext,f;
  bool  consec=false;
  
  f=d;
  do{
	if(f == b){
	  b=d;
	  d=f;
	  break;
	} 
	f=f->next;
  }while(f != clefts);


  if(b->next == d || b->prev == d){consec=true;}

  dprev=d->prev;
  dnext=d->next;

  d->next->prev = d->prev;
  d->prev->next = d->next;

  if(consec==false){
	b->next->prev = b->prev;
	b->prev->next = b->next;
	
	d->prev=b->prev;
	d->next=b->next;
	b->prev->next=d;
	b->next->prev=d;
	
	dprev->next=b;
	dnext->prev=b;
	
	b->next=dnext;
	b->prev=dprev;
  }else{
	b->prev->next=d;
	b->prev=d;
	d->prev=b->prev;
	d->next=b;
  }

  return;
}
/***********************************************************************/
/*        1         2         3         4         5         6         7*/
/*234567890123456789012345678901234567890123456789012345678901234567890*/
/*        1         2         3         4         5         6         7*/
/***********************************************************************/
void merge_clefts(){
    tCleft  c,d;
    tSphere r,s,t;
    bool connected;
	// twe

    c=clefts;
    do{
        connected=false;
        //PAUSE;
        d=c->next;
        do{
            r=c->start;
            do{
                s=d->start;
                do{
                    if(dist3d(r->center,s->center) <= r->radius+s->radius){
                        connected=true;

                        // 0. reassign d spheres to new cleft c
                        t=d->start;
                        do{
                            t->cleft=c;
                            t=t->next;
                        }while(t != d->end->next);

                        // 1. add d spheres to c
                        // re-routing sphere connections
                        d->end->next->prev=d->start->prev;
                        d->start->prev->next=d->end->next;
                        c->end->next->prev=d->end;
                        d->end->next=c->end->next;
                        d->start->prev=c->end;
                        c->end->next=d->start;
                        c->end=d->end;
        
                        // 2. updating cleft volume and number of spheres
                        //c->volume += d->volume;
                        c->num_spheres += d->num_spheres;
                        d->next->prev=d->prev;
                        d->prev->next=d->next;

			d->prev = NULL;
			d->next = NULL;
                        //free(d); 
                    }
                    s=s->next;

                }while(s != d->end->next && !connected);
                r=r->next;

            }while(r != c->end->next && !connected);

            if(!connected){
                d=d->next;
            }else{
                d=c->next;
                connected=false; 
            }
        }while(d != c);

        c=c->next;
    }while(c != clefts->prev);

    return;
}
/***********************************************************************/
/*        1         2         3         4         5         6         7*/
/*234567890123456789012345678901234567890123456789012345678901234567890*/
/*        1         2         3         4         5         6         7*/
/***********************************************************************/
void AddNewSphere(float center[], float radius){
   tSphere s,t;
   tCleft  c;
   bool connect=false;
   int i;

   NEW(s, tsSphere);
   s->inum=num_spheres;

   if(!clefts){             // no clefts yet
	 //printf("first cleft and sphere\n");
	 ADD(spheres, s);       // 1st sphere
	 c=MakeNullCleft();     // 1st cleft
	 c->start=c->end=s;
	 s->cleft=c;
   }else{                   // for each existing cleft
	 c=clefts;
	 do{                    
	   t=c->start;          // for each sphere starting from start
	   do{                  // see if new sphere should be connected to it
	 if(dist3d(t->center,center) <= t->radius + radius){
	   //printf("adding new sphere to existing cleft %d\n",c->label);
	   c->num_spheres++;
	   connect=true;               // if yes, then

	   s->next=t->next;
	   t->next->prev=s;
	   t->next=s;
	   s->prev=t;

	   s->cleft=c;                 // assign s to cleft c
	   if(c->end == t){c->end=s;}  //if there was only one before
	 }
	 t=t->next;
	   }while(t != c->end->next && !connect);
	   c=c->next;
	 }while(c != clefts && !connect);  // stop if !connect to any or all were tested

	 if(!connect){          // if !connect means all clefts were checked
	   ADD(spheres, s);
	   c=MakeNullCleft();
	   c->start=c->end=s;
	   s->cleft=c;
	 }
   }

   for(i=0;i<3;i++){s->center[i]=center[i];}
   s->radius=radius;
   //s->num_gpoints=0;

   return;
}
/***********************************************************************/
/*        1         2         3         4         5         6         7*/
/*234567890123456789012345678901234567890123456789012345678901234567890*/
/*        1         2         3         4         5         6         7*/
/***********************************************************************/
tCleft MakeNullCleft( void ){
   tCleft  c;
   
   NEW( c, tsCleft );
   ADD( clefts, c );

   if(c->next == c){
	 c->label=0;
   }else{
	 c->label=c->prev->label +1;
   }
   c->volume=0.0;
   c->num_spheres=1;

   return c;
}

/***********************************************************************/
/*        1         2         3         4         5         6         7*/
/*234567890123456789012345678901234567890123456789012345678901234567890*/
/*        1         2         3         4         5         6         7*/
/***********************************************************************/
void print_cleft_surface(char filename[], tCleft c){

  //int count=0;
  FILE *outfile_ptr;
  //char type[]={'N','C','O','S'};

  printf("filename for cleft coordinates: %s\n",filename);
  outfile_ptr = fopen(filename, "w");
  if(outfile_ptr == NULL){
    fprintf(stderr, "ERROR: unable to open input file: %s\n",filename);
    exit(8);
  }

  fprintf(outfile_ptr,"REMARK  Volume of Cleft %d: %7.3f\n",c->label,c->volume);
  //if(count < 9999){count++;}
  //fprintf(outfile_ptr,"HETATM %4d   %c CFT A %3d    %8.3f%8.3f%8.3f\n",
  //count,type[ctype],c->label,grd[l].coor[0],
  //grd[l].coor[1],grd[l].coor[2]);


  fclose(outfile_ptr);

  return;
}

/***********************************************************************/
/*        1         2         3         4         5         6         7*/
/*234567890123456789012345678901234567890123456789012345678901234567890*/
/*        1         2         3         4         5         6         7*/
/***********************************************************************/
float dist3d(float a[], float b[]){
  float dist=0.0f;
  int i;

  for(i=0;i<3;i++){
    dist += (a[i]-b[i])*(a[i]-b[i]);
  }
  dist = (float)sqrt(dist);

  return dist;
}
/***********************************************************************/
/*        1         2         3         4         5         6         7*/
/*234567890123456789012345678901234567890123456789012345678901234567890*/
/*        1         2         3         4         5         6         7*/
/***********************************************************************/
void read_commandline(int argc, char *argv[]){
  int   i,j,k;
  char  usage[3000];
  char  helper[10];
  char  chain,alt;
  char  outbase[100];
  char  anchor_nam_copy[4];

  strcpy(usage,"Usage:\n      Get_Cleft [obligatory-arguments list] ");
  strcat(usage,"[optional-arguments list]\n\n");
  strcat(usage,"Obligatory Arguments:\n");
  strcat(usage,"-p <filename> : PDB filename\n");
  strcat(usage,"Optional Arguments:\n");
  strcat(usage,"-l [1.50]        : min sphere radius.\n");
  strcat(usage,"-u [4.00]        : max sphere radius.\n");
  strcat(usage,"-h               : output hetero group atoms in cleft.\n");
  strcat(usage,"-H               : output all atoms of hetero groups in cleft.\n");
  strcat(usage,"-c chain_id      : Chain ID to be considered, all if none.\n");
  strcat(usage,"-t <value>       : Number of clefts to be modelled (0=all).\n");
  strcat(usage,"-[a|i|b] RESNUMCA: Where (RES=name NUM=number C=chain A=Alt) is a protein residue\n");
  strcat(usage,"                   in contact with the cleft or a hetero molecule within the cleft\n");
  strcat(usage,"                   to be printed out.\n");
  strcat(usage,"                   \'-a\' prints out all atoms in the selected cleft.\n");
  strcat(usage,"                   \'-i\' prints out all atoms in contact (see option \'-k\') with\n");
  strcat(usage,"                        RESNUMCA in the selected cleft. This option is appropriate\n");
  strcat(usage,"                        if RESNUMCA is a hetero group.\n");
  strcat(usage,"                   \'-b\' produces two output files corresponding to each of the above.\n");
  strcat(usage,"                   (Use \'-\' for Blank in chain and Alt loc).\n");
  strcat(usage,"-ca              : include C-alpha of residues.\n");
  strcat(usage,"-cb              : include C-beta of residues.\n");
  strcat(usage,"-r               : include all otoms of the residue.\n");
  strcat(usage,"-s               : output cleft spheres (centre coordinates and radii).\n");
  strcat(usage,"-k [5.0]         : Threshold distance for contact definition.\n");
  strcat(usage,"-o outfile       : Output filename.\n");

  // assignment of default values to optional parameters
  top_clefts=0;
  sphere_lwb=1.5;
  sphere_upb=4.0;
  chn_counter=-1;
  output_het=0;
  het_whole=0;
  anchor_flg=0;
  contact_threshold=5.0;
  output_im_cleft=0;
  output_omim_clefts=0;
  output_extra_atoms=0;
  calpha=0;
  cbeta=0;
  complete_residue=0;
  output_spheres=0;
  strcpy(outbase,".");
  // copy argv values to the respective global arguments
  
  if(argc == 1){
	  printf("%s", usage);
	  exit(0);
  }

  for(i=0; i<argc; ++i){
	  
	  if(strcmp(argv[i],"-p")==0){
		  strncpy(pdb_base,argv[i+1],4);
		  pdb_base[4]='\0';
		  
		  strcpy(pdb_file,argv[++i]);
		  if(strstr(pdb_file,".pdb") == NULL) strcat(pdb_file,".pdb");
		  
	  }else if(strcmp(argv[i],"-t")==0){
		  sscanf(argv[++i],"%d",&top_clefts);

	  }else if(strcmp(argv[i],"-l")==0){
		  sscanf(argv[++i],"%f",&sphere_lwb);

	  }else if(strcmp(argv[i],"-u")==0){
		  sscanf(argv[++i],"%f",&sphere_upb);

	  }else if(strcmp(argv[i],"-h")==0){
		  output_het=1;

	  }else if(strcmp(argv[i],"-H")==0){
		  output_het=1;
		  het_whole=1;

	  }else if(strcmp(argv[i],"-c")==0){
		  chn_counter++;
		  chn[chn_counter]=argv[++i][0];

	  }else if(strcmp(argv[i],"-a")==0 || strcmp(argv[i],"-i")==0 || strcmp(argv[i],"-b")==0){
		  // XXXNNNNAB
		  // XXX = required %3s 3 letter code (- to add space for shorter residue name)
		  // NNNN = integer, any length residue number
		  // A = chain id, '-' for blank required
		  // B = altloc id, '-' for blank required
		  for(j=0;j<3;j++) anchor_nam[j]=argv[i+1][j];
		  anchor_nam[3]='\0';
		  
          // The following block of code is used to manage spaces in residue names (spaces must be before characters !)
          int nSpaces = 0;
          for(j=0;j<3;j++)
          {
            if(anchor_nam[j]=='-')
            {
                nSpaces++;
                anchor_nam[j]=' ';
            }
          }
          switch(nSpaces)
          {
            case 0: break;
            case 1: 
                if(anchor_nam[0] != ' ' && anchor_nam[2] == ' ')
                {
                    anchor_nam[2] = anchor_nam[1];
                    anchor_nam[1] = anchor_nam[0];
                    anchor_nam[0] = ' ';
                }
                break;
            case 2:
                j=0;
                while(anchor_nam[j] == ' ') j++;
                anchor_nam[2] = anchor_nam[j];
                anchor_nam[1] = ' ';
                anchor_nam[0] = ' ';
                break;
          }
          
		  
		  for(j=0;j<(int)(strlen(argv[i+1])-1);j++) helper[j]=argv[i+1][3+j];
		  helper[j]='\0';
		  sscanf(helper,"%d",&anchor_num);
		  anchor_chn=argv[i+1][strlen(argv[i+1])-2];
		  
		  if(anchor_chn == '-') anchor_chn=' ';
		  anchor_flg=1;
		  anchor_alt=argv[i+1][strlen(argv[i+1])-1];
		  if(anchor_alt == '-') anchor_alt=' ';
		  
		  //printf("name='%s'\tnum='%d'\tchn='%c'\talt='%c'\n",anchor_nam,anchor_num,anchor_chn,anchor_alt);
		  if(strcmp(argv[i],"-i")==0){
			  output_im_cleft=1;
		  }else if(strcmp(argv[i],"-b")==0){
			  output_omim_clefts=1;
		  }
		  
		  i++;
		  //printf("ANCHOR: nam=%s num=%d chn=<%c> alt=<%c>\n",anchor_nam, anchor_num,
		  //     anchor_chn,anchor_alt);
		  //PAUSE;
	  }else if(strcmp(argv[i],"-k")==0){
		  sscanf(argv[++i],"%f",&contact_threshold);
	  }else if(strcmp(argv[i],"-r")==0){
		  output_extra_atoms=1;
		  complete_residue=1;
	  }else if(strcmp(argv[i],"-ca")==0){
		  output_extra_atoms=1;
		  calpha=1;
	  }else if(strcmp(argv[i],"-cb")==0){
		  output_extra_atoms=1;
		  cbeta=1;
	  }else if(strcmp(argv[i],"-s")==0){
		  output_spheres=1;
	  }else if(strcmp(argv[i],"-o")==0){
		  strcpy(outbase,argv[++i]);
	  }
    
  }
  
  if(strcmp(outbase,"") == 0){
    strcpy(outbase,pdb_base);
  }

  if(anchor_flg == 1){
    top_clefts=0;

    if(anchor_chn == ' '){
      chain='-';
    }
    else {
      chain=anchor_chn;
    }
      
    if(anchor_alt == ' '){
      alt='-';
    }
    else {
      alt=anchor_alt;
    }
    for(k=0;k<4;k++)
    {
        if (anchor_nam[k] == ' ') anchor_nam_copy[k] = '-';
        else anchor_nam_copy[k] = anchor_nam[k];
    }
    sprintf(out_base,"%s_%s%d%c%c",outbase,
	    anchor_nam_copy,anchor_num,chain,alt);
  }else{
    sprintf(out_base,"%s",outbase);
    if(chn_counter > -1){
      //printf("chn_counter=%d\n",chn_counter);
      strcat(out_base,"_");
      for(i=0;i<=chn_counter;i++){
	//printf("%d <%c>\n",i,chn[i]);
	sprintf(out_base,"%s%c",out_base,chn[i]);
      }
    }
  }

  //printf("out_base=<%s>",out_base);
  //PAUSE;

  return;
}
/***********************************************************************/
/*        1         2         3         4         5         6          */
/*234567890123456789012345678901234567890123456789012345678901234567890*/
/*        1         2         3         4         5         6         7*/
/***********************************************************************/
void read_pdb(char file[]){
  FILE  *infile_ptr;        // pointer to input file
  char  buffer[81];         // a line from the INPUT file
  char  coor_char[10];      // string used to read the coordinates
  char  field[7];
  int   i,j,k;
  tAtom a,b;
  tRes  r;
  int  rnum;
  char rnam[4];
  char chain;
  char altloc;

  //for(i=0;i<3;i++){
  //pdb_min[i]=1000;
  //pdb_max[i]=-1000;
  //}
  int chn_flag;
  int chn_c;
  int read_line;
  int prot_type;
  //int count;

  /*
    char    line[82];   // PDB line
    char    name[5];    // atom name
    float   radius;     // atomic radius
    float   coor[3];    // coordinates
    tSphere ofSphere;  // pointer to sphere to which atom belongs
    int     isprot;     // 1 if ATOM - 0 if HETATM
    int     inum;
    tAtom   next,prev;  // connects the list of all atoms
    tAtom   next_inres; // link to next in list of all atoms in the same residue
    tAtom   prev_inres; // link to prev in list of all atoms in the same residue
    tRes    ofres;
  */


  infile_ptr = fopen(file, "r");

  if(infile_ptr == NULL){
    fprintf(stderr, "ERROR: unable to open input file: %s\n",file);
    exit(8);
  }
  while (fgets(buffer, sizeof(buffer),infile_ptr)){
    read_line=0;
    prot_type=0;
    for(i=0;i<6;i++) field[i]=buffer[i];
    field[6]='\0';

    if(strcmp(field,"ATOM  ") == 0 || strcmp(field,"HETATM") == 0){

      for(i=0;i<3;i++) rnam[i]=buffer[i+17];
      rnam[3]='\0';
      for(i=0;i<4;i++){coor_char[i]=buffer[22+i];}
      coor_char[4]='\0';

      sscanf(coor_char,"%d",&rnum);
      chain=buffer[21];
      altloc=buffer[16];

      read_line=1;
      if(strcmp(field,"ATOM  ") == 0){
	// check if it is the correct chain
	prot_type=1;
	if(chn_counter > -1){
	  chn_flag=0;
	  chn_c=0;
	  while(chn_c  <= chn_counter && chn_flag == 0){
	    if(buffer[21] == chn[chn_c]) chn_flag=1;
	    chn_c++;
	  }
	  if(chn_flag == 0) read_line=0;
	}
	
	// check if it is the correct altloc, unless it is the anchor
	// either ' ' or 'A' is accepted, otherwise, same anchor_alt
	if(strcmp(rnam,anchor_nam)==0 && 
	   rnum == anchor_num && chain == anchor_chn){
	  if(altloc != anchor_alt) read_line=0;
	}else if(altloc != ' ' && altloc != 'A'){
	  read_line=0;
	}
      }
      else if(strcmp(field,"HETATM") == 0){
	// check if it is the correct altloc, unless it is the anchor
	// either ' ' or 'A' is accepted, otherwise, same anchor_alt
	if(strcmp(rnam,anchor_nam)==0 && 
	   rnum == anchor_num && chain == anchor_chn){
	  if(altloc != anchor_alt) read_line=0;
	}else if(altloc != ' ' && altloc != 'A'){
	  read_line=0;
	}
      }
    }
    if(read_line == 1){
      NEW( a, tsAtom );
      ADD( atoms, a );      
      if(a->next == a){
	a->inum=0;
      }else{
	a->inum=a->prev->inum +1;
      }
      //a = MakeNullAtom();

      for(i=0;i<=3;i++){a->name[i]=buffer[i+12];}
      a->name[4]='\0';

      a->isprot=prot_type;

      a->radius=assign_radius(a->name);

      for(i=0;i<81;i++) a->line[i]=buffer[i];
      a->line[81]='\0';

      a->out=0;
      //printf("%s\n",a->line);
      //PAUSE;
      //printf("ATOM <%s> %7.3f\n",a->name,a->radius);
      //PAUSE;

      for (j=0;j<=2;j++){
	for(i=0;i<=7;i++){coor_char[i]=buffer[30+j*8+i];}
	coor_char[8]='\0';
	//printf("j=%d coor_char=<%s>\n",j,coor_char);
	sscanf(coor_char,"%f",&a->coor[j]);
	//if(a->coor[j] < pdb_min[j]){pdb_min[j]=a->coor[j];}
	//if(a->coor[j] > pdb_max[j]){pdb_max[j]=a->coor[j];}
      }
      r=assign_residue(rnam,rnum,chain,altloc,a);
      a->ofres=r;

      a->ofres->out=0;
      //printf("r->(name=%s num=%d chain=%c): a->name=%s\n",a->ofres->name,
      //   a->ofres->num,a->ofres->chain,a->ofres->atom_anchor->name);
      b=r->atom_anchor;
      //count=0;
      //do{
	//count++;
	//printf("b->(name=%s num=%d chain=%c): a->name=%s\n",b->ofres->name,
	//    b->ofres->num,b->ofres->chain,b->name);
      //PAUSE;
      //b=b->next_inres;
      //}while(b != r->atom_anchor);
      //printf("count=%d\n",count);
      //PAUSE;
    }
  }
  fclose(infile_ptr);    

  return;
}
/***********************************************************************/
/*        1         2         3         4         5         6         7*/
/*234567890123456789012345678901234567890123456789012345678901234567890*/
/*        1         2         3         4         5         6         7*/
/***********************************************************************/
tRes assign_residue(char rnam[],int rnum, char chn, char alt, tAtom a){
  tRes res = NULL,new_ = NULL;
  int make_new_res=0;
  int old=0;
  /*
    LIST OF MEMBERS OF STRUCTURE RESIDUE:
    char  name[4];    // residue name
    int   num;        // residue number
    char  chain;      // residue chain
    tAtom atom_anchor;     // first atom read found to be part of this residue
    int   hetero;     // 0 for aminoacids, 1 for hetero groups
    tRes  next,prev;  // link to other residues in the list
  */

  /*
    printf("-------------------------------\n");
    printf("rnam=%s rnum=%d chn=%c:\n",rnam,rnum,chn);
    printf("line=%s\n",a->line);
    printf("name=%s\n",a->name);
    printf("num =%d\n\n\n",a->inum);
    PAUSE;
  */

  if(resids == NULL){
    make_new_res=1;
  }else{
    res=resids;
    do{
      if(strcmp(res->name,rnam) == 0 && res->num == rnum && res->chain == chn){
	old = 1;
	new_ = res;
      }
      res=res->next;
    }while(old==0 && res != resids);

    if(old == 0) make_new_res=1;
  }

  if(make_new_res==1){
    NEW( new_, tsRes );
    ADD( resids, new_ );
    strcpy(new_->name,rnam);
    new_->chain=chn;
    new_->alt=alt;
    new_->num=rnum;
    new_->atom_anchor=a;
    new_->cleft_anchor=0;
    new_->out=0;
    if(a->isprot == 1){
      new_->hetero=0;
    }else{
      new_->hetero=1;
    }
    a->next_inres=a->prev_inres=a;
  }else{
    a->next_inres=new_->atom_anchor;
    new_->atom_anchor->prev_inres->next_inres=a;
    a->prev_inres=new_->atom_anchor->prev_inres;
    new_->atom_anchor->prev_inres=a;
  }

  /*
    res=resids;
    do{      
    printf("res->(name=%s num=%d chain=%c): a->name=%s\n",res->name,
    res->num,res->chain,res->atom_anchor->name);
    res=res->next;
    }while(res != resids);
    PAUSE;
  */

  return new_;
}
/***********************************************************************/
/*        1         2         3         4         5         6         7*/
/*234567890123456789012345678901234567890123456789012345678901234567890*/
/*        1         2         3         4         5         6         7*/
/***********************************************************************/
float assign_radius(char atm[]){
  float r;
  char  c;

  c=atm[1];
  
  switch(c){
  case 'N':
    r=1.70f;
    break;
  case 'O':
    r=1.50f;
    break;
  case 'C':
    r=1.90f;
    break;
  case 'S':
    r=1.90f;
    break;
  case 'P':
    r=1.23f;
    break;
  case 'L':
    r=1.75f;
    break;
  case 'F':
    r=1.70f;
    break;
  case 'R':
    r=1.85f;
    break;
  case 'I':
    r=2.00f;
    break;
  case 'E':
    r=1.50f;
    break;
  case 'A':
    r=1.50f;
    break;
  default :
    r=1.50f;
    break;
  }

  return (r);
}
