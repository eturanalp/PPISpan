// goinpin2.cpp : Defines the entry point for the console application.
//
/*
THE ALGORITHM FOR PINSPAN:
ARGS:	pinf: tab delimited data file that contains PIN info [p1,p2,confidence] e.g. LeeConfidentNet.net.txt
		pf: protein file, one protein per line with label e.g. LeePRoteinsWithAspectFandP_tabbed.txt
		goasf: GeneOntology Associations file in unfiltered form
		goobof: GeneOntology terms and OBO file(gene_ontology_edit.obo.txt)
				(for obo format see www.geneontology.org)
        argv[4]= test file containing a subG , the contents of this file is read into a list of subGs.
		         If "rdp:filename" is used then the duplicate patterns in the given file
		         are removed and the result is written to rdp_filename. 
				 If "rdp:" is used then the maxresult_file is cleared of duplicates
				 If "bonfer:" is used then the Z-scores in rdp_maxresultfile are corrected
		argv[5]= support  #15
		argv[6]= maxdepth  #8
		argv[7]= tlevel  (target level in GTT) or GO subset name as in obo file (e.g. goslim_yeast) #3
		argv[8]= exclude or not  (Certain edges were excluded:(3735-3735)(3723-3723)(7046-7046)) # not
		argv[9]= use_embeddings or not  #not
		argv[10]=the file of list of proteins to be removed prior to discovery (one protein per line) #not
		         or not
        argv[11]= Which Statistical significance score to use for ordering of results: "P" or "Z"  #Z
		          (default is Z)
        argv[12]= Random Graph Model (RGM=F for frequency preservation, RGM=L for locality or RGM=NONE) #RGM=F
		srgv[13]= FIARN(find in a random network) or not: Finds the patterns in one of random networks 
PRE: Each protein name and label is less than 10 chars long. Input files are all tabbed text
        pf and pinf contains consistent and complete data 

OUTPUT:
   1- maxresult_file: The file that contains maximal frequent subgraphs
   2- resfile: This file contains all frequent subgraphs along with their embeddings
   3- rankfile: This file contains all frequent subgraphs in resfile sorted by Z-scores
   4- positions: All embeddings appended
   5- rdp_file: Contents of maxresult_file where duplicates are removed

SAMPLE INVOCATION:
goinpin7_sp2.exe human_ophid_go_annotations.txt human_ophid_PIN_annotated.txt gene_ontology_edit.obo.txt
                 subG3.txt 100 8 0 exclude use_embeddings toberemoved.txt Z RGM=F
goinpin7_sp2.exe GOannotations_for_SGD_thru_EXPASY_swissprot.txt SGD_PIN_from_DIP_swissprot.txt 
                 gene_ontology_edit.obo.txt subG3.txt 50 6 goslim_yeast not not not Z NONE
goinpin7_sp2.exe GOAnnotations_yeast_tolgadan.txt PIN_DIP_04112007_tolgadan.txt 
                 gene_ontology_edit.obo rdp: 30 8 goslim_yeast not use_embeddings not Z RGM=F

STEPS :
1. Read each file and create an eficient data structure suitable for graph operations
2. Read into GoTermTable, the obo GO terms file, create the tree structure (GTT)
3. Set the label of each node to its parent in the levelth level in GTT
4. Calculate label pair frequencies and sort.
5. Run frequent pattern discovery algorithm:
    For each label pair in increasing freq order
	  generate all subgraphs that can be produced from it using gSpan like enumeration technique
	      Run our novel subgraph isomorphism test algorithm using DFS tree of the subGraph

NOTES:
** If you want to find all embeddings (some of which may be overlapping) call mapped_or_false_connections()
   function with parameter allow_overlap==true
** "use_embeddings" ve "not" pgm parametrelerinin maxsubs dosya çiktilari arasinda olusan fark(bazi freq.
   subgraphlar "not" çiktisinda varken use_embeddings çiktisinda yoklar) maximum cardinality set'i bulmaya
   çalismamamizdan kaynaklanmaktadir. Embeddingleri deneme ve bulus siramiz sonucu etkilemektedir.
** use_embeddings oldugunda, child patternlarin embeddinglerini parentlerinin embedding listelerini 
   kisaltarak buluyoruz.
 
UPDATE HISTORY:
DATE		UPDATE
12/22/06    Overlaping vs. non overlaping embeddings
12/30/06	Added subG.embeddings linked list to facilitate the use of exact locations of
			embeddings of s in the search for a child of s. 
			Added find_subG2() , updated class subG , mapG2() 
01/07/07    Added argv[8] and argv[9]
01/12/07    Added argv[10] for the proteins to ve removed from pht prior to processing
03/29/07    Added create_randomPINS(pht,10,20000) , modified read_subG_find(argv[4],pht)
03/31/07    Added random_graph_statistics() , compute_Z_score(s) and find_subG3(temp,pin)
04/10/07    Added listing of results according to the Zscore, 
            append_to_file() now runs calculate_zscore() ,
			modified terrible hash function= 10x improvement with new line:r=r*(x<<(x%8));
04/16/07    Support for multiple-labelled proteins is added
11/09/07    Support for GO subsets is implemented. GO subsets such as goslim_yeast or 
            goslim_generic are cut-down subsets of terms in GO. They give a broad overview 
			of the ontology content without the detail of the specific fine grained terms. 
			Now, a subset name such as goslim_yeast may be specified instead of target GO level.
11/26/07    subG.is_subgraph() is added , remove_duplicate_patterns() added
02/19/08    find_embeddings_of_non_labelled_patterns() is added
02/21/08    bonferroni_correct_Z_score() is added, bonfer: parameter (agv[4]) is added,
            now rdp: runs with Z-scores bonf. corrected


TODO LIST:
Priority/
Severity(1-5)            ToDo
-------------  -------------------------------------------------
-- Memory deallocation after min_s(), correct permutation in min_s() 
-- maxsubs dosyasinda maximal olmayan örüntüler neden var?
-- Output embeddings info for patterns in maxsubs file, not just in results file
-- Farkli min_support parametreleri için neden ayni örüntünün supportu farkli bulunuyor?
   Bunun için find_subG2() incelenmeli..
-- find_subG2() bütün olasiliklar için subG2()'yi çagiriyor mu? Compare it with subG.is_subgraph()

*/

#include <iostream>
//#include <stdafx.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include "pValue.h"
//#include <limits>

unsigned int hts;
unsigned int embed[11];
int support=10;
int maxdepth=6; // maximum depth of DFS CODE TREE to explore
int pinn; // protein interaction number
int tlevel=5;  // target level
char maxresult_file[100];  // the output file that  contains the maximal freq subgraps
char resfile[100];  // result file that contains all freq patterns
char rankfile[100];   // rank sorted result file
int use_embeddings=0;  //boolean
int exclude=0;    // boolean parameter , see exclude()
char go_subset[40]; //[40]="goslim_yeast";  // the target GO subset such as "goslim_yeast" instead of a target GO level
long s_id=1;  // static id of subG

struct protein{
	char name[10];
	float conf;  // confidence
	protein * next; };
	
struct emb{
    unsigned int embed[11];
	struct emb* next;
	};

class subG{
   public:
   int labels[10];   // maximum 10 node'lu graph buluyoruz 
   int nol; // number of labels
   short int adjm[10][10];  // adjacency matrix
   subG* next;
   struct emb* embeddings;  // link to the root of linked list of found embeddings of this subG
   int grow;   // the index of extending node from which we grow a an edge
   int status;
   long id;
   subG(){
      nol=0; next=NULL; embeddings=NULL; grow=0; Zscore=0;freq=0; pValue=0; status=0; id=s_id; s_id++; 
      }
   ~subG(){  // delete the embedings linked list
       struct emb* t;
	   for (struct emb* temp=embeddings; (temp) ; temp=t ) {t=temp->next; delete temp;}
      }
   float Zscore;
   int freq;   // frequency of subG in pht // number of embeddings?
   float pValue;
   int nedges(){int c=0; for (int i=0; (i<nol) ;i++) for (int j=i+1; (j<nol) ;j++) if (adjm[i][j]==1) c++; return c;}

   bool is_subgraph(subG* s){  // is s a subgraph of this
	   unsigned int mapx[11];
	   for(int i=0 ; i<nol ; i++){
		  for (int j=0; j<s->nol ; j++) {
			  if (labels[i]==s->labels[j]) {
                 for (int k=0; k<11 ; k++) mapx[k]=0; 
     	         if (mapGr(j,mapx,i,s))  {
					//printf("\nOne subgraph map: ");
  		            //for (int t=0; t<s->nol ; t++)  printf("%d ",embed[t]);            
					return true;
			        }
		         }
	          }
	      }
       return false;
      }

   int mapped_or_false_connections(int x, int t,  unsigned int map[11],subG* s){
    // is protein_no already mapped or is mapping it with x will induce incorrect edges?
    // it is known that they have the same label
    // overlapping between embeddings may be allowed using last_subG field of PROINFO in pht
	for (int i=0; ((i< s->nol)) ; i++ )  {
		if (t==(map[i]-1)) return 1;  //already mapped
	   }
	// if for each neighbour j of x in s  
	for (int j=0; ((j< s->nol)) ; j++ )  
		 if ((map[j] > 0)&&(s->adjm[x][j] == 1))  {
              // see if map[j] is neighbour of map[x] and if the label of map[j]=label of j
			 if (adjm[map[j]-1][t]!=1)  // false connection: not neighbours
				  return 1;   
			 else if (s->labels[j]!=labels[map[j]-1]) return 1; // false label
		     }
	return 0;
    }

   unsigned int mapGr( int v, unsigned int map[11], int vthis, subG* s){
		// v: the index of new node in subGraph s that is going to be mapped to vthis'th node of this
		// map[] is the set of visited nodes(node indexes) of this ordered by their respective indexes in s->labels
	    //       ;for example, map[2]==5 means the 2nd node of s is mapped to 5th node of this
		// the result set of proteins is copied to embed[]
	    // returns 0 if no match was found
	    // returns 1 otherwise

       unsigned int *newmap;
	   unsigned int p;
	   unsigned int *maplist[200]; 
	   unsigned int *newmaplist[200]; 
       int nom=0;  // the number of maps in maplist
	   int newnom=0;
	   int c=0;
	   bool first=true; bool insidefirst=false; bool nextneighbour=false;
	   map[v]=vthis+1;
	   map[10]++;
       if (map[10]== s->nol)  { 
		   memcpy(embed,map,11*sizeof(int)); 
		   //printf("\nFound embedding"); print_protein_array(embed); 
		   return 1;
	       }
	   for (int i=0; i< s->nol ; i++ )  {  //for each neighbour i of v (who is not yet mapped)
		   if ( (s->adjm[v][i]==1) &&(map[i]==0)&&(v!=i)) {  
			  for (int t=0; ((t<nol)&&((nom+newnom)<=150)) ; t++ )  {
				   // for each neighbour t of vthis (who is not yet mapped) whose label is equal to s->labels[i]           
				   if ((labels[t])&&(labels[t]==s->labels[i])){		// if i and t have the same label
					   if (first) {
					        if (!mapped_or_false_connections(i,t,map,s)) { 
		                       newmap=(unsigned int*)malloc(11*sizeof(unsigned int));
			                   newmap=(unsigned int*)memcpy(newmap, map, 11*sizeof(unsigned int));    
                               maplist[nom]=newmap; nom++;  // add new mapping to the list
							   if (nom>25) { //printf("\nToo many candidates for a single label %dth=%s..",i,s->labels[i] );
								             break;}
						       if (mapGr(i,newmap,t,s)) {for(int k=0;k<nom;k++) free(maplist[k]); 
							                                     for(int k=0;k<newnom;k++) free(newmaplist[k]);
							                                     return 1;		}
							   insidefirst=true;
						       }
					       }
					   else {
							   for(int k=0;k<nom;k++){
								   if (maplist[k][i]!=0) {maplist[k][i]=0;maplist[k][10]--;}     
								   if (!mapped_or_false_connections(i,t,maplist[k],s)) {
		                                newmap=(unsigned int*)malloc(11*sizeof(unsigned int));
			                            newmap=(unsigned int*)memcpy(newmap, maplist[k], 11*sizeof(unsigned int));    
                                        newmaplist[newnom]=newmap; newnom++;  // add new mapping to the list
										if (newnom>25) { //printf("\nToo many candidates for a single label %dth=%s..",i,s->labels[i] );
											             break;}
						                if (p=mapGr(i,newmap,t,s)) {for(int k=0;k<nom;k++) free(maplist[k]); 
										                                  for(int k=0;k<newnom;k++) free(newmaplist[k]);
																		  return 1;}             
								        nextneighbour=true;
								       }
								    }  // for

				                }
			            }
			      }
			  if (insidefirst) 
				   first=false;
			  if (nextneighbour) {
			      nextneighbour=false;
			      for(int k=0;k<nom;k++) if (maplist[k]) free(maplist[k]);  
				  nom=newnom;
			      newnom=0;
			      for(int k=0;k<nom;k++) { 
				     maplist[k]=newmaplist[k];
			         }
			      }
		      } 
	      }
	   if (first) return 0;
	   else if (nom>0) {
           for(int k=0;k<nom;k++) {
             if (map[10]<maplist[k][10])  memcpy(map, maplist[k], 11*sizeof(unsigned int));  
		     }  
	       }
	   for(int k=0;k<nom;k++) if (maplist[k]) free(maplist[k]);
	   return 0;
	  }

   };

class PROINFO {  // protein info class , one object per protein 
public:
	char *id;
	int label;
 	int otherlabels[5];
	int non;  // number of neighbours
	struct protein* nl;  // neighbour list
	long last_subG;   // points to the last freq subG that this node was mapped to
	                   // used to detect overlapings
	PROINFO(){
	  id=NULL;
	  label=0;
	  non=0;
	  nl=NULL;
	  last_subG=0;
	  for(int i=0; i<5 ; i++) otherlabels[i]=0;
 	  }

	// copy constructor
	PROINFO(const PROINFO &p){
       *this=p;
	   }
   
	int add_other_label(int ilabel){
		if (ilabel==label) return 0;
		for(int i=0; i<5 ; i++)  if (otherlabels[i]==0) {otherlabels[i]=ilabel; return ilabel;}
		                         else if (ilabel==otherlabels[i]) return 0;
		return 0; 
	}

	int is_my_label(int ilabel){  // is ilabel one of my other labels // what about the primary label?
        for(int i=0; i<5 ; i++)  if (otherlabels[i]==ilabel) {return ilabel;}
		return 0;
	}

};

class goterm {
public:
	unsigned int id;
	char *name;
	char namesp;  // F for molecular_function, P for biological_process and C for Cellular_compartment 
	char *def;
	unsigned int is_a[11];
	unsigned int part_of[11];
	int level_plus1;
	char *subset;
	goterm(){
	  name=NULL;def=NULL;subset=NULL;namesp='\0';
	  id=0;
	  level_plus1=0;
	  for (int k=0;k<11;k++) is_a[k]=part_of[k]=0;
	  }
};

typedef struct labelpair{
   int label1;
   int label2;
   int freq;
} labp;

class subGlist{  // list of subG's ordered by Z-score
public:
subG* list;
subG* listend;  // points to subG that points to the last element
int maxlength;
int length;

subGlist(int max){
   list=listend=NULL;
   maxlength=max;  // this is of no use
   length=0;
   }
~subGlist(){
  subG* t;
  for (subG* temp=list; (temp) ; temp=t ) {t=temp->next; delete temp;}
}
subG* add_to_simple_list(subG* ss){
   subG* temp;
   subG* oldtemp;
   subG* s=new subG(*ss);  // bu copy constructor çalisinca s->id ne olur? Bil bakalim.
   if (!s) return NULL;
   for (temp=list,oldtemp=NULL; ((temp)&&(temp->Zscore < s->Zscore)) ; temp=temp->next) oldtemp=temp;
   s->next=temp;
   if (oldtemp) oldtemp->next=s;
   else list=s;
   if (!temp) listend=s; // NOT TESTED: we need to update listend as well
   return s;
}
};  // class subGlist


typedef struct DFScodeEdge{
	int v1;  // the discovery order number of the first vertex of the edge
 	int v2;  // the discovery order number of the second vertex of the edge
	int ei;  // the edge index in soalp;
} DFScode_edge;

class DFScode{
public:
	DFScode_edge cedges[50];
	int noedges;  // number of edges
	DFScode(){
      noedges=0;
	  //for (int i=0; i<50 ; i++) {cedges[i].ei =0;cedges[i].v1=0; cedges[i].v2=0;}
	  }
};

labp* lpt;  // label pair table  
goterm* gtt[100000];  // go term table
PROINFO* pht;  // protein hash table
int* soalp;    // sorted array of label pairs (holds indexes to lpt)(in asc order)
int nolp;  // number of label pairs
int ei;    // edge index in soalp that points to the edge that is being processed
int cc=0;  // count of mapG2 calls
unsigned int emptyarray[11];
PROINFO * randomPIN[1000];
int randomPINcount=0;
subGlist* toplist; // of highest Z-score patterns
char command[300]="";
pValue pV;
const int ZSCORE=1;
const int PVALUE=2;
int ordering_value=ZSCORE;
char* RGM;  // Random Graph Model (RGM=F for frequency preservation, RGM=L for locality or RGM=NONE)
int pruned=0;

void display_subG(subG* s,unsigned int p[11], FILE* ff);
void print_proteins(unsigned int p[11]);
void print_protein_array(unsigned int p[11]);
void merge_maps(unsigned int**maplist, int nom, unsigned int *map, int v, int nol);
int verify_sub(unsigned int map[11],subG* s, bool show);
void read_subG(subG* s,unsigned int p[11],char *sf);
bool neighbours(unsigned int p1, unsigned int p2);
int append_to_file(char* rfname ,subG* s ,unsigned int pr[11], unsigned int positions[1000][11]);
void statistics();
subG* random_subgraph_generator(int size, int* v,unsigned int* p,unsigned int proteins[11],int* nnodes, int* nedges);
PROINFO* copyPIN(PROINFO* pht);
bool randomizePIN(PROINFO* pin, int itn);
int create_randomPINS(PROINFO* pht, int count, int rand_count);
float compute_Z_score(subG* s);
long double bonferroni_correct_Z_score(double Zscore, int experiments);

unsigned int hashf(char* s, unsigned int hts) {  // s is 10 char string
  //unsigned char *i=(unsigned char *)s;
  //int len=strlen(s);
  unsigned int r=1; int x;
  //printf("%c _ ",(char )(*s));
  for (char *i=(char *)s; (((*i)!='\0')); i++) {
	  x=(*i);
	  r=r+x*x; 
	  r=r*(x<<(x%8)) ;  // shift left at most 7 times
	  //r=r+((*i)-40)*((*i)-40); n++; 
	  //printf("%d_",x);
     }
  return (r % (hts-1))+1;
}
unsigned int  add_protein(PROINFO* p,char* p1, char* label, unsigned int hts){
	//inserts protein p1 into hash table p with its label . hts is hash table size , p1  is max 10 chars
	//  if p1 is already in p then label is inserted into other labels of p1.
  unsigned int i; 
  int ilabel=atoi(label);
  for(i=hashf(p1,hts); ((p[i].id!=NULL)&&(strcmp(p1,p[i].id))) ; i = (i==(hts-1)) ? 1 : i+1) ;
  if (p[i].id)  // if p1 is already in p at position i then
	p[i].add_other_label(ilabel); 
  else {   // else if p1 is new to p
    char *temp=new char[10];
    p[i].id=strcpy(temp,p1);
    //temp=new char[10];
    p[i].label=ilabel;
    }
  return i;
}

unsigned int get_protein(PROINFO* p,char* p1,  unsigned int hts){
 unsigned int i;
 for(i=hashf(p1,hts); ((p[i].id!=NULL)&&(strcmp(p[i].id,p1))) ; i = (i==(hts-1)) ? 1 : i+1) ;
 if (p[i].id==NULL) return hts;  // bu durum pinf içinde olan bazi proteinlerin pf içinde olmamasindan kaynaklanir
                                                // o nedenle verilerin tutarliligi saglandiktan sonra bu programam yüklenmelidir.
 else if (i==0) printf ("\nHash index value=0000000000");
 else return i;
 return hts;
}

int insert_assoc(PROINFO* p,char *p1, char *p2, float conf){
  unsigned int x1=get_protein(p,p1, hts);
  unsigned int x2=get_protein(p,p2, hts);
  struct protein* temp;
  if ((x1==hts)||(x2==hts)) return 0;
  if (x1==x2) return 0;  // eliminate self loops
  // else insert protein p1
  temp=new struct protein;
  strcpy(temp->name,p2);
  temp->conf=conf;
  struct protein *temp2=p[x1].nl;
  p[x1].nl=temp;
  p[x1].nl->next=temp2;
  p[x1].non ++;
  //if (x2==hts) return 0;
  // else insert protein p1
  temp=new struct protein;
  strcpy(temp->name,p1);
  temp->conf=conf;
  temp2=p[x2].nl;
  p[x2].nl=temp;
  p[x2].nl->next=temp2;       
  p[x2].non ++;
  return 1;
}

int insert_assoc(PROINFO* p,unsigned int x1, unsigned int x2, float conf){
  struct protein* temp;
  temp=new struct protein;
  strcpy(temp->name,p[x2].id );
  temp->conf=conf;
  struct protein *temp2=p[x1].nl;
  p[x1].nl=temp;
  p[x1].nl->next=temp2;
  p[x1].non ++;
  temp=new struct protein;
  strcpy(temp->name,p[x1].id);
  temp->conf=conf;
  temp2=p[x2].nl;
  p[x2].nl=temp;
  p[x2].nl->next=temp2;       
  p[x2].non ++;
  return 1;
}

int mapped_or_false_connections(int x, unsigned int protein_no,  unsigned int map[11],subG* s, bool allow_overlap, PROINFO* pht){
 // is protein_no already mapped or is mapping it with x will induce incorrect edges?
 // it is known that they have the same label
 // overlapping between embeddings may be allowed using last_subG field of PROINFO in pht

	for (int i=0; ((i< s->nol)) ; i++ )  {
		//if (map[i]==0) continue;
		if (map[i]==protein_no) return 1;  //already mapped
	   }

	if (!allow_overlap) {
       if (pht[protein_no].last_subG==s->id) return 1;
	   }
	
	// if for each neighbour j of x in s  
	for (int j=0; ((j< s->nol)) ; j++ )  
		 if ((map[j] > 0)&&(s->adjm[x][j] == 1))  {
              // see if map[j] is neighbour of map[x]=protein_no and if the label of map[j]=label of j
			 if (!neighbours(map[j],protein_no))  // false connection: not neighbours
				  return 1;   
			  else if ((s->labels[j]!=pht[map[j]].label)&&(!pht[map[j]].is_my_label(s->labels[j]))) return 1; // false label
		      }

	 return 0;
  }


	unsigned int mapG2( int v, unsigned int map[11], unsigned int protein_no, subG* s, int caller){
		// v: the index of new node in subGraph s that is going to be mapped to protein_no
		// map[] is the set of visited nodes(protein hash table indexes) ordered by their respective indexes in s->labels
		// the result set of proteins is copied to embed[]

       unsigned int *newmap;
	   unsigned int t,p;
	   unsigned int *maplist[200]; 
	   unsigned int *newmaplist[200]; 
       int nom=0;  // the number of maps in maplist
	   int newnom=0;
	   int c=0;
	   bool first=true; bool insidefirst=false; bool nextneighbour=false;
	   map[v]=protein_no;
	   map[10]++;
	   //print_protein_array(map); printf("v=%d c=%d ",v,caller);
       if (map[10]== s->nol)  { 
		   memcpy(embed,map,11*sizeof(int)); 
		   //printf("\nFound embedding"); print_protein_array(embed); 
	       for(unsigned int k=0;k<embed[10];k++)  pht[embed[k]].last_subG = s->id ;  //record embedding to prevent overlapings
		   //struct emb* tempem=new struct emb; 
		   //memcpy(tempem->embed ,map,11*sizeof(int)); 
		   //tempem->next=s->embeddings; 
		   //s->embeddings=tempem;
		   return protein_no;
	       }
	   for (int i=0; i< s->nol ; i++ )  {  //for each neighbour i of v (who is not yet mapped)
		   if ( (s->adjm[v][i]==1) &&(map[i]==0)&&(v!=i)) {  
			  for (struct protein* j=pht[protein_no].nl; ((j)&&((nom+newnom)<=150)) ; j=j->next )  {
				   // for each neighbour t of protein (who is not yet mapped) whose label is equal to s->labels[i]           
				   //printf("neighbour name  and labels: %s and %s - %s \n",j->name,pht[get_protein(pht,j->name,hts)].label,s->labels[i]);
				   if ((pht[get_protein(pht,j->name,hts)].label)&&((pht[get_protein(pht,j->name,hts)].label==s->labels[i])||
					                                               (pht[get_protein(pht,j->name,hts)].is_my_label(s->labels[i])))){
					   // if i and j have the same label
					   if (first) {
					        if (!mapped_or_false_connections(i,t=get_protein(pht,j->name,hts),map,s,false,pht)) { 
		                       newmap=(unsigned int*)malloc(11*sizeof(unsigned int));
			                   newmap=(unsigned int*)memcpy(newmap, map, 11*sizeof(unsigned int));    
                               maplist[nom]=newmap; nom++;  // add new mapping to the list
							   if (nom>25) { //printf("\nToo many candidates for a single label %dth=%s..",i,s->labels[i] );
								             break;}
						       //printf("\n [-2-]Added to maplist in first nom= %d protein %d=%d",nom, i,t);
						       if (p=mapG2(i,newmap,t,s,v) > 0) {for(int k=0;k<nom;k++) free(maplist[k]); 
							                                     for(int k=0;k<newnom;k++) free(newmaplist[k]);
							                                     return p;		}
							   insidefirst=true;
						       }
					       }
					   else {
							   for(int k=0;k<nom;k++){
								   if (maplist[k][i]!=0) {maplist[k][i]=0;maplist[k][10]--;}     
								   if (!mapped_or_false_connections(i,t=get_protein(pht,j->name,hts),maplist[k],s,false,pht)) {
								   	    //printf("\nhere comes the fork at protein %d=%d , nom=%d , c=%d, v=%d",i,t,nom,c,v);
		                                newmap=(unsigned int*)malloc(11*sizeof(unsigned int));
			                            newmap=(unsigned int*)memcpy(newmap, maplist[k], 11*sizeof(unsigned int));    
                                        newmaplist[newnom]=newmap; newnom++;  // add new mapping to the list
										if (newnom>25) { //printf("\nToo many candidates for a single label %dth=%s..",i,s->labels[i] );
											             break;}
						                //printf("\n [-3-] Added to maplist %d" , nom);
						                if (p=mapG2(i,newmap,t,s,v) > 0) {for(int k=0;k<nom;k++) free(maplist[k]); 
										                                  for(int k=0;k<newnom;k++) free(newmaplist[k]);
																		  return p;}             
								        nextneighbour=true;
								       }
								    }  // for

				                }

						   // find a way to delegate the maximal found subtree to upper level function call
					      
			            }
			      }
		      //if ((nom+newnom)>150) ;
			   //	  printf("\nSearch Space is huge: %d + %d = %d >150",nom,newnom,(nom+newnom));
			  //printf("\n[-4-]Finished one iteration i=%d nextneighbour=%d",i,nextneighbour);
			  if (insidefirst) 
				   first=false;
			  if (nextneighbour) {
				  //printf("\nNext neighbour count=%d",++c);
			      nextneighbour=false;
			      for(int k=0;k<nom;k++) if (maplist[k]) free(maplist[k]);  
				  nom=newnom;
			      newnom=0;
			      for(int k=0;k<nom;k++) {
				     // if (maplist[k]!=NULL) delete maplist[k];  // Need to do proper garbage collection 
				     maplist[k]=newmaplist[k];
			         }
			      }
		      } 
	      }
       //printf("\n[-5-]");
	   if (first) return 0;
	   else if (nom>0) {
           for(int k=0;k<nom;k++) {
             if (map[10]<maplist[k][10])  memcpy(map, maplist[k], 11*sizeof(unsigned int));  
		     }  
	       }
	   for(int k=0;k<nom;k++) if (maplist[k]) free(maplist[k]);
	   return 0;
	  }

unsigned int mapG3( int v, unsigned int map[11], unsigned int protein_no, subG* s, int caller, PROINFO* pht){
		// v: the index of new node in subGraph s that is going to be mapped to protein_no
		// map[] is the set of visited nodes(protein hash table indexes) ordered by their respective indexes in s->labels
		// the result set of proteins is copied to embed[]

       unsigned int *newmap;
	   unsigned int t,p;
	   unsigned int *maplist[200]; 
	   unsigned int *newmaplist[200]; 
       int nom=0;  // the number of maps in maplist
	   int newnom=0;
	   int c=0;
	   bool first=true; bool insidefirst=false; bool nextneighbour=false;
	   map[v]=protein_no;
	   map[10]++;
	   //print_protein_array(map); printf("v=%d c=%d ",v,caller);
       if (map[10]== s->nol)  { 
		   memcpy(embed,map,11*sizeof(int)); 
		   //printf("\nFound embedding"); print_protein_array(embed); 
	       for(unsigned int k=0;k<embed[10];k++)  pht[embed[k]].last_subG = s->id;  //record embedding to prevent overlapings
		   //struct emb* tempem=new struct emb; 
		   //memcpy(tempem->embed ,map,11*sizeof(int)); 
		   //tempem->next=s->embeddings; 
		   //s->embeddings=tempem;
		   return protein_no;
	       }
	   for (int i=0; i< s->nol ; i++ )  {  //for each neighbour i of v (who is not yet mapped)
		   if ( (s->adjm[v][i]==1) &&(map[i]==0)&&(v!=i)) {  
			  for (struct protein* j=pht[protein_no].nl; ((j)&&((nom+newnom)<=150)) ; j=j->next )  {
				   // for each neighbour t of protein (who is not yet mapped) whose label is equal to s->labels[i]           
				   //printf("neighbour name  and labels: %s and %s - %s \n",j->name,pht[get_protein(pht,j->name,hts)].label,s->labels[i]);
				   if ((pht[get_protein(pht,j->name,hts)].label)&&((pht[get_protein(pht,j->name,hts)].label==s->labels[i])||
					                                               (pht[get_protein(pht,j->name,hts)].is_my_label(s->labels[i])))){
					   // if i and j have the same label
					   if (first) {
					        if (!mapped_or_false_connections(i,t=get_protein(pht,j->name,hts),map,s,false,pht)) { 
		                       newmap=(unsigned int*)malloc(11*sizeof(unsigned int));
			                   newmap=(unsigned int*)memcpy(newmap, map, 11*sizeof(unsigned int));    
                               maplist[nom]=newmap; nom++;  // add new mapping to the list
							   if (nom>25) { //printf("\nToo many candidates for a single label %dth=%s..",i,s->labels[i] );
								             break;}
						       //printf("\n [-2-]Added to maplist in first nom= %d protein %d=%d",nom, i,t);
						       if (p=mapG3(i,newmap,t,s,v,pht) > 0) {for(int k=0;k<nom;k++) free(maplist[k]); 
							                                     for(int k=0;k<newnom;k++) free(newmaplist[k]);
							                                     return p;		}
							   insidefirst=true;
						       }
					       }
					   else {
							   for(int k=0;k<nom;k++){
								   if (maplist[k][i]!=0) {maplist[k][i]=0;maplist[k][10]--;}     
								   if (!mapped_or_false_connections(i,t=get_protein(pht,j->name,hts),maplist[k],s,false,pht)) {
								   	    //printf("\nhere comes the fork at protein %d=%d , nom=%d , c=%d, v=%d",i,t,nom,c,v);
		                                newmap=(unsigned int*)malloc(11*sizeof(unsigned int));
			                            newmap=(unsigned int*)memcpy(newmap, maplist[k], 11*sizeof(unsigned int));    
                                        newmaplist[newnom]=newmap; newnom++;  // add new mapping to the list
										if (newnom>25) { //printf("\nToo many candidates for a single label %dth=%s..",i,s->labels[i] );
											             break;}
						                //printf("\n [-3-] Added to maplist %d" , nom);
						                if (p=mapG3(i,newmap,t,s,v,pht) > 0) {for(int k=0;k<nom;k++) free(maplist[k]); 
										                                  for(int k=0;k<newnom;k++) free(newmaplist[k]);
																		  return p;}             
								        nextneighbour=true;
								       }
								    }  // for

				                }

						   // find a way to delegate the maximal found subtree to upper level function call
					      
			            }
			      }
		      //if ((nom+newnom)>150) ;
			   //	  printf("\nSearch Space is huge: %d + %d = %d >150",nom,newnom,(nom+newnom));
			  //printf("\n[-4-]Finished one iteration i=%d nextneighbour=%d",i,nextneighbour);
			  if (insidefirst) 
				   first=false;
			  if (nextneighbour) {
				  //printf("\nNext neighbour count=%d",++c);
			      nextneighbour=false;
			      for(int k=0;k<nom;k++) if (maplist[k]) free(maplist[k]);  
				  nom=newnom;
			      newnom=0;
			      for(int k=0;k<nom;k++) {
				     // if (maplist[k]!=NULL) delete maplist[k];  // Need to do proper garbage collection 
				     maplist[k]=newmaplist[k];
			         }
			      }
		      } 
	      }
       //printf("\n[-5-]");
	   if (first) return 0;
	   else if (nom>0) {
           for(int k=0;k<nom;k++) {
             if (map[10]<maplist[k][10])  memcpy(map, maplist[k], 11*sizeof(unsigned int));  
		     }  
	       }
	   for(int k=0;k<nom;k++) if (maplist[k]) free(maplist[k]);
	   return 0;
	  }



int find_subG( int v, unsigned int *pr, unsigned int protein_no, subG* s){
// find and returns the number of embeddings of s in pht
int label=s->labels[v];
//char *label="0045324";
//printf("\nlabel=%s",label);
int count=0;
unsigned int x[11];
unsigned int positions[1000][11];  // of embedings in pht (the index of s->label[v])

for (unsigned int i=1;((i<hts)&&(count<1000));i++) {
    if ((pht[i].id)&&((pht[i].label==label)||(pht[i].is_my_label(label)))) {
	   for (int j=0; j<11 ; j++) x[j]=0;
	   //count++;
	   if (mapG2(v, x, i, s, 11)>0)  {
           //positions[count]=i;
		   for (int j=0; j<11 ; j++)  positions[count][j]=embed[j];
		   count++; //printf("%d ",count);
	       }
	   }
    }  // for i

s->freq=count;
if (count>=support) { append_to_file(resfile,s,pr,positions); 
                      toplist->add_to_simple_list(s);}
//else if (v==1) append_to_file(resfile,s,pr,positions);  // if a test call
return count;
}

int find_subG3(subG* s,PROINFO* pin){
// find and returns the number of embeddings of s in pin, does not populate embeddings list
int label=s->labels[0];
//char *label="0045324";
//printf("\nlabel=%s",label);
int count=0;
unsigned int x[11];
unsigned int positions[1000];  // of embedings in pht (the index of s->label[v])

for (unsigned int i=1;((i<hts)&&(count<1000));i++) {
    if ((pin[i].id)&&((pin[i].label==label)||(pin[i].is_my_label(label)))) {
	   for (int j=0; j<11 ; j++) x[j]=0;
	   //count++;
	   if (mapG3(0, x, i, s, 11,pin)>0)  {
           positions[count]=i;
		   count++; //printf("%d ",count);
	       }
	   }
    }  // for i
if (pin==pht) { /*printf("true...");*/  s->freq=count;}
return count;
}

int find_subG2(unsigned int *pr, subG* s){
// find and returns the number of embeddings of s in pht, uses embeddings list of s
//int label=s->labels[v];
//char *label="0045324";
//printf("\nlabel=%s",label);
int v;
int count=0;
unsigned int x[11];
unsigned int positions[1000][11];  // of embedings in pht (the index of s->label[0])
unsigned int p;
struct emb* temp;
struct emb* prev=NULL;
if (s->embeddings){  // if there is an embeddings list of s 
for (struct emb* i=s->embeddings ;((i)&&(count<1000)); ) {
  // for each embedding, i->embed will be called with v=s->grow and p=i->embed[s->grow]
	//memcpy(x ,i->embed ,11*sizeof(int));
    v=s->grow;
	p=i->embed[s->grow];
	i->embed[s->grow]=0;
	i->embed[10]--;
	if (mapG2(v, i->embed, p, s, 11)>0)  {
           //positions[count]=i->embed[0];
		   for (int j=0; j<11 ; j++)  positions[count][j]=embed[j];
		   count++; //printf("%d ",count);
           memcpy(i->embed ,embed,11*sizeof(int));
		   prev=i;
		   i=i->next;
	       }
	else {
       if (prev) prev->next=i->next;
	   else s->embeddings=i->next;
       temp=i;
	   i=i->next;
       delete temp;
	   }
    }  // for i
}
else {  // there is no embeddings list of s , therefore do a search on pht
int label=s->labels[0];
for (unsigned int k=1;((k<hts)&&(count<1000));k++) {
    if ((pht[k].id)&&((pht[k].label==label)||(pht[k].is_my_label(label)))) {
	   for (int j=0; j<11 ; j++) x[j]=0;
	   //count++;
	   if (mapG2(0, x, k, s, 11)>0)  {
           //positions[count]=k;
		   for (int j=0; j<11 ; j++)  positions[count][j]=embed[j];
		   
		   count++; //printf("%d ",count);
	       struct emb* tempem=new struct emb; 
		   memcpy(tempem->embed ,embed,11*sizeof(int)); 
		   tempem->next=s->embeddings; 
		   s->embeddings=tempem;
	       }
	   }
    }  // for k
}  // else

s->freq=count;
if (count>=support) { append_to_file(resfile,s,pr,positions); 
                      toplist->add_to_simple_list(s);
                      }
return count;
}

void clean_graph(){
     // make statistics and check the graph for consistency
	int dc=0; int self=0;
	int isolated=0;
	int fisolated=0;
	statistics();
    int nep=0;
	int dnep=0;  // distinct nonexisting proteins
	bool neptrue;
	for(unsigned int i=1;i<hts;i++){
		neptrue=false;
		if (pht[i].id) {
			if (pht[i].non==0) isolated++;
	  	  for (struct protein*	tmp=pht[i].nl; (tmp!=NULL) ; tmp=tmp->next)
            if (hts==get_protein(pht,tmp->name,hts))
			  {nep++; neptrue=true;}
          if (neptrue) dnep++;
		  }
	   }
	 // detects self loops
	for(unsigned int i=1;i<hts;i++){
		if (pht[i].id) {
	  	  for (struct protein*	tmp=pht[i].nl; (tmp!=NULL) ; tmp=tmp->next)
            if (i==get_protein(pht,tmp->name,hts))
			  {self++; printf("\nProblem    [%d]%s--%s,%d ",i,pht[i].id,tmp->name,get_protein(pht,tmp->name,hts));}
		  }
	   }
	printf("\nself:%d  dnep:%d isolated:%d",self,dnep,isolated); // detects self loops


/*	
	t=get_protein(pht,"ydr404c",hts);
	printf("\n ydr404c=%d label=%d",t,pht[t].label );
	for (struct protein*	tmp=pht[t].nl; (tmp!=NULL) ; tmp=tmp->next)
		printf("%s ",tmp->name);
	t=get_protein(pht,"ydl140c",hts);
	printf("\n ydl140c=%d label=%d",t,pht[t].label );
	for (struct protein*	tmp=pht[t].nl; (tmp!=NULL) ; tmp=tmp->next)
		printf("%s ",tmp->name);
	printf("\n pht[3445].id=%s",pht[3445].id );
*/
}


int insert_label_pair(int il1, int il2){
  // into lpt
  int i, temp;
  long int x;
  if (il1) {
     if (il2) {     
		if (il1>il2) {temp=il2;il2=il1;il1=temp;}  // because il1 is the lesser one 
		x=il1+il2;
        //x=(long int)((long int)il2*(long int)il2)+((long int)il1*(long int)il1); // sum of  squares
	    x=x % (pinn*2);
		for(i=(int)x; ((lpt[i].label1)&&(lpt[i].label1!=il1)&&(lpt[i].label2!=il2)) ; i = (i==((pinn*2)-1)) ? 1 : i+1) ;
		if (!(lpt[i].label1)){ lpt[i].label1=il1;
		                       lpt[i].label2=il2;}
		lpt[i].freq++; 
	    }
     else 
	    printf("\n Can not convert label2=%d into int.",il2);
     }
  else 
	  printf("\n Can not convert label1=%d into int.",il1);

  return 0;
}

int frekans(int il1,int il2){
 int x;int i; int temp;
 if (il1>il2) {temp=il2;il2=il1;il1=temp;}  // because il1 is the lesser one 
 x=il1+il2;
 x=x % (pinn*2);
 for(i=(int)x; ((lpt[i].label1)&&(lpt[i].label1!=il1)&&(lpt[i].label2!=il2)) ; i = (i==((pinn*2)-1)) ? 1 : i+1) ;
 if (!(lpt[i].label1)) return 0;
 return lpt[i].freq;  
}

void calc_label_pair_freq(){
   // PRE: lpt is allocated,
	// calculates the frequencies of each label pair
  int l1;
  int l2;
  int temp,j;
  // first insert all pairs into lpt and set their frequencies
  for(unsigned int i=1 ; i<hts ; i++){
	 if (pht[i].id){
       j=0;
	   for (l1=pht[i].label; (l1) ; l1=pht[i].otherlabels[j], j++){  // for the label of pht[i] and otherlabels
  	     for (struct protein* tmp=pht[i].nl; (tmp!=NULL) ; tmp=tmp->next){    // scan the neighbourhood list
           l2=pht[get_protein(pht,tmp->name,hts)].label ;
		   insert_label_pair(l1,l2);
		   } // for tmp
		 if (j>=5) break;
		 }  // for j
	   } // if id
 	}
  // then, insertion sort lpt by freq
  nolp=0;
  for (int i=0; i<(pinn*2) ; i++) if (lpt[i].label1) nolp++;
  soalp=new int[nolp+1];
  int p=0;
  for (int i=0; i<(pinn*2) ; i++) if (lpt[i].label1) {
    soalp[p]=i;
    p++;
    }
  for (int i=0; i<nolp ; i++)
	  for(int j=i+1; j<nolp; j++){
		  if (lpt[soalp[i]].freq > lpt[soalp[j]].freq) {  // swap
             //if (soalp[i]==55262) printf("\n55262.freq=%d , going from i=%d to j=%d",lpt[soalp[i]].freq,i,j);
		     temp=soalp[j];
			 soalp[j]=soalp[i];
			 soalp[i]=temp;
		    }
  	    }
  //for (int i=0; i<nolp ; i++) printf("\ni=%d index=%d freq=%d label1=%d label2=%d",i,soalp[i],lpt[soalp[i]].freq, lpt[soalp[i]].label1, lpt[soalp[i]].label2 ); 
  printf("\nNumber of label pairs:%d",nolp);
  // then set ei to the least freq edge more than support
  for (ei=nolp-1; ei>=0 ; ei--) 
    if (lpt[soalp[ei]].freq < (2*support))	   
	   break;
  if (ei<0) {
     printf("\nNo edge has frequency more than support=%d",support);
     exit(1);
     }
  printf("\nei=%d , We have %d edges with more than support=%d ",ei,(nolp-ei),support);

}

void remove_edge(int e){
   // remove the edge pointed by the e'th element in soalp from pht
	// is it necessary to ?
    int l1=lpt[soalp[e]].label1 ;
	int l2=lpt[soalp[e]].label2 ;
	unsigned int other;
	struct protein* oldnext;
	for(unsigned int i=1;i<hts;i++){
		if(pht[i].id){
			if ((pht[i].label == l1)||(pht[i].is_my_label(l1))) {
				oldnext=NULL;
				for (struct protein* next=pht[i].nl ; (next) ; next=next->next ) {
				   other=get_protein(pht,next->name,hts);
				   if (other>=hts) break;
				   if ((pht[other].label == l2)||(pht[other].is_my_label(l2))) {
				       if (oldnext) oldnext->next =next->next ; 
					   else pht[i].nl= next->next ;
					   pht[i].non --;
				       }
				   else 
				       oldnext=next;
				   }
			   }
			else if ((pht[i].label == l2)||(pht[i].is_my_label(l2))) {
				oldnext=NULL;
				for (struct protein* next=pht[i].nl ; (next) ; next=next->next ) {
				   other=get_protein(pht,next->name,hts);
				   if (other>=hts) break;
				   if ((pht[other].label == l1)||(pht[other].is_my_label(l1))) {
				       if (oldnext) oldnext->next =next->next ; 
					   else pht[i].nl= next->next ;
					   pht[i].non --;
				       }
				   else 
				       oldnext=next;
				   }			   }
		   }

	   }

   }


void remove_proteins(char* fname){  // from pht
   // remove the proteins in the file ff (one line per protein)
	FILE* ff=fopen(fname,"r");
	char buffer[20];
	char protein[10];
    struct protein* oldnext;
    struct protein* temp;
	unsigned int other;
	while (!feof(ff)) {
	   fgets( buffer, 20, ff ); 
       sscanf( buffer, "%10s ", protein );
	   printf("\n%s - %d",protein, strlen(protein));
	   unsigned int index=get_protein(pht,protein,hts);
	   if ((index<hts)) {
	      for (struct protein* next=pht[index].nl ; (next) ; ) {
		     other=get_protein(pht,next->name,hts);
		     if (other>=hts) continue;
			 oldnext=NULL;
			 for (struct protein* onext=pht[other].nl ; (onext) ; onext=onext->next ) {
		        if (strcmp(onext->name,protein)==0) {
		          if (oldnext) oldnext->next =onext->next ; 
			      else pht[other].nl= onext->next ;
			      pht[other].non --;
                  // need to do proper deletion of onext here 
		          }
		        else 
		          oldnext=onext;
		        }  // for onext
			 temp=next;
			 next=temp->next;
             delete temp;
		     }  // for next 
		  pht[index].nl=NULL;
		  pht[index].non=0;
		  // we do not completely delete the pht entry here because it causes mayhem if we do so, 
		  //  hash values may become invalid
	      }  // if index

	   }
	//printf("\nid:%s,non:%d",pht[get_protein(pht,"ypl012w",hts)].id ,pht[get_protein(pht,"ypl012w",hts)].non);
	fclose(ff);
   }

void record_s(subG* s){}
/*
int *find_rightmost_path(subG *s){
   // returns the right-most path as an array of int[11] where int[10] holds the length
   // the rigthmost child is the one in int[10]'th posisiton, therefore root is int[0] 
   // the root of s is the 0'th node
   // algorithm: 
   //   follow the rightmost path until either you hit the end or you hit the root!
   int *rmpath= new int[11];
   rmpath[0]=0; rmpath[10]=1;
   bool have_choice=false;
   short int pivot=0;
   short int np;  // next pivot
   for (np=(s->nol-1) ; !(s->adjm[pivot][np]) ; np-- );
   int pc=1;  // pivot count
   while ((pc<=10)&&(not_loop(rmpath,np)&&(!leaf(np))) {  // while there is no loop
      pivot=np;
	  rmpath[pc]=pivot;
	  pc++;
      rmpath[10]++;
      for (np=(s->nol-1) ; !(s->adjm[pivot][np]) ; np-- );
      }
   if (pc>10) {
      rmpath[10]=1;
      }
   else if leaf(np) {

      }
   else { // loop
      
 
      }

   }
*/
/*
int *find_rightmost_path(subG *s,int* rmpath){
	// pivot is the last entry into the rmpath
	// the right most path of the tree rooted at pivot is found and returned
	// rmpath's length is at least 2,  rmpath[0]=0,
	if (rmpath[10]>=10) { 
	    rmpath[10]=1;
		return rmpath;
	   }
	int pivot=rmpath[rmpath[10]-1];
	int parent=rmpath[rmpath[10]-2];
    int np,i;  // next pivot
	for (np=(s->nol-1) ; ((np>=0)&&((!(s->adjm[pivot][np]))||(np==parent))) ; np-- );
	if (np<0) return rmpath;  // pivot is a leaf node
	else{
	   //if is_loop(rmpath,np) {}
       for(i=0; i< rmpath[10] ; i++) 
		   if (rmpath[i]==np) 
			   break;
	   if (i< rmpath[10]) {  // there is a loop starting from ith element in rmpath (rmpath[i]==np)
          rmpath[10]=i+1;
		  return rmpath;
	      }
	   else {
          rmpath[rmpath[10]]=np;
		  rmpath[10]++;
		  return find_rightmost_path(s,rmpath);
	      }
	   }

}

*/
class tnode{
  public:
	 int index; // index of the node in the subg
     tnode* child;  // points to first child
	 tnode* next;  // points to next sibling
     int noc;  // nmber of children 
	 int diorno;   // vertex discovery order number v0,v1,v2 during create_DFS_tree()..
	 //subG* s;  // subG that this tree is constucted from
	 static int visited[11];  // the visited nodes of a subgraph
     tnode(){noc=0;child=NULL;index=0;next=NULL;diorno=0;}
	 tnode(int ix,int no){noc=0;child=NULL;index=ix;diorno=no;}
	 tnode* add_child(int i,int dio)
	    { tnode* t=new tnode(); t->next = child; child=t; noc++; t->index=i;t->diorno=dio; return t;}
	 void del_tree(tnode *tree){
		 tnode* temp;
		 for (tnode* t=tree->child ; (t) ;  t=temp) {temp=t->next ; del_tree(t);}
		 delete tree;	
	    }
	 tnode* create_DFS_tree(subG* s,int* nextno,DFScode* code,int diornoa[10]){
		 // visited[10]=0 must be set before calling
		 // this also creates the DFS code and puts it into code,
		 // diornoa is used to hold discovery order numbers of all vertices
		 // nextno is the next diorno to be assigned to the next vertex
		 visited[visited[10]]=index;
		 visited[10]++;
		 diorno=(*nextno);   // diorno starts at 0
 		 (*nextno)++;  
		 diornoa[index]=diorno;
		//index=ix; 
	    //tnode* t=new tnode();
        // starting with the 0th node in s as the root
        int np;  // next pivot
		tnode* c;
		for (np=(s->nol-1) ; (np>=0) ; np-- ) 
			if (s->adjm[index][np]==1) 
				if (not_visited(np)) {
					c=add_child(np,(*nextno));
					diornoa[np]=(*nextno);
					// add np to DFS code
					code->cedges[code->noedges].v1=diorno;
                    code->cedges[code->noedges].v2=c->diorno;
                    code->cedges[code->noedges].ei=get_ei(s->labels[index],s->labels[c->index]);
					code->noedges++;
					// add backward links from c to DFS code
		            for (int j=(s->nol-1) ; (j>=0) ; j-- ) 
						if ((s->adjm[c->index][j]==1)&&(j!=index)&&(!not_visited(j))){
					       code->cedges[code->noedges].v1=c->diorno;
                           code->cedges[code->noedges].v2=diornoa[j];
						   code->cedges[code->noedges].ei=get_ei(s->labels[c->index],s->labels[j]);
					       code->noedges++;
						  }
					c->create_DFS_tree(s,nextno,code,diornoa);
					}
        return this;
	    }
	 bool not_visited(int np){
		 for (int i=0; (i<visited[10]) ; i++) 
			 if (visited[i]==np) return false;
	     return true;
	     }
	 void rmpath(){  // printd the right most path
          // the visited[] must be empty before first call
        //printf("||%d ",index);
		visited[visited[10]]=index;
        visited[10]++;
		tnode* oldt=NULL;
        for (tnode *t=child; t ; t=t->next) oldt=t;
		if (oldt) oldt->rmpath(); 
	    }

	 void print_tree(){  // rooted at this
         printf("\n%d -> ",index);
		 tnode* oldt=NULL;
		 for (tnode *t=child; t ; t=t->next) {printf("%d ",t->index);}
		 for (tnode *t=child; t ; t=t->next) t->print_tree(); 
	    }

	 int get_ei(int l1,int l2){
        int i;
        for (i=ei; i<nolp ; i++)
           if (((lpt[soalp[i]].label1 == l1)&&(lpt[soalp[i]].label2 ==  l2))|| 
               ((lpt[soalp[i]].label2 == l1)&&(lpt[soalp[i]].label1 ==  l2)) )
               break;
        if (i<nolp)
		       return i;
		   else 
		      printf("\nCould not find a matching edge for DFScode edge %d-%d , ei:%d",l1,l2,ei);
        return i;
	    }
/*
	 DFScode_edge* DFScode(subG* s, DFScode* code){  // DFS code of s rooted at this
         // code is the accumulated DFC code up to this,
		 // diorno is the discovery order number as v1,v2,v3.. in gSpan
		 printf("$%d$%d# ",index,diorno); 
		 for (tnode *t=child; t ; t=t->next){
			 
			 t->DFScode(s,code); 
		     }
		 return NULL;
	    }
*/
	 int permutations_of_not_yet_visited_neighbours(int perms[1000][11],int wisited[11],subG* s){
        // fills perms with permutations of not yet visited neighbours of this
        // according to the wisited array
        // returns the number of permutations
		int items[11];
		int n=0;  // number of items to be permuted
		for (int np=(s->nol-1) ; (np>=0) ; np-- ) 
			if ((s->adjm[index][np]==1)&&(index!=np)&&(!wisited[np])) {
				items[n]=np;
		        n++;
	            }
        if (n==0) return 0;
		else if (n==1) {
            perms[0][0]=items[0];
			perms[0][10]=1;
			return 1;
		    }
		else {  // there are more then one items to generate permutations of
            items[10]=n;
			// next_permutation( items, items);
            return generate_perms(items,perms);
		    }
	    }  // permutations_of_not_yet_visited_neighbours

	int generate_perms(int items[11],int perms[1000][11]){
        // aii: add item index

        for(int i=0; (i<items[10]) ; i++) perms[0][i]=items[i];
		for(int i=items[10]-1; (i>=0) ; i--) perms[1][items[10]-i-1]=items[i];
		perms[0][10]=items[10];
        perms[1][10]=items[10];
		return 2;

		int temp[1000][11];
		int nop=1; // number of perms
        perms[0][0]=items[0]; 
		perms[0][0]=1;
		int n=0;
		for (int aii=1 ; (aii<items[10]) ; aii++) { // for each item 
            n=0;
			for (int i=0; (i<nop) ; i++) {  // add this item into existing permutations
				for (int j=0; (j<aii) ; j++) { // for each existing permutation
					                                    // insert item[aii] into positions
                   temp[n][j]=items[aii];		   
				   }  // for j
			   }
		   }
	    } 

	bool min_s(subG* s, DFScode* code, int *pivot, tnode* root, int wisited[11],int diornoa[10],int* nextno){
       // is s==min(s) ?
       // pivot: that holds the index up to which the DFS code of this is identical to that of subG's 
       // code is the DFS code of s
       // diornoa is the array of discovery order number of vertices
       // need to set wisited to 0's before calling
	   //printf("\nVisiting node at index %d , discovery order no:%d",index,(*nextno));
       wisited[index]=1;
	   diorno=(*nextno);   // diorno starts at 0
	   if (diorno>=s->nol) return true;
 	   (*nextno)++;  
	   diornoa[index]=diorno;
       int perms[1000][11];
	   perms[0][10]=0;
       int n=permutations_of_not_yet_visited_neighbours(perms,wisited,s);
	   tnode* c;  // child
	   (*pivot)++; 
	   if ((*pivot)>=code->noedges) return true;
	   int tpivot=(*pivot);
	   int wwisited[11];
	   memcpy(wwisited,wisited,11*sizeof(int));
       int ddiornoa[10];
	   memcpy(ddiornoa,diornoa,10*sizeof(int));
	   for (int i=0 ; i<n ; i++) {  // for each permutation of nonvisited neighbours
		  // for each of the nonvisited children c, check c.min_s
          nextperm:
          memcpy(wisited,wwisited,11*sizeof(int));
		  memcpy(diornoa,ddiornoa,10*sizeof(int));
          (*nextno)=diorno+1;
		  (*pivot)=tpivot;
 		  for (int j=0; j<perms[i][10] ; j++){
		 	 if (!wisited[perms[i][j]]){       
                // check if forward edge from index to perms[i][j] has link equality
                if (diornoa[index]<code->cedges[(*pivot)].v1) return false;
		        else if (diornoa[index]==code->cedges[(*pivot)].v1){
			        if ((*nextno)<code->cedges[(*pivot)].v2) return false;
			        else if ((*nextno)==code->cedges[(*pivot)].v2) {
				       int eii=get_ei(s->labels[index],s->labels[perms[i][j]]);
				       if ( eii < code->cedges[(*pivot)].ei)
					      return false;
					   else if (eii == code->cedges[(*pivot)].ei){  // next DFC code lines are equal 
						  // check for backward link equality
                          // Her bir backward link permütasyonu ayri bir DFS koddur aslinda
						  diornoa[perms[i][j]]=(*nextno);
						  for (int k=(s->nol-1) ; ((k>=0)) ; k-- ){ 
						     if ((s->adjm[perms[i][j]][k]==1)&&(k!=perms[i][j])&&(index!=k)&&(wisited[k])){
                               (*pivot)++;
                               if (diornoa[perms[i][j]]<code->cedges[(*pivot)].v1) return false;
		                       else if (diornoa[perms[i][j]]==code->cedges[(*pivot)].v1){
			                       if (diornoa[k]<code->cedges[(*pivot)].v2) return false;
			                       else if (diornoa[k]==code->cedges[(*pivot)].v2) {
				                      int eii=get_ei(s->labels[perms[i][j]],s->labels[k]);
				                      if ( eii < code->cedges[(*pivot)].ei)
					                     return false;
					                  else if (eii == code->cedges[(*pivot)].ei){  // next DFC code lines are equal 

									     }
									  else {goto nextperm;}
									  }
								   else {goto nextperm;}
							       }
							   else {goto nextperm;}
							   }
						     }  // for k
                          // Go deeper
						  c=add_child(perms[i][j],(*nextno));
			              if (c->min_s(s,code,pivot,root,wisited,diornoa,nextno))
							 return true;
						  else
				             return false;
					      }  
                          
					   else break;
					   }
					else break;
				    }
				else break;
			    }  // if not wisited
		     }  // for j
	      }  // for i
	   return true;
	   }  // min_s

} ;  // tnode

bool minCanonnicalCode(subG* s){
     tnode* tt=new tnode();
	 DFScode* dc=new DFScode();
	 int no=0; 
	 int diornoa[10];
	 tt->create_DFS_tree(s,&no,dc,diornoa);
	 //for(int i=0;i<dc->noedges; i++) 
     //		 printf("\nCode edge:%d,%d,%d",dc->cedges[i].v1 ,dc->cedges[i].v2, dc->cedges[i].ei);
     
	 for (int m=0 ; (m<s->nol) ; m++){
	   tt=new tnode(m,0);
	   int pivot=-1;
	   int wisited[11]={0,0,0,0,0,0,0,0,0,0,0};
	   int diornoaa[10]={0,0,0,0,0,0,0,0,0,0};
	   int nextno=0;
	   bool mins=tt->min_s(s,dc,&pivot,tt,wisited,diornoaa,&nextno);
	   if (mins) ; //printf("\nTrue");
	   else if (!mins) { return false;}
	 }	
	return true;
  }

int tnode::visited[11]; 

tnode* rmpath(subG* s){
   // output: tnode->visited holds the rm path
   tnode *tt=new tnode();
   tt->visited[10]=0;
   tt->index =0;
   int no=0;
   DFScode* dc=new DFScode();
   int a[10];
   tt->create_DFS_tree(s,&no,dc,a);
   tt->visited[10]=0;
   tt->rmpath();  // put the result in visited[]
   return tt;
   }

int copy_embeddings(subG* news, subG* s, int bindex, int index){
    // copies the liked lists from s to news while erasing the index and bindex entries 
    struct emb* t;
	news->embeddings = NULL;
    for (struct emb* temp=s->embeddings; (temp) ; temp=temp->next  ) {
	 	t=new struct emb;
        memcpy(t->embed ,temp->embed ,11*sizeof(int));
		//t->embed[index]=0; t->embed[10]--;
		if (bindex!=index) {t->embed[bindex]=0; t->embed[10]--;}
		else {t->embed[news->nol-1]=0; }
		t->next=news->embeddings;
		news->embeddings=t;
	    }	
	return 1;
   }

subG* grow(subG* s,int index){   // for node at index; generate all possible forward extensions 
                                 // with all valid labels and append it to the children list 
 //  returns a linked list of subG's whom are the childre of s
  int extlabel;  // extending label
  subG* ss=NULL;  // head of linked list  
  subG* temp=NULL;
  
  for(int i=nolp-1; i>=ei ; i--){
     // first find the label of the extending node
	 int label=s->labels[index];
	 if (lpt[soalp[i]].label1==label){
        extlabel=lpt[soalp[i]].label2;
	    }
	 else if (lpt[soalp[i]].label2==label){
        extlabel=lpt[soalp[i]].label1;
	    }
	 else
		 continue;
     long temp;
     // then create the child tree
	 subG *news=new subG;
	 temp=news->id;
     memcpy(news,s,sizeof(subG));
	 news->id=temp;
     news->labels[news->nol]=extlabel;
     news->adjm[news->nol][index]=1;
     news->adjm[index][news->nol]=1;
     news->nol++;
	 news->grow=index;
	 copy_embeddings(news,s,index,index);   
     // prepend news to ss
	 news->next=ss;
	 ss=news;
	 }

  return ss;

}

subG* back_link(subG* s,int bindex, int index){
 // create a backward link: link bindex and index (right most node)
 //  returns a linked list of subG's whom are the childre of s
 //  if there are edge labels then all possible labels must be tested
 //  for now, there is no edge label therefore only one subG wil be created
   long temp;
   int fr=frekans(s->labels[bindex],s->labels[index]);
   if (fr<(2*support)) return NULL;
   if (lpt[soalp[ei]].freq > fr) return NULL;
   if (s->adjm[index][bindex]==1) return NULL;
   if (index==bindex) return NULL;
   if ((index>10)||(bindex>10)||(bindex<0)||(index<0)) return NULL;
	  subG *news=new subG;
	  temp=news->id;
      memcpy(news,s,sizeof(subG));
	  news->id=temp;
      news->adjm[index][bindex]=1;
      news->adjm[bindex][index]=1;
	  news->next=NULL;
	  news->grow=index;
	  copy_embeddings(news,s,bindex,index);
      return news;

   }

subG* generate_children(subG* s){  // returns a linked list of subGs
	if (s->nol >= 10) {printf("\n Size of subG is 10 or more"); return NULL;}
	subG* ss=NULL;
    subG* children=NULL;
	subG* news;
	subG* temp;
	subG* oldtemp=NULL;
	tnode* tt=rmpath(s);
	
	for(int i=0 ; i < (tt->visited[10]) ; i++) {  // for each node on rmpath starting from root
       // create a forward link: grow an edge from the node visited[i]
	   if (news=grow(s,tt->visited[i])) { // generate all possible extensions with all valid labels
	      for (temp=news,oldtemp=NULL; temp ; temp=temp->next) oldtemp=temp;  // oldtemp points to the last node
          if (oldtemp) oldtemp->next=children;  // and append it to the children list
	      else printf("\nSomething Wrong with growing edge");
	      children=news;
	      }       
	   }
	
	for(int i=tt->visited[10]-3 ; i >= 0 ; i--) {  // for each node on rmpath starting from rightmost node
       // create a backward link: link visited[i] and visited[tt->visited[10]-1] (right most node)
	   if (news=back_link(s,tt->visited[i],tt->visited[tt->visited[10]-1])) { // generate all backlinks with all valid edge labels
	      for (temp=news,oldtemp=NULL; temp ; temp=temp->next) oldtemp=temp;  // oldtemp points to the last node
          if (oldtemp) oldtemp->next=children;  // and append it to the children list
	      else printf("\nSomething Wrong with back growing edge");
	      children=news;
	      }       
	   }
	tt->del_tree(tt);
	return children;
    }


void subgraphs(subG* s,int depth,int count) {
	//find all frequent subgraphs whose root is s in the candidate tree
	// s is a frequent subgraph 
	// this procedure is similar to Subprocedure 2 (Subgraph Mining ) in gSpan
   int count_child;
   //subG* sxs=new subG(*s);
   if (depth>maxdepth) { 
	   append_to_file(maxresult_file ,s , emptyarray,NULL);
	   return;
       }
   if (!minCanonnicalCode(s)){
	   pruned++;
       }
   record_s(s);  // record s as a legitimate result
   subG* oldchild=NULL;
   bool maximal=true;
   for (subG* child=generate_children(s); child ; ) { 
     // for each child graph
       //printf("\ns->id:%d",child->id );
	   if (use_embeddings) count_child=find_subG2(emptyarray,child);  // count the embeddings
	   else	               count_child=find_subG(0 ,emptyarray,0 ,child);  // count the embeddings
	   if (count_child >= support) { 
	      printf("-*-%d",count_child);
          subgraphs(child,depth+1,count_child);
	      maximal=false;
	      }
	   oldchild=child;
	   child=child->next;
	   delete oldchild;
       }
   if ((maximal)&&(s->nol>=3)) 
       //compute_Z_score(s);
       //printf("\n Zscore: %f , Freq:%d", s->Zscore , s->freq );
	   //sxs= new subG(); (*sxs)=(*s);
	   append_to_file(maxresult_file ,s , emptyarray,NULL); // this is a maximal freq subgraph
       //toplist->add_to_list(s);
	   //toplist->add_to_simple_list(s);
   }

void all_subgraphs(){
  // for each frequent edge ei; starting with the lowest frequent one
  //    find all frequent subgraphs that contain ei
  //    remove ei from pht. PRE: ei is set.

  for (ei++; ei<nolp ; ei++) {  // for each frequent edge ei; starting with the lowest frequent one
     // convert ei into subG
     printf("\nChecking edge[%d,%d]",lpt[soalp[ei]].label1,lpt[soalp[ei]].label2);
	 subG *s=new subG;
	 s->nol = 2;  s->freq=lpt[soalp[ei]].freq/2;
	 s->labels[0]=lpt[soalp[ei]].label1 ;
	 s->labels[1]=lpt[soalp[ei]].label2 ;	   
	 s->adjm[0][1]=1; s->adjm[0][0]=0;
     s->adjm[1][0]=1; s->adjm[1][1]=0;
	 subgraphs(s,2,lpt[soalp[ei]].freq); //find all frequent subgraphs whose root is s in the candidate tree
	 remove_edge(ei); // remove ei from pht, check if there is any prospect for more freq subgraphs
	 	
     }
  // print results
  FILE* plf=fopen(rankfile, "a"); fprintf(plf,"\n%s\n-------------",command);
  for (subG* temp=toplist->list ; (temp) ; temp=temp->next) {
	  display_subG(temp,new unsigned int[11],plf);
	  fprintf(plf,"\n------------pruned:%d----",pruned);
      }
  fclose(plf);

}

int read_subG_find(char* ff, PROINFO* pin){ 
   int x=0;
   subG* olds=NULL;
   unsigned int ppr[11]; for(int k=0; k<11;k++)  ppr[k]=0;   
   subG* s=new subG();
   read_subG(s,ppr,ff);  // reads all subgraphs in file ff and returns them in a linked list pointed to by s
   if (s->is_subgraph(s->next)) printf("\nYes Subgraph");
   else printf("\nNot Subgraph");
   //int perms[1000][11]; int wisited[11]; for(int k=0;k<11;k++)wisited[k]=0;
   //int nn;
   //return find_subG(1 ,ppr,ppr[0] ,s);
   for (subG* temp=s; (temp) ; temp=temp->next) {
     if (olds) delete olds;
     x=find_subG3(temp,pin);
	 printf("\n This subG has %d embeddings and its Z-score is %f .",x,temp->Zscore=compute_Z_score(temp));
	 display_subG(temp,ppr,stdout);
	 tnode *tt=new tnode();
	 olds=temp;
	 
	 tt->visited[10]=0;
	 DFScode* dc=new DFScode();
	 int no=0; int diornoa[10];
	 tt->create_DFS_tree(temp,&no,dc,diornoa);
	 int m=0;
	 for(int i=0;i<dc->noedges; i++) 
		 printf("\nCode edge:%d,%d,%d",dc->cedges[i].v1 ,dc->cedges[i].v2, dc->cedges[i].ei);
     
	 for (m=0 ; (m<temp->nol) ; m++){
	 tt=new tnode(m,0);
	 int pivot=-1;
	 int wisited[11]={0,0,0,0,0,0,0,0,0,0,0};
	 int diornoaa[10]={0,0,0,0,0,0,0,0,0,0};
	 int nextno=0;
	 bool mins=tt->min_s(temp,dc,&pivot,tt,wisited,diornoaa,&nextno);
	 if (mins) printf("\nTrue");
	 else if (!mins) printf("False");


	 }	 
	 //		        lpt[soalp[dc->cedges[i].ei]].label1 ,
	 //			lpt[soalp[dc->cedges[i].ei]].label2
	 //tt->DFScode(temp,NULL);
/*
	 wisited[1]=1;
	 nn=tt->permutations_of_not_yet_visited_neighbours(perms,wisited,temp);
	 printf("\nPermutations:%d",nn);
	 printf("\n%d %d",perms[0][0],perms[0][1]);
	 printf("\n%d %d",perms[1][0],perms[1][1]);
*/
     }

   return x;
}

void exclusions(){
  exclude=1;
  for (int i=ei; i<nolp ; i++){
	 if ((lpt[soalp[i]].label1 == 3735)&&(lpt[soalp[i]].label2 == 3735)) {
        remove_edge(i);
	    }
	 else if ((lpt[soalp[i]].label1 == 3723)&&(lpt[soalp[i]].label2 == 3723))
		remove_edge(i);
	 else if ((lpt[soalp[i]].label1 == 3674)&&(lpt[soalp[i]].label2 == 3674))
		remove_edge(i);
     }
  printf("\nCertain edges were excluded:(3735-3735)(3723-3723)(7046-7046)");

}

void test7(){
/*
	subG s;
s.nol=2;
s.next =NULL;
s.labels[0]=1111;
s.labels[1]=2222;
s.adjm [0][1]=1;
s.adjm [1][0]=1;
int rp[11];
rp[10]=2;
rp[0]=0;
rp[1]=1;
find_rightmost_path(&s,rp);
for (int i=0; i<rp[10] ; i++) printf("|%d",rp[i]);
*/
/*  
	subG* s=new subG;
    unsigned int pr[11]; for(int k=0; k<11;k++)  pr[k]=0; 
	read_subG(s,pr,"subG.txt");
	int rp[11];
    rp[10]=2;
    rp[0]=0;
    int np,i;  // next pivot
	for (np=(s->nol-1) ; (!(s->adjm[0][np])) ; np-- );
    rp[1]=np;
//	find_rightmost_path(s,rp);
    for (int i=0; i<rp[10] ; i++) printf("|%d",rp[i]);
*/
/*
tnode *tt=new tnode();
tt->visited[10]=0;
//tt.s=s;
tt->index =0;
tt->create_DFS_tree(s);
tt->visited[10]=0;
tt->rmpath();
*/
  
/*	
	for (ei=nolp-1; ei>=0 ; ei--) 
    if (lpt[soalp[ei]].freq < (2*support))	   
	   break;
  if (ei<0) {
     printf("\nNo edge has frequency more than support=%d",support);
     exit(1);
     }
  printf("\nWe have %d edges with more than support=%d ",(nolp-ei),support);

tnode* tt=rmpath(s);
tt->print_tree();
//tt->del_tree(tt);
int count=0;
//subG *temp=grow(s,tt->visited[0]);
subG *temp=generate_children(s);
while (temp) {
   display_subG(temp,pr,stdout);
   temp=temp->next;
   count++; 
   if (count>10) break;
   }
printf("\nCount:%d",count);

printf("\nFrekans(3682,16563)=%d",frekans(3682,16563));
printf("\nFrekans(7064,16563)=%d",frekans(7064,16563));
printf("\nFrekans(5515,16563)=%d",frekans(5515,16563));
*/
/*
for(int j=0; j<(pinn*2) ; j++)
   if ((lpt[j].label1 == 3682)||(lpt[j].label2 == 3682))
	   printf("{label1:%d label2=%d}",lpt[j].label1 , lpt[j].label2);

*/



//exit(1);

}

void test6(){
 int histo[100];
int v; unsigned int p; 
 subG* ss;
 int found=0;
 int total_found=0;
 int not_found=0;
 int not_verified=0;
 int ggfailed=0;
 int j;
 int nn;
 int res,itn;
 for (int i=0;i<100;i++) histo[i]=0;
 printf("\nGive the number of iterations, nodes and support:");
 res=scanf("%d %d",&itn,&nn, &support);
 printf("\n nn=%d itn=%d",nn,itn);
 for (j=0 ; total_found<itn ; j++) {
   unsigned int x[11]; for (int i=0; i<11 ; i++) x[i]=0;
   unsigned int pr[11]; for(int k=0; k<11;k++)  pr[k]=0;   
   unsigned int ppr[11]; for(int k=0; k<11;k++)  ppr[k]=0;   
   ss=random_subgraph_generator(nn, &v,&p,pr,&v,&v);			
   if (ss!=NULL) {
	        //display_subG(ss,pr,stdout);
			//FILE* subG3=fopen("subG3.txt","w");
			//display_subG(ss,pr,subG3);
			//fclose(subG3);
            //subG* s=new subG;
            //read_subG(s,ppr,"subG2.txt");
      //display_subG(ss,pr,stdout);
	  //printf("\nDo you want to explore this pattern? (Y or N)");
      //scanf("%s",choice);
	  //if (strcmp(choice,"Y")) continue; 
	  //else {
	    //if (!verify_sub(pr,s,false)) {printf("\nCould not verify generated embedding "); not_verified++;}    
//		if (found=find_subG(1 ,ppr,ppr[0] ,s)>0)  {
		found=find_subG3(ss,pht);
	    if (found>0)  {
			histo[found]++;
		    //printf("Found %d embeddings in pht",found);
			 printf("\nHISTOGRAM: ");
             for (int i=0;i<100;i++) if (histo[i]) printf("%d=%d ",i,histo[i]);
             printf("\n");
			total_found++;
		    if (!verify_sub(embed,ss,false)) {printf("\nCould not verify found embedding "); not_verified++;}
	       }
	    else {
		    not_found++;
		    printf("\n Could not find any embedding. Shame on you."); 
			//display_subG(ss,pr,stdout);
			//exit(1);
			/*
			getc(stdin);
			FILE* subG2=fopen("subG2.txt","w");
			display_subG(ss,pr,subG2);
			fclose(subG2);
            subG* s=new subG;
			for(int k=0; k<11;k++)  ppr[k]=0;   
	        read_subG(s,ppr,"subG2.txt");
			for (int i=0; i<11 ; i++) x[i]=0;
			if (mapG2(v, x , ppr[v] ,s, 11)>0)
				printf("\nSuccessfull in the second run");
			else 
			    exit(1);
			getc(stdin);
			*/
	        }
	     //}  // ELSE
      }
   else {
	   //printf("\nRandom Graph Generator failed.."); 
	   ggfailed++;}
 }  // for j
 printf("\n%d Runs resulted with total found=%d , not_verified=%d , not_found=%d ggfailed=%d\n",j,total_found, not_verified,not_found,ggfailed);
 printf("HISTOGRAM: ");
 for (int i=0;i<100;i++) if (histo[i]) printf("%d=%d ",i,histo[i]);
 printf("\n");
}

void test5(){
    subG* s=new subG;
    unsigned int x[11]; for (int i=0; i<11 ; i++) x[i]=0;
    unsigned int pr[11]; for(int k=0; k<11;k++)  pr[k]=0; 
    //int nom=0;
    //unsigned int *mappings[200];
	read_subG(s,pr,"subG3.txt");
    if (s!=NULL) {
	 for (int i=0; i<s->nol ; i++){
	   for (int j=0; j<11 ; j++) x[j]=0;
	   if (mapG2(i, x, pr[i] ,s, 11)>0)  {
		  if (verify_sub(embed,s,true)==0) {printf("\nCould not verify found embedding "); }
	      }
	   else {
		  printf("\n Could not find embedding at the specified point..");  exit(1);
	      }
	   }
     }
    // Test Neighbout list
    /* 
	xx=get_protein(pht,"ybr189w",hts);
	if (xx<hts) printf("label of %s is %s \n",pht[xx].id, pht[xx].label);
	printf("Protein %s's one neighbur is %s \n",pht[xx].id , pht[xx].nl->name);
	int c=0;
	for (struct protein* t=pht[xx].nl  ; t ; t=t->next , c++) printf("%s:%4f  ",t->name, t->conf);
	printf("\n neighbour count for %s is %d\n",pht[xx].id,c);	
	*/
	// Test OK   
	exit(1);
  }

void random_graph_statistics(){ 
	// counts subgraphs randomly selected from pht
 int histo[10][100][100]; //  histogram of subgraphs[# of nodes][# of edges][# of embeddings]
 int v; unsigned int p; 
 int nnodes, nedges;
 subG* ss;
 int found=0;
 int total_found=0;
 int not_found=0;
 int not_verified=0;
 int ggfailed=0;
 int j;
 int nn;
 int res,itn;
 for (int s=0;s<10;s++) for (int q=0;q<100;q++)   for (int r=0;r<100;r++) histo[s][q][r]=0;
 printf("\nGive the number of iterations, nodes and support:");
 res=scanf("%d %d",&itn,&nn, &support);
 printf("\n nn=%d itn=%d",nn,itn);

 for (j=0 ; total_found<itn ; j++) {
   unsigned int x[11]; for (int i=0; i<11 ; i++) x[i]=0;
   unsigned int pr[11]; for(int k=0; k<11;k++)  pr[k]=0;   
   unsigned int ppr[11]; for(int k=0; k<11;k++)  ppr[k]=0;   
   ss=random_subgraph_generator(nn, &v,&p,pr,&nnodes,&nedges);	
   if (ss!=NULL) {
		found=find_subG3(ss,pht);
	    if (found>0)  {
 			 histo[nnodes][nedges][found]++;
			 total_found++;
		     if (!verify_sub(embed,ss,false)) {printf("\nCould not verify found embedding "); not_verified++;}			 
			 if (found > 1) {
				FILE* subG3=fopen("random_graphs.txt","a");
				fprintf(subG3,"\nThis subG contains %d nodes, %d edges and has %d embeddings in PIN:",nnodes,nedges,found);
			    display_subG(ss,pr,subG3);
			    fclose(subG3);
			    }
			 if (((total_found % 100)==0)||(total_found>=itn)){
				 FILE* subG3=fopen("random_graphs_histo.txt","a");
				 fprintf(subG3,"\niteration:%d , max nodes=%d , level=%d",itn,nn, tlevel);
                 fprintf(subG3,"\n%d Runs resulted with total found=%d , not_verified=%d , not_found=%d ggfailed=%d",j,total_found, not_verified,not_found,ggfailed);
                 fprintf(subG3,"\nHISTOGRAM:\n ");
                 for (int p=0;p<10;p++) 
                    for (int q=0;q<100;q++)
                       for (int r=0;r<100;r++) 
						   if (histo[p][q][r]) fprintf(subG3,"histo[\t%d\t][\t%d\t][\t%d\t]=\t%d\n",p,q,r,histo[p][q][r]);
			     fclose(subG3);
				//printf("\nHISTOGRAM: \n");
				//for (int i=0;i<10;i++) { 
				//   printf("\n%d-node graphs: ",i);
				//   for (int j=0;j<100;j++) if (histo[i][0][j]) printf("%d=%d ",j,histo[i][0][j]);
				//   }
                //printf("\n");
				}
	        }
	    else {
		    not_found++;
		    printf("\n Could not find any embedding. Shame on you."); 
	        }
      }
   else {  ggfailed++;}
 }  // for j
 printf("\n%d Runs resulted with total found=%d , not_verified=%d , not_found=%d ggfailed=%d\n",j,total_found, not_verified,not_found,ggfailed);
 printf("\nHISTOGRAM: \n");
 for (int p=0;p<10;p++) 
   for (int q=0;q<100;q++)
     for (int r=0;r<100;r++) if (histo[p][q][r]) printf("histo[\t%d\t][\t%d\t][\t%d\t]=\t%d\n",p,q,r,histo[p][q][r]);
}

PROINFO* copyPIN(PROINFO* pht){
  // creates a copy of pin pht of size hts

  PROINFO* pin=new PROINFO[hts];
  memcpy(pin,pht,sizeof(PROINFO)*hts);
  for (unsigned int i=1;((i<hts));i++) {
    //char *tempc=new char[10];
    if (pht[i].id==NULL) continue;
	pin[i].id=strcpy(new char[10],pht[i].id);
	pin[i].label = pht[i].label ;
    pin[i].non = pht[i].non ;
	pin[i].last_subG=0;
	memcpy(pin[i].otherlabels,pht[i].otherlabels,sizeof(int)*5);
	// copy neighbour list
	struct protein*	tmpd;
    struct protein*	oldtmpd;
	if (pht[i].nl==NULL) pin[i].nl=NULL;
	else {  
      pin[i].nl=tmpd=new struct protein; 
	  for (struct protein*	tmp=pht[i].nl; (tmp!=NULL) ; tmp=tmp->next){
		strcpy(tmpd->name,tmp->name);
		tmpd->conf = tmp->conf ;
		oldtmpd=tmpd;
	    tmpd=new struct protein; 
		oldtmpd->next=tmpd;
	    }
	  oldtmpd->next=NULL;
	  delete tmpd;
	  } // else 

    } // for 
  return pin;
}

bool find_neighbour(PROINFO* pin, unsigned int i1, unsigned int j1, unsigned int* ii2, unsigned int* jj2){
  // returns an edge (i2-j2) randomly from pin where the closer i2-j2 to i1-j1, the higher their chance of being selected.
  // Close neighbours are favored.
  // First determine a random distance between i1-j1 and i2-j2.
  // target distance d is an exponentialy distributed random variable
  int RADIUS=10;
  double lambda=0.2; 
  int d=(((int)(-1*(log(1.0-((double)((rand()*rand())%1000))/1000)/lambda)))%RADIUS)+1;  
  unsigned int i=i1;
  unsigned int x;
  bool node_exists;
  int size=d;
  int n;
  struct protein* temp;
  struct protein* oldtemp;
  unsigned int proteins[11];
  int MAX_ptrial=5;
  int ptrial=0;
  int j=0; // the length of the path found
  while ((j<d)&&(ptrial<MAX_ptrial)){  // try to find a length d path MAX_ptrial times
	 for(j=0;(j<d);j++){  // add size nodes to subG by following a possibly cyclic path
        node_exists=true;	
		for (int trial=0; ((node_exists)&&(trial<(pin[i].non+1))) ; trial++){  // try a neighbour which does not exist on the path
		  n=(rand() % (pin[i].non));
		  temp=oldtemp=pin[i].nl;
		  for(int k=0 ;((k<n)&&(temp));k++) {oldtemp=temp;temp=temp->next;}
		  if (temp) x=get_protein(pin,temp->name,hts);
		  else if (oldtemp) x=get_protein(pin,oldtemp->name,hts);
		  else return false;
		  if (!(x<hts)) return false;
		  node_exists=false;
		  if (i1==x) node_exists=true;
          else for(int k=0;k<=j;k++) if (proteins[k]==x) node_exists=true;
		  }  // for node_exists
		if (node_exists==false) i=x;
		else break;
		proteins[j]=i; 
	    }  // int j to size
	 proteins[10]=j;
     ptrial++;
     }
  if (ptrial==MAX_ptrial) return false;
  (*ii2)=proteins[d-1];
  // assign a random neighbour of ii2 to jj2		  
  n=(rand() % (pin[*ii2].non));
  temp=oldtemp=pin[*ii2].nl;
  for(int k=0 ;((k<n)&&(temp));k++) {oldtemp=temp;temp=temp->next;}
  if (temp) x=get_protein(pin,temp->name,hts);
  else if (oldtemp) x=get_protein(pin,oldtemp->name,hts);
  else return false;
  if (!(x<hts)) return false;
  else (*jj2)=x;
  //printf(" d: %d",d);
  return true;
  }

bool find_matching_edge(PROINFO* pin, unsigned int i1, unsigned int j1, unsigned int* ii2, unsigned int* jj2){
  // returns an edge (i2-j2) randomly from pin. i2-j2 that has matching labels as i1-i2 .
  unsigned int i2,j2;
  struct protein * temp;
  int start=(unsigned int)((rand()*rand()*13191)%(hts-1));
  for(i2=start+1; ((i2!=start)) ;i2=(i2+1)%hts ) {
    if ((pin[i2].id==NULL)||(pin[i2].non==0)) continue;
	if ((pin[i2].label==pin[i1].label)||(pin[i2].is_my_label(pin[i1].label))) {
      if (i1==i2) continue;
	  for(temp=pin[i2].nl ; ((temp)&&((pin[j2=get_protein(pin,temp->name,hts)].label!=pin[j1].label)||(pin[j2].is_my_label(pin[j1].label)))); temp=temp->next ) ;        
	    if (temp){
           (*ii2)=i2; (*jj2)=j2;   
		   if (j2>=hts) {printf("\nSomething wrong with hash table & get_protein!!"); continue;}
		   return true;
	       }
	   }   
    }
  //printf("\nCould not find a matching edge for %d-%d.",pin[i1].label, pin[j1].label);
  return false;
}

void test9(){
  int sum=0;
  int RADIUS=10;
  double lambda=0.2; 
  double dd;
  for (int i=0 ; i<100 ; i++) {
    dd=((double)((rand()*rand())%1000))/1000;
	printf("\ndd=%f\tF-1=\t%d",dd,(int)(-1*(log(1.0-((double)((rand()*rand())%1000))/1000)/lambda)));
    }
}

bool randomizePIN(PROINFO* pin, int itn){
  // perturb pin itn times using a well defined randomization process
  // that reflects known mechanisms (such as evolution)
  unsigned int i1,i2,x,j1,j2,tt;
  struct protein * temp;
  struct protein * oldtemp;
  int n1,n2,l1,l2,r1;
  float conf1,conf2;
  bool j1_found,i1_found;
  int count=0; int nme_count=0;
  srand( (unsigned)time( NULL ) );
/*
  printf("\n%d : ",pin[21413].id );
  for(temp=oldtemp=pin[21413].nl ; ((temp)); temp=temp->next ) {oldtemp=temp; printf("%d ",get_protein(pin,temp->name,hts));} 
  printf("\n%d : ",pin[19740].id );
  for(temp=oldtemp=pin[19740].nl ; ((temp)); temp=temp->next ) {oldtemp=temp; printf("%d ",get_protein(pin,temp->name,hts));} 
  printf("\n%d : ",pin[21004].id );
  for(temp=oldtemp=pin[21004].nl ; ((temp)); temp=temp->next ) {oldtemp=temp; printf("%d ",get_protein(pin,temp->name,hts));} 
*/
  for(int j=0; j<itn ; j++){
      // do swap between two random edges i1-j1 and i2-j2
	  // First, find these edges
	  // First find edge i1-j1
	  r1=((rand()*rand()*13191)%(nolp-ei))+ei;
	  l1=lpt[soalp[r1]].label1 ;
      l2=lpt[soalp[r1]].label2 ;
	  i1=tt=(unsigned int)((rand()*rand()*13191)%(hts-1));
	  j1_found=false;
	  i1_found=false;

	  do {
	    i1=(i1+1)%hts; 
		if (pin[i1].id){
 		  if ((pin[i1].label==l1)||(pin[i1].is_my_label(l1))) { 
			  i1_found=true;
			  for(temp=pin[i1].nl ; (temp); temp=temp->next ) {
                j1=get_protein(pin,temp->name,hts);
				if ((pin[j1].label==l2)||(pin[j1].is_my_label(l2))) {
					j1_found=true;
					break; }
			    }
		    }  // label

		  } // id
	  } while ((!j1_found)&&(i1!=tt));
	  if (!j1_found)  {  
          // randomly select i1     
	      for(i1=(unsigned int)((rand()*rand()*13191)%hts); ((pin[i1].id==NULL)||(pin[i1].non==0)) ;i1=(unsigned int)((rand()*rand()*13191)%hts));   
          // randomly select node j1 of edge i1-j1
	      n1=((rand()*rand()*13191)%pin[i1].non );
	      temp=oldtemp=pin[i1].nl;
	      for(int k=0 ;((k<n1)&&(temp));k++) {oldtemp=temp;temp=temp->next;}
	      if (temp)
		    x=get_protein(pin,temp->name,hts);
	      else if (oldtemp)
		    x=get_protein(pin,oldtemp->name,hts);
	      else {
			return false;}
	      if (!(x<hts)) return false;      
	      j1=x;
	      }
	  else { count++;
	      }

      // randomly select i2     
	  for(i2=(unsigned int)((rand()*rand()*13191)%hts); ((pin[i2].id==NULL)||(pin[i2].non==0)) ;i2=(unsigned int)((rand()*rand()*13191)%hts));   
	  if (i2>=hts) { printf("\ni2>=hts"); return false;} 
      // randomly select node j1 of edge i1-j1
      /*
	  n1=((rand()*rand()*13191)%pin[i1].non );
	  temp=oldtemp=pin[i1].nl;
	  for(int k=0 ;((k<n1)&&(temp));k++) {oldtemp=temp;temp=temp->next;}
	  if (temp)
		    x=get_protein(pin,temp->name,hts);
	  else if (oldtemp)
		    x=get_protein(pin,oldtemp->name,hts);
	  else {
			return false;}
	  if (!(x<hts)) return false;      
	  j1=x;
	  */
	  // find node j2 of edge i2-j2
      n2=((rand()*rand()*13191)%pin[i2].non );
	  temp=oldtemp=pin[i2].nl;
	  for(int k=0 ;((k<n2)&&(temp));k++) {oldtemp=temp;temp=temp->next;}
	  if (temp)
		    x=get_protein(pin,temp->name,hts);
	  else if (oldtemp)
		    x=get_protein(pin,oldtemp->name,hts);
	  else {
			return false;}
	  if (!(x<hts)) return false;
      j2=x;
	  if (strcmp(RGM,"RGM=F")==0)
  	    if (find_matching_edge(pin,i1,j1,&i2,&j2)) {
		  //printf("\nFound matching edges: (%d-%d)(%d-%d)(%d-%d)",i1,j1,i2,j2,pin[i1].label ,pin[j1].label);
	      }
	    else nme_count++;

	  if (strcmp(RGM,"RGM=L")==0)
	    if (find_neighbour(pin,i1,j1,&i2,&j2))  {
		  //printf("\nFound neighbour edges: (%d-%d)(%d-%d)(%d-%d)",i1,j1,i2,j2,pin[i1].label ,pin[j1].label);
	      }
	    else nme_count++;
	  if ((i1==i2)||(i1==j2)||(j1==i2)||(j1==j2)) { continue;}

	  //printf("\n%s : Swapping %d-%d & %d-%d  ",pin[i1].id , i1,j1,i2,j2);	
  /*
	  printf("\ni1-->%s:%d -->",pin[i1].id , i1 );
      for(temp=oldtemp=pin[i1].nl ; ((temp)); temp=temp->next ) {oldtemp=temp; printf("%d ",get_protein(pin,temp->name,hts));}	      
      printf("\nj1 -->%s:%d -->",pin[j1].id , j1 );
      for(temp=oldtemp=pin[j1].nl ; ((temp)); temp=temp->next ) {oldtemp=temp; printf("%d ",get_protein(pin,temp->name,hts));}	      

      printf("\ni2 -->%s:%d -->",pin[i2].id , i2 );
      for(temp=oldtemp=pin[i2].nl ; ((temp)); temp=temp->next ) {oldtemp=temp; printf("%d ",get_protein(pin,temp->name,hts));}	      

	  printf("\nj2 -->%s:%d -->",pin[j2].id , j2 );
      for(temp=oldtemp=pin[j2].nl ; ((temp)); temp=temp->next ) {oldtemp=temp; printf("%d ",get_protein(pin,temp->name,hts));}	      
*/
	  // remove i1-j1
	  for(oldtemp=NULL , temp=pin[i1].nl ; ((temp)&&(strcmp(temp->name,pin[j1].id))); temp=temp->next ) {oldtemp=temp;}  
	  if (temp==NULL) return false;
	  if (oldtemp==NULL) {
         pin[i1].nl=temp->next ;
	     }
	  else {
         oldtemp->next=temp->next;
	     }
	  pin[i1].non--;
	  conf1=temp->conf ;
	  delete temp;

	  // remove j1-i1
	  for(oldtemp=NULL , temp=pin[j1].nl ; ((temp)&&(strcmp(temp->name,pin[i1].id))); temp=temp->next ) {oldtemp=temp;}  
	  if (temp==NULL) return false;
	  if (oldtemp==NULL) {
         pin[j1].nl=temp->next ;
	     }
	  else {
         oldtemp->next=temp->next;
	     }
	  pin[j1].non--;
	  delete temp;

	  // remove i2-j2
	  for(oldtemp=NULL , temp=pin[i2].nl ; ((temp)&&(strcmp(temp->name,pin[j2].id))); temp=temp->next ) {oldtemp=temp;}  
	  if (temp==NULL) return false;
	  if (oldtemp==NULL) {
         pin[i2].nl=temp->next ;
	     }
	  else {
         oldtemp->next=temp->next;
	     }
	  conf2=temp->conf ;
	  pin[i2].non--;
	  delete temp;
	  // remove j2-i2
	  for(oldtemp=NULL , temp=pin[j2].nl ; ((temp)&&(strcmp(temp->name,pin[i2].id))); temp=temp->next ) {oldtemp=temp;}  
	  if (temp==NULL) return false;
	  if (oldtemp==NULL) {
         pin[j2].nl=temp->next ;
	     }
	  else {
         oldtemp->next=temp->next;
	     }
	  pin[j2].non--;
	  delete temp;

	  // Now , add edges i1-j2 and i2-j1
      insert_assoc(pin,i1,j2,conf1);
      insert_assoc(pin,i2,j1,conf2);
      //printf("inserted..");
/*	
	  printf("\n%s : Swapped %d-%d & %d-%d  ",pin[i1].id , i1,j1,i2,j2);	
      printf("\ni1-->%s:%d -->",pin[i1].id , i1 );
      for(temp=oldtemp=pin[i1].nl ; ((temp)); temp=temp->next ) {oldtemp=temp; printf("%d ",get_protein(pin,temp->name,hts));}	      
      printf("\nj1 -->%s:%d -->",pin[j1].id , j1 );
      for(temp=oldtemp=pin[j1].nl ; ((temp)); temp=temp->next ) {oldtemp=temp; printf("%d ",get_protein(pin,temp->name,hts));}	      

      printf("\ni2 -->%s:%d -->",pin[i2].id , i2 );
      for(temp=oldtemp=pin[i2].nl ; ((temp)); temp=temp->next ) {oldtemp=temp; printf("%d ",get_protein(pin,temp->name,hts));}	      

	  printf("\nj2 -->%s:%d -->",pin[j2].id , j2 );
      for(temp=oldtemp=pin[j2].nl ; ((temp)); temp=temp->next ) {oldtemp=temp; printf("%d ",get_protein(pin,temp->name,hts));}	      
*/
	  }  // for j<itn
   
	  printf("\nrandomization complete with %d iterations %d of which is j1_found, nme_count:%d.",itn,count,nme_count);  
   return true;
  }  

 int create_randomPINS(PROINFO* pht, int count, int rand_count){
    // create count randomized PINS from pht with rand_count shuffles/perturbations
    PROINFO* pin1;
	for (int i=0; i<count ; i++){
       pin1=copyPIN(pht);
	   if (!randomizePIN(pin1,rand_count)) {printf("\nrandomizePIN failed!");break;}
	   randomPIN[i]=pin1;
	   }
	return count;
  }

float compute_Z_score(subG* s){
   int counts[1000];
   int sum=0;
   float stddev=0;
   float avg=0;
   for(int i=0; i<randomPINcount ; i++){
      counts[i]=find_subG3(s,randomPIN[i]);
      //printf("\nEmbeddings counts[%d]=%d",i,counts[i]);
 	  }
   for(int i=0; i<randomPINcount ; i++){
      sum+=counts[i];
 	  }   
   avg=(float)sum/randomPINcount;
   for(int i=0; i<randomPINcount ; i++){
      stddev+=(fabs((avg-counts[i]))*fabs((avg-counts[i])));
 	  }   
   stddev=stddev/randomPINcount;
   stddev=sqrt(stddev);
   s->Zscore = (fabs(((s->freq)-avg)))/stddev;
   s->pValue = pV.computePValue(s->nol , s->nedges() , tlevel , s->freq);
   if (ordering_value==PVALUE)
      s->Zscore = s->pValue;
   printf("\n Z-score: %f , freq:%d",s->Zscore, s->freq);
   //display_subG(s,new unsigned int[11],stdout);
   printf("\nZ-score= %d - %f / %f",s->freq,avg,stddev);
   return (fabs(((s->freq)-avg)))/stddev; 
 }

void print_protein_array(unsigned int p[11]){
    printf("\nProtein Array: ");
	for (unsigned int i=0; i<10 ; i++) printf("%d ",p[i]);
  }

void print_proteins(unsigned int p[11]){
  printf("\nProtein set with %d nodes:",p[10]);
  for (unsigned int i=0; i<10 ; i++) {
	  if (p[i]>0){
	     printf("\np[%d]=%d  pht[%d]=%s,%d ->",i,p[i],p[i],pht[p[i]].id ,pht[p[i]].label);
         for (struct protein* next=pht[p[i]].nl ; (next) ; next=next->next ) printf("%s ",next->name );
	     }
      }
  }

void display_subG(subG* s,unsigned int p[11],FILE* ff){
	fprintf(ff,"\nsubgraph with %d nodes, embeddings:%d , Z-score:%6f ______\n",s->nol,s->freq ,s->Zscore );
	for(int i=0;i<s->nol;i++){
        fprintf(ff,"%d  ",s->labels[i]);
	   }
	for(int i=0;i<s->nol;i++){
		if (p[i]<hts) fprintf(ff,"\n%s.%d -->",pht[p[i]].id , pht[p[i]].label);
		else fprintf(ff,"\n(null).0 -->");
        for(int j=0;j<s->nol;j++){
		   if (s->adjm[i][j]==1) fprintf(ff," 1 ",s->adjm[i][j] );
		   else fprintf(ff," 0 " );
		   }
	   }
	//fflush(ff);
	
	//printf("\n%d  %s", p[9],s->labels[9]);
  }

void read_subG(subG* s,unsigned int p[11],char *sf){
	// Read  subgraphs from file sf, form and return a linked list of subGs
	// Its format is the same as the one output by display_subG()
	// The first line may be anything including blank line
	FILE* ff=fopen(sf,"r");
	char temp[20];
	char temp1[20];
    char temp2[20];
    char temp3[20];
	char buffer[200];
	subG* olds=NULL;
    char* cp;
    int count=0;  // subg count

	while (!feof(ff)) {

	fgets( buffer, 200, ff ); 
	if (strstr(buffer,"subgraph")){
	sscanf( buffer, "%s %s %d",temp,temp,&(s->nol));
	if (cp=strstr(buffer,"embeddings:"))
      sscanf( cp+11, "%d",&(s->freq));
	if (cp=strstr(buffer,"embeddings:"))
      sscanf( cp+11, "%d",&(s->freq));
	if (cp=strstr(buffer,"Z-score:"))
      sscanf( cp+8, "%f",&(s->Zscore));
	//fscanf(ff,"%s %s %d",temp,temp,&(s->nol));
	for(int i=0;i<s->nol;i++){
        fscanf(ff,"%s",temp3 );
	    s->labels[i]=atoi(temp3);
	    }
	fprintf(stdout,"\nsubgraph with %d nodes, embeddings:%d , Z-score:%f ______\n",s->nol,s->freq ,s->Zscore );
    //for(int i=0;i<s->nol;i++){
    //    fprintf(stdout,"%d  ",s->labels[i] );
    //    }
	//fprintf(stdout,"\n" );
	p[10]=s->nol;
	for(int i=0;i<s->nol;i++){
		fscanf(ff,"%s %s",temp1,temp2);
		(*(strchr(temp1,'.')))='\0';
		p[i]=get_protein(pht,temp1,hts);
		//if (p[i]>=hts) fprintf(stdout,"\n%s . %d %s\t",temp1, p[i], temp2 );
		//fgets( buffer, 200, ff ); 
        for(int j=0;j<s->nol;j++){
           fscanf(ff,"%d",&(s->adjm[i][j]));
		   //if (s->adjm[i][j]==1) fprintf(stdout," 1 ");
    	   //else if (s->adjm[i][j]==0) fprintf(stdout," 0 " );
		   }
        //fprintf(stdout,"\n" );
	   }
	s->next=new subG();
	s->next->next=NULL;
	olds=s;
	s=s->next;
    count++;  // subG count
	}  // if buffer


	}  // while
    if (olds) delete olds->next;
	olds->next=NULL;
    fclose(ff);
  }

bool neighbours(unsigned int p1, unsigned int p2){ // are these neighbours
		struct protein* temp = NULL;
		for (temp=pht[p1].nl ; ((temp)&&(strcmp(pht[p2].id,temp->name ))) ; temp=temp->next);
		if (temp) return true;
		for (temp=pht[p2].nl ; ((temp)&&(strcmp(pht[p1].id,temp->name ))) ; temp=temp->next);
		if (temp) return true;
		return false;
	  }

int verify_sub(unsigned int map[11],subG* s, bool show){
	int result=1;
	for (int i=0;i<s->nol ; i++) 
		    if (!(map[i]) ) {
			    if (show) printf("protein in map[%d]=0",i);  
				return 0;}
			else if ((pht[map[i]].label!=s->labels[i])&&(!pht[map[i]].is_my_label(s->labels[i])))  {
				if (show) printf("label mismatch: label of map[%d]!=%d",map[i],s->labels[i]);  
				return 0;}
	if (show) printf("\nsubgraph with %d nodes:______\n",s->nol);
	for(int i=0;i<s->nol;i++){
        if (show) printf("%d  ",s->labels[i] );
	   }
	for(int i=0;i<s->nol;i++){
		if (show) printf("\n%s : %d\t-->",pht[map[i]].id , map[i]);
        for(int j=0;j<s->nol;j++){
			   if (neighbours(map[i],map[j]))
				   if (s->adjm[i][j]==1){
					   if (show) printf(" 1 " );}
				   else {
					   if (show) printf(" X " );}
			   else
				   if (s->adjm[i][j]==1) {
			           if (show) printf(" - " );
					   result=0;
				       }
				   else {
					   if (show) printf(" 0 " );}

		   }
	   }        
    return result;   
   }

char* graph_text_drawing(subG* s,unsigned int position[11]){
  int isec;  // intersetion point,count
  char *td=(char*)malloc(101*s->nol);  // text drawing box
  memset(td, ' ', 101*s->nol );
  for (int i=0 ; i<s->nol ; i++) {
	  strncpy(td+(i*101)+(i*10),pht[position[i]].id,strlen(pht[position[i]].id));  
      strncpy(td+((i+1)*101)-1,"\n",1);
      }
  for(int i=0;i<s->nol;i++){
      for(int j=i+1;j<s->nol;j++){
		if (s->adjm[i][j]==1){
          isec=(i*101)+((j*10)+3);
          td[isec]='|';  // boxing char 0xBF in MSDOS code page is a better fit here 
		  for(int k=isec-1; k>(i*101)+(i*10) ; k--)
			  if (td[k]==' ') td[k]='-';
		  for(int k=isec+101, c=0; c<j-i ; k+=101,c++)
			  if (td[k]==' ') td[k]='|';		  
		
		  }
		}
	  }

  td[(101*s->nol)-1]='\0';
  return td;

}

int append_to_file(char* rfname ,subG* s , unsigned int pr[11], unsigned int positions[1000][11]){
    // open rf file to append
	//printf("opening result file %s ...", rfname);
	//if (compute_Z_score(s)!=compute_Z_score(s)) printf("\nThey are not equal anam.");
	char* gtd;
	if (s->Zscore==0) compute_Z_score(s);
	FILE* rf=fopen(rfname, "a");
   	FILE* pf=fopen("positions.txt", "a");
	display_subG(s,pr,rf);
	if (positions){
 	  fprintf(rf,"\n [IDs]:\n");
	  for (int i=0;i<s->freq ; i++) {
		for(int j=0 ; j<s->nol ; j++)  {
		  fprintf(rf," %s,\t",pht[positions[i][j]].id );
		  fprintf(pf,"\n%s ",pht[positions[i][j]].id );
		  }
	    fprintf(rf,"\n");
	    }
	  fprintf(rf,"\n%s",gtd=graph_text_drawing(s,positions[0])); free(gtd);
	  }
	fprintf(rf,"\n------------");
	fclose(rf);
	fclose(pf);

    // write the rank list to file
  if ((rand()%10000)>9990){  // with 0.001 probability
    FILE* plf=fopen(rankfile, "a"); fprintf(plf,"\n%s\n-------------",command);
    for (subG* temp=toplist->list ; (temp) ; temp=temp->next) {
	  display_subG(temp,new unsigned int[11],plf);
	  fprintf(plf,"\n------------");
      }
    fclose(plf);
    } // if rand()
	return 1;
   }



bool read_goterms(char *obofile, unsigned int noterms){ 
	// reads all obo terms from obofile and places them
	// into the array of term pointers gtt.
	// m=number of go terms read
	printf("\nReading OBO file %s with param %s...", obofile,go_subset);
	char temp1[100];
	char temp2[100];
	int n=0; int m=0;
	char buffer[2000];
    char* name;
	char namesp;
	bool last_stanza_is_term=false;
	unsigned int id;
	unsigned int isaid;
	unsigned int maxisa_id;
    unsigned int maxisa=0;
    unsigned int partofid;
	unsigned int maxpartof_id;
    unsigned int maxpartof=0;
    unsigned int altid=0;

	for(int i=0;i<100000;i++) gtt[i]=NULL;
	FILE* obof=fopen(obofile, "r");
	if (obof==NULL)
		printf("\nCould not read OBO file %s with param %s...", obofile,go_subset);

    while (!feof(obof)) {
		if (fgets( buffer, 2000, obof )==NULL) {printf("\nfgets error!,%d",id); continue;}
		//printf("\nAfter fgets! %s",buffer); 
		if (strstr(buffer,"[Term]")==buffer) {
			last_stanza_is_term=true;
			//if (altid>0) gtt[altid]=gtt[id];
            id=0; name=NULL; namesp='\0';
			altid=0;
		    }
		else if (strstr(buffer,"[")==buffer)
            last_stanza_is_term=false;
		if (last_stanza_is_term){
		    if (strstr(buffer,"id:")==buffer) {  // if ID tag
			  if (id) {
				printf("\nUnexpected ID tag at %dth line:%s",n,buffer);
				return false; 
			    }
			  sscanf( buffer, "%10s%3s%d ", temp1, temp2, &id );
              if ((id<0)||(id>noterms)) return false;
			  gtt[id]=new goterm;   
			  m++;
			  }			
			else if (strstr(buffer,"name:")==buffer) {  // if name tag       
			  if (id){  
			     gtt[id]->name =(char*)malloc(100);
				 //printf("\n Inside id gtt[id]->name: %s",gtt[id]->name);
			     }
			  else { 
				  return false; }
			  if (strlen(buffer)<8) printf("\nname of term missing.."); 
			  else {
				  strncpy(gtt[id]->name,buffer+6,99); //gtt[id]->name[99]='\0';
				  }
			  name=gtt[id]->name;
			  //printf("\n Name At %dth line:%s",n,name);
			  }
			else if (strstr(buffer,"namespace:")==buffer) {  // if namespace tag
			  //printf("\n Name At %dth line:%s",n,buffer);              
			  if (!id) return false;
			  if (strlen(buffer)<13) printf("\nnamespace of term missing.."); 
			  else {
                  if (strstr(buffer,"function")) gtt[id]->namesp='F';
				  else if (strstr(buffer,"process")) gtt[id]->namesp='P';
				  else if (strstr(buffer,"cellular")) gtt[id]->namesp='C';
				  else gtt[id]->namesp='\0';		      
			      }
			  namesp=gtt[id]->namesp;
			  }
			else if (strstr(buffer,"def:")==buffer) {  // if definition tag
			  //printf("\n Name At %dth line:%s",n,buffer);              
              if (id)
			    gtt[id]->def =(char*)malloc(300);
			  else 
				  return false;
			  if (strlen(buffer)<8) printf("\ndefinition of term missing.."); 
			  else {
				  strncpy(gtt[id]->def,buffer+6,299);
				  gtt[id]->def[299]='\0';
			      }
			  }
            // subset: 
			else if ((strstr(buffer,"subset:")==buffer)&&(strcmp(go_subset,""))) {  // if subset tag such as goslim_yeast
                //printf("\n Subset tag at %dth line, with GO ID[%d]:%s",n,id,buffer);
				if (strstr(buffer,go_subset)) {
                  printf("Subset tag at %dth line, with GO ID[%d]:%s",n,id,buffer);
				  if (id)
			        gtt[id]->subset =new char[20];
			      else 
				      return false;
			      if (strlen(buffer)<10) printf("\nsubset of term missing.."); 
			      else {
				      strncpy(gtt[id]->subset,go_subset,19);
				      gtt[id]->subset[19]='\0';
			          }
			      }  // if go_subset
			  }
            
		    if (strstr(buffer,"is_a:")==buffer) {  // if is_a tag
			  sscanf( buffer, "%10s%3s%d ", temp1, temp2, &isaid );
              if ((isaid<0)||(isaid>noterms)) return false;
			  if (gtt[id]->is_a[10]<=10) {
				  gtt[id]->is_a[gtt[id]->is_a[10]]=isaid;
			      gtt[id]->is_a[10]++;}
			  else return false;  
			  if (gtt[id]->is_a[10]>maxisa) {maxisa=gtt[id]->is_a[10]; maxisa_id=id;}
			  }
		    if (strstr(buffer,"relationship: part_of")==buffer) {  // if part_of tag
			  sscanf( buffer, "%s%s%3s%d ", temp1, temp1, temp2, &partofid );
              if ((partofid<0)||(partofid>noterms)) return false;
			  if (gtt[id]->part_of[10]<=10) {
				  gtt[id]->part_of[gtt[id]->part_of[10]]=partofid;
			      gtt[id]->part_of[10]++;}
			  else return false;  
			  if (gtt[id]->part_of[10]>maxpartof) {maxpartof=gtt[id]->part_of[10]; maxpartof_id=id;}
			  }
			if (strstr(buffer,"alt_id:")==buffer) {  // if alternate id tag
			  sscanf( buffer, "%10s%3s%d ", temp1, temp2, &altid );
              if ((altid<0)||(altid>noterms)) altid=0;
			  else gtt[altid]=gtt[id];
			  }

		    }
        n++;
	    }
      //printf("\nlast: %s  %d name:%s",temp1,id,name);
	  printf("\nmax_isa: %d   max_isa_id:  %d ",maxisa,maxisa_id);
	  printf("\nmax_partof: %d   max_partof_id:  %d ",maxpartof,maxpartof_id);
	  printf("\nInserted GO Terms:%d   Lines Read:%d\n",m,n);
	  fclose(obof);   
      return true;
   }


unsigned int select_a_parent_term(unsigned int x){ // of x
   // favor 'F'=molecular Function
   unsigned int spare=0;
   unsigned int i;
   for(i=0; ((i<gtt[x]->is_a[10])&&(gtt[x]->is_a[i])) ; i++ )
	   if (gtt[gtt[x]->is_a[i]]->namesp == 'F') return (gtt[x]->is_a[i]);  // ilk buldugun F'i dönder
	   else if (gtt[gtt[x]->is_a[i]]->namesp == 'P') spare=gtt[x]->is_a[i];
	   else if (!spare) spare=gtt[x]->is_a[i] ;  // 'C'
   if (spare) return spare;  
   else if (gtt[x]->part_of[10]){
      for(i=0; ((i<gtt[x]->part_of[10])&&(gtt[x]->part_of[i])) ; i++ )
	    if (gtt[gtt[x]->part_of[i]]->namesp == 'F') return (gtt[x]->part_of[i]);
	    else if (gtt[gtt[x]->part_of[i]]->namesp == 'P') spare=gtt[x]->part_of[i];
	    else if (!spare) spare=gtt[x]->part_of[i] ;  // 'C'
      if (spare) return spare;
      }
   return 0;   // this is a root term      
}

void set_gtt_levels(unsigned int root){ 
	// sets the level_plus1 field of each 
 int hor_count=0;
 gtt[root]->level_plus1=1;
 bool deeper=true;
 for(int levels=1; (deeper) ; levels++){
    deeper=false;
	hor_count=0;
	for (unsigned int x=0; (x<100000) ; x++) {
		for(unsigned int i=0; ((gtt[x])&&(i<gtt[x]->is_a[10])&&(gtt[x]->is_a[i])) ; i++ ){
		   if (gtt[gtt[x]->is_a[i]]->level_plus1==levels) {
			   gtt[x]->level_plus1=levels+1;
			   deeper=true;
			   }
		   } // i
        if ((gtt[x])&&(gtt[x]->level_plus1 == (levels+1)))  hor_count++;
		//if ((gtt[x])&&(gtt[x]->level_plus1 >= 13)) printf("\nThe level_plus1 of GO term %d is %d.",x,gtt[x]->level_plus1);
	    }  // x
    printf("\nHorizontal count for level %d is %d.",levels-1,hor_count);
    }  // levels
}
void update_labels2(int level){  // set the level of all labels to level
	// it replaces each label in pht with its (grand)parent whose distance to the root GO term is level
    // it first finds path to the root and then, sets the pht[i].label to  the levelth term on this path
    // level zero is the root level
	int n,j;
	unsigned int p;
	unsigned int path[100];
	unsigned int anc[5]; for (int k=0;(k<5);k++) anc[k]=0; // ancestors
	
	for(unsigned int i=0;i<hts;i++){
		if (pht[i].id){
			//printf("%s->",pht[i].label);
			n=0;
			for (p=(unsigned int)pht[i].label; ((p!=0)&&(n<100)) ; p=select_a_parent_term(p)) {path[n]=p; n++;}
			if (n<100)
			  if (level<n) 
				  pht[i].label=path[n-level-1];
			  else ; // don't change the label because it's depth is less than level
			else 
              printf("\nGTT has maxdepth>100. Can not set labels. ");

			for (j=0 ; ((j<5)&&(pht[i].otherlabels[j])) ; j++ ) {
			  n=0;
			  for(int k=0 ; (k<100) ; k++) path[k]=0;
			  for (p=(unsigned int)pht[i].otherlabels[j]; ((p!=0)&&(n<100)) ; p=select_a_parent_term(p)) {path[n]=p; n++;}
			  if (n<100)
			    if (level<n) 
		  		   pht[i].otherlabels[j]=path[n-level-1];
			  } // j
			/*
			for(p=(unsigned int)pht[i].label; ((p!=0)&&(n<up)) ; p=gtt[p]->is_a[0]) { n++;}
			if (p!=0)
              pht[i].label=(int)p;
            */
			/*
			if (p==0)
				itoa(atoi(pht[i].label),pht[i].label,10); 
			else
				itoa(p,pht[i].label,10); 
            */ 
			//printf("%s ",pht[i].label);
		    }
	   }

   }


unsigned int* get_ancestors_of_GO_term(unsigned int gterm,unsigned int* anc,int level) {
unsigned int ancgrid[50][101];  // ancgrid[2][] will contain all ancestors of gterm at level 2 etc.
for(unsigned int i=0;(i<50) ; i++)	  for (unsigned int j=0 ; (j<101) ; j++) ancgrid[i][j]=0;
unsigned int m;
ancgrid[gtt[gterm]->level_plus1-1][0]=gterm; ancgrid[gtt[gterm]->level_plus1-1][100]=1; 
bool grid_changed=true;
while (grid_changed) {
  grid_changed=false;
  for(unsigned int i=0;(i<50) ; i++)
	  for (unsigned int j=0 ; ((j<ancgrid[i][100])&&(j<100)) ; j++)
		  if (ancgrid[i][j]) {
			  for(unsigned int k=0; (k<gtt[ancgrid[i][j]]->is_a[10]) ; k++)  { // for each ancestor of ancgrid[i][j]
				  for(m=0; (m<ancgrid[gtt[(gtt[ancgrid[i][j]]->is_a[k])]->level_plus1-1][100]) ; m++ ) // check if in ancgrid
                       if (ancgrid[gtt[(gtt[ancgrid[i][j]]->is_a[k])]->level_plus1-1][m]==gtt[ancgrid[i][j]]->is_a[k]) break;
			      if ((m==ancgrid[gtt[(gtt[ancgrid[i][j]]->is_a[k])]->level_plus1-1][100])&&(m<100)) { // if not in ancgrid
                     ancgrid[gtt[(gtt[ancgrid[i][j]]->is_a[k])]->level_plus1-1][m]=gtt[ancgrid[i][j]]->is_a[k];
                     ancgrid[gtt[(gtt[ancgrid[i][j]]->is_a[k])]->level_plus1-1][100]++;
                     grid_changed=true;
			         }
			      }
		     }  // if ancgrid[][]


  }  // while 
if (ancgrid[level][100]==0) return NULL;
for(unsigned int j=0 ; j<ancgrid[level][100] ; j++) 
  anc[j]=ancgrid[level][j] ;
anc[100]=ancgrid[level][100];

if (strcmp(go_subset,"")){  //  if we want to obtain terms from a GO subset then, find them in ancgrid[]  
  unsigned int k=0;
  for(unsigned int i=0;(i<50) ; i++)	  
	for (unsigned int j=0 ; ((j<ancgrid[i][100])&&(j<100)) ; j++) 
		if ((gtt[ancgrid[i][j]]->subset)&&(strcmp(gtt[ancgrid[i][j]]->subset,go_subset)==0)) {
			unsigned int q;
			for(q=0; (q<k) ; q++) if (anc[q]==ancgrid[i][j]) break;
			if (q==k) {
			  anc[k]=ancgrid[i][j]; 
			  k++;
			  }
		    }
  if (k==0) {
	  printf("\nThe go term %d has no GO subset term as ancestor.",gterm);
      anc[0]=gterm; anc[100]=1;
      }
  else if (k==1) {anc[0]=gterm;anc[100]=k;}  // the only subset term is the root 3674, 
                         // then keep the gterm, there are 200 such proteins for yeast
  else {  // (k>1) then remove the root go terms from anc
	  for(unsigned int t=0; (t<k); t++){
        for(unsigned int tt=0; (tt<ancgrid[0][100]) ; tt++)
			if (anc[t]==ancgrid[0][tt]) {  // remove anc[t]
               for(unsigned int ttt=t+1; (ttt<k) ; ttt++)   anc[ttt-1]=anc[ttt];
			   k--;
			   }
	    }
	  if (k<1)  {anc[0]=gterm;anc[100]=1; }
	  else { anc[100]=k;}
      }  // else

  //if (k>1) printf("\n protein with go_slim  %d from %d",anc[1],gterm);
  }
/*
for(unsigned int i=0;(i<50) ; i++)	  {
	printf("\n%d : ",i);
	for (unsigned int j=0 ; ((j<ancgrid[i][100])&&(j<100)) ; j++) printf(" %d",ancgrid[i][j]);
  }
*/
return anc;
}

void update_labels3(int level){  // set the level of all labels to level
	// it replaces each label in pht with its (grand)parent whose distance to the root GO term is level
    // level zero is the root level
	int j,k,n,p;
	unsigned int *res=new unsigned int[101];
    unsigned int temp[5];
	for(unsigned int i=0;i<hts;i++){
		if (pht[i].id){
			memcpy(temp,pht[i].otherlabels,5*sizeof(int));
			// first get the ancestors of label and set it
			if (get_ancestors_of_GO_term((unsigned int)pht[i].label,res,level)) {
                pht[i].label=(int)res[0];
				for(j=1;((j<6));j++) {
					if (j<((int)res[100])) pht[i].otherlabels[j-1]=res[j];
					else pht[i].otherlabels[j-1]=0;
				    }
			   }
			// then get the ancestors of otherlabels (if exist) and set them
			for(j=0; (j<5) ; j++) if (!pht[i].otherlabels[j]) break;
			for(k=0; ((k<5)&&(j<5)) ; k++){
			   if (!temp[k]) continue;
			   if (get_ancestors_of_GO_term(temp[k],res,level)) {
				  for (n=0 ; ((n<5)&&(n<(int)res[100])) ; n++){  // for each ancestor term in res[]
                     if (res[n]==pht[i].label) continue;
					 for(p=0 ; (p<5) ; p++) if (pht[i].otherlabels[p]==res[n]) break;
					 if (p<5) continue; // same label
                     for(p=0 ; (p<5) ; p++) if (!pht[i].otherlabels[p]) break;                   
                     if (p>=5) break;
                     pht[i].otherlabels[p]=res[n];
				     }
			      }
			   for(j=0; (j<5) ; j++) if (!pht[i].otherlabels[j]) break;
			   }  // k
			
		    }  // if id
	   }  // for i
   printf("\nAll Labels updated !");
   }

void test8(){
unsigned int *res=new unsigned int[101];
/*
//printf("\n31263---is a:%d ,%d",gtt[31263]->is_a[0], gtt[31263]->is_a[1]);
if (get_ancestors_of_GO_term(31459,res,10)) 
  for (unsigned int i=0; i<res[100] ; i++) 
	  printf("GO:%d",res[i]);
*/


}

bool pht_vs_gtt(){
	for(unsigned int i=0;i<hts;i++){
       //if (i>35800) printf(" -%d -%d |",i,pht[i].label);
  	   if (pht[i].id)
       if ((gtt[pht[i].label])==NULL ) 
		   printf("\nFirst label %d of protein at index %d is not found in obo file.",pht[i].label,i);
	   for (int j=0; ((pht[i].otherlabels[j])&&(j<5)) ; j++) 
		   if ((gtt[pht[i].otherlabels[j]])==NULL ) 
		     printf("\nlabel %d of protein at index %d is not found in obo file.",pht[i].otherlabels[j],i);
	   }
    // calculate the average distance of a protein's label(go term) to the root term.
	int n; int total=0; int pr=0;
	//printf("\n	gtt[4]->is_a[0]=%d",gtt[4]->is_a[0]);
	//printf("gtt[8150]=%d gtt[4]=%d ",gtt[8150] , gtt[4] );
	for(unsigned int i=0;i<hts;i++){
		if (pht[i].id){
			pr++;
			n=0;
			//printf("\n%d : ",pht[i].label);
			for(unsigned int p=gtt[pht[i].label]->is_a[0]; (p!=0) ; p=gtt[p]->is_a[0]) { 
			   n++;
			   }
            total += n;     
			//printf("%d",n);
			}
	   }
	float average=(float)((float)total/(float)pr);
	printf("\n The average level of a protein in PIN is: %f",average);
    
	// calculate average number of labels per protein
	int count=0;
    int sum=0;
    int j;
    for(unsigned int i=0;i<hts;i++){
	  if (pht[i].id) {
		  //printf("\nProtein:%s  Labels:%d",pht[i].id,pht[i].label);  // label yazdir
	    count++;
		for(j=0; ((j<5)&&(pht[i].otherlabels[j])) ; j++) {
		  if (pht[i].label==pht[i].otherlabels[j]) printf("\nLAbel % d repeated.",pht[i].label);
		  //printf(" %d",pht[i].otherlabels[j]);  // ek label yazdir
		  }
	    sum=sum+1+j;
	    }
      }
    printf("\nProtein Count:%d , Label sum:%d , LAbel per protein:%f",count, sum, (float)((float)sum/(float)count));

	return true;
  }

void statistics(){
// about pht

	// First see if non and nl values are consistent
for(unsigned int i=0;i<hts;i++){
  if (pht[i].id){
    if (((pht[i].non ==0)&&(pht[i].nl!=NULL))||((pht[i].non)&&(pht[i].nl==NULL))) 
		 printf("%d th protein has inconsistent non and nl values",i);
    }
  }
}

void end_pgm(char** argv){
  FILE* f;
  f=fopen(resfile,"a");
  fprintf(f,"\nThe folowing command executed successfully\n%s",command);
  fclose(f);
  f=fopen(maxresult_file,"a");
  fprintf(f,"\nThe folowing command executed successfully\n%s",command);
  fclose(f);
  printf("\nThe folowing command executed successfully\n%s",command);
  printf("\nEnd of GOINPIN..");
}

subG* random_subgraph_generator(int size, int* v,unsigned int* p,unsigned int proteins[11],int* nnodes, int* nedges){
    // extracts a size-node small connected graph from pht , 
	// p and v are the matched nodes in subgraph and pht respectively     
	// proteins[] is the protein indexes of the subgraph following the same node order of subG
	// all parameters except size are outputs. 
	 float conp=(float)0.5;  // probability of a node making a connection to another node in the subgraph
	 bool node_exists;
	 int n,z;
	 unsigned int i,pt,x;
	 n=((rand()*rand()*13191)%hts );
	 for(i=(unsigned int)((rand()*rand()*13191)%hts); (!pht[i].id) ;i=(unsigned int)((rand()*rand()*13191)%hts));   
	 if (i>=hts) return NULL;
	 (*p)=i;
	 (*v)=0;  // the first one in s->labels corresponds to the starting protein
	 subG *s=new subG;
	 struct protein* temp;
	 struct protein* oldtemp;
	 //unsigned int proteins[11];
	 if (pht[i].non==0) return NULL;
	 int j=0;
	 for(j=0;(j<size);j++){  // add size nodes to subG by following a possibly cyclic path
		s->labels[j]=pht[i].label;
		for (z=0; ((pht[i].otherlabels[z])&&(z<5)) ; z++);
		if (z) if (((rand())%(z+2))<z) if (z==1) s->labels[j]=pht[i].otherlabels[0];
		                               else      s->labels[j]=pht[i].otherlabels[rand()%z];
		proteins[j]=i;
        node_exists=true;	
		for (int trial=0; ((node_exists)&&(trial<size)) ; trial++){  // try a neighbour which does not exist on the path
		  n=(rand() % (pht[i].non));
		  temp=oldtemp=pht[i].nl;
		  for(int k=0 ;((k<n)&&(temp));k++) {oldtemp=temp;temp=temp->next;}
		  if (temp)
		    x=get_protein(pht,temp->name,hts);
		  else if (oldtemp)
		    x=get_protein(pht,oldtemp->name,hts);
		  else {
			return NULL;}
		  if (!(x<hts)) return NULL;
		  node_exists=false;
          for(int k=0;k<=j;k++) if (proteins[k]==x) node_exists=true;
		  }  // for node_exists
		if (node_exists==false)		i=x;
		else break;
 	    s->adjm[j][j+1]=1;
            s->adjm[j+1][j]=1;
	    }  // int j to size
	 size=j;
	 proteins[10]=size;
	 (*nedges)=size-1;
/*
	 for(int k=0 ;k<size;k++)  { // remove duplicates
		 //printf("[-3-]-id:%s -label:%s  ",pht[proteins[k]].id,s->labels[k]);
		 for(int m=k+1 ;m<size;m++) 
			 if (proteins[k]==proteins[m]){
				  //printf("\nReduced size to %d for protein pht[%d]=%s ",size-1,proteins[m],pht[proteins[m]].id );
				  for(int t=m; t<size ; t++) {
					  proteins[t]=proteins[t+1]; 
					  s->labels[t]=s->labels[t+1];}
                  size--; m--;
			      }
			 else if (m==(k+1)){  // connect the nodes on the path
				   s->adjm[k][m]=1;
                   s->adjm[m][k]=1;                  
			      }
	     }
*/
    if (size<=1) return NULL;
    for(int k=0 ;k<size;k++)  { // add edges to the path with probability conp
		for (struct protein*	tmp=pht[proteins[k]].nl; (tmp!=NULL) ; tmp=tmp->next){
            pt=get_protein(pht,tmp->name,hts);
			//if (pt==hts) return NULL;
			for (int d=0; d<size ; d++) 
			   if ((k!=d)&&(pt==proteins[d])&&(rand()<(RAND_MAX*conp))){
				   if (s->adjm[k][d] !=1) {
				     s->adjm[k][d]=1;
                     s->adjm[d][k]=1;
				     (*nedges)++;
				     }
				   //printf("k=%d,d=%d   ",k,d);
			      }
		   }
	   }
	int kk; int min=0;
	for (kk=0; kk<size ; kk++) if (pht[proteins[kk]].non<pht[proteins[min]].non) min=kk;
    (*p)=proteins[min];
	(*v)=min;
	s->nol=size;
	proteins[10]=size;
	(*nnodes)=size;
	//printf("\nEnd of random_subgraph_generator()");
	return s;
}

void remove_duplicate_patterns(char* ff,bool bonfer){
	   //read subGs from file argv[4], remove duplicates and write to file "rdp_xxx"
	   // also bonberroni corrects the Z-scores if bonfer==true
       unsigned int ppr[11]; for(int k=0; k<11;k++)  ppr[k]=0;   
       subG* s=new subG();
       read_subG(s,ppr,ff);  // reads all subgraphs in file ff and returns them in a linked list pointed to by s
       subG* olds=NULL;
	   int nn=0; int nnn=0;
	   subG* oldtempp=NULL;
	   subG* temp;
	   for (temp=s; (temp) ; temp=temp->next) {
          //if (olds) delete olds;
		  for  (subG* tempp=temp->next; (tempp) ; tempp=tempp->next){
			  if (tempp->is_subgraph(temp)){
				  nn++; //display_subG(temp,ppr,stdout); display_subG(tempp,ppr,stdout); printf("\n");
                  temp->status=1;
			      }
			  else if (temp->is_subgraph(tempp)) {
				  nn++;//display_subG(tempp,ppr,stdout);display_subG(temp,ppr,stdout);printf("\n");
			      tempp->status=1;
			      } 
		     }
		  olds=temp;
		  if (temp->status!=1) nnn++;
	      }
	   printf("\n%d subGs equivalent from file %s.",nn,ff);
	   olds=NULL;nn=0;
	   subG* tobedeleted=NULL;
	   int number_of_patterns=nnn;
	   nnn=0;
	   char rdpf[100]="rdp_";rdpf[94]='\0';
       FILE* plf=fopen(strncat(rdpf,ff,90), "w"); fprintf(plf,"\n%s\n-------------",command);

	   for (temp=s; (temp) ; temp=temp->next) {
          if (olds) delete(olds);
		  nnn++;
		  if (temp->status==1) // remove from list
              nn++;
		  else {
             if (bonfer)
				 temp->Zscore=(float)bonferroni_correct_Z_score((double)temp->Zscore ,number_of_patterns);
	         display_subG(temp,ppr,plf);
		     fprintf(plf,"\n------------");
		     }
		  olds=temp; 
		  /*
		  if (tobedeleted){
             delete(tobedeleted);
			 tobedeleted=NULL;
		     }
		  if (temp->status==1){ // remove from list
              nn++;
			  if (olds){
                 olds->next=temp->next;
			     }
			  else {
                 s=s->next; 
			     }
			  tobedeleted=temp;
		      }
		  else {
	         //display_subG(temp,new unsigned int[11],plf);
	         //fprintf(plf,"\n------------");
             // fflush(plf);
		     }
		  */
	      }
	   fclose(plf);
	   printf("\nRemoved %d duplicate patterns out of %d",nn,nnn);
}

int find_patterns_in_both(char* ff,char* gg){
	   //read subGs from file ff and gg, compare all, and report intersection in file "x_ff_gg"
       unsigned int ppr[11]; for(int k=0; k<11;k++)  ppr[k]=0;   
       subG* s=new subG();
       read_subG(s,ppr,ff);  // reads all subgraphs in file ff and returns them in a linked list pointed to by s
	   for(int k=0; k<11;k++)  ppr[k]=0;   
       subG* s2=new subG();
       read_subG(s2,ppr,gg);
       subG* olds=NULL;
	   int nn=0; int nnn=0;
	   subG* oldtempp=NULL;
	   subG* temp;
	   char ff_gg[210]="x_";ff_gg[202]='\0';
	   strncat(strncat(ff_gg,ff,100),gg,100);
       FILE* plf=fopen(ff_gg, "a"); fprintf(plf,"\n%s\n-------------",command);
	   for (temp=s; (temp) ; temp=temp->next) {
          //if (olds) delete olds;
		  for  (subG* tempp=s2; (tempp) ; tempp=tempp->next){
			  if ((tempp->is_subgraph(temp))){
				  nn++; display_subG(temp,ppr,plf); //display_subG(tempp,ppr,plf); printf("\n");
                  temp->status=2; tempp->status=2;
				  break;
			      }
		      }
		  olds=temp;
		  nnn++; 
	      }
       fprintf(plf,"\nFile %s contains %d patterns. \n%d of them are in file %s . ",ff,nnn,nn,gg);
       fprintf(plf,"\n*****************************************************");
	   nn=0;nnn=0;   
	   for (temp=s2; (temp) ; temp=temp->next) {
		  for  (subG* tempp=s; (tempp) ; tempp=tempp->next){
			  if ((tempp->is_subgraph(temp))){  //(tempp->nol==temp->nol)&&(tempp->nedges()==temp->nedges())&&
				  nn++; display_subG(temp,ppr,plf); //display_subG(tempp,ppr,plf); printf("\n");
                  temp->status=2; tempp->status=2;
				  break;
			      }
		      }
		  olds=temp;
		  nnn++;
	      }
	   
	   fprintf(plf,"File %s contains %d patterns.\n %d of them are in file %s . ",gg,nnn,nn,ff);
	   fclose(plf);
       return nn;
}

void bonferroni_correct_file(char *fff){
	FILE* ff=fopen(fff,"r");
	char ff_gg[210]="bonf_";ff_gg[202]='\0';
	strncat(ff_gg,fff,100);
	float z;
	double zz;
	char buffer[200];
	char* cp;
	int nn=0;  //
	while (!feof(ff)) {
	  fgets( buffer, 200, ff ); 
	  if (cp=strstr(buffer,"Z-score:")) nn++;
	  }
	fclose(ff);
	ff=fopen(fff,"r");
    FILE* gg=fopen(ff_gg, "w");
	while (!feof(ff)) {
	  fgets( buffer, 200, ff ); 
	  if (cp=strstr(buffer,"Z-score:")){
        sscanf( cp+8, "%f",&z);
        zz=(double)bonferroni_correct_Z_score((double)z,nn);
	    //fscanf(ff,"%s %s %d",temp,temp,&(s->nol));
	    //fseek(ff,fp+(cp-buffer)+8,SEEK_SET);
		sprintf(cp+8,"%6f ______\n",zz);
	    }
      fprintf(gg,"%s",buffer );
	  }
	fclose(ff);
	fclose(gg);
	printf("\nBonf corrected based on %d experiments",nn);
}

void find_embeddings_of_non_labelled_patterns(char* ff){
	// Finds all of the embeddings(up to 1000) of the non-labelled patterns in ff.
	// Writes them to the result file. Pgm should end after this procedure because all lables are set.
    // 1st step: set all labels to a unique value e.g. 3333
	int nn=0;
	for(unsigned int i=0;i<hts;i++){
      if (pht[i].id)
		if ( pht[i].label>0 ){
          pht[i].label=3333;
          for (int j=0; ((pht[i].otherlabels[j])&&(j<5)) ; j++) 
            pht[i].otherlabels[j]=0;
		  nn++;
		  }
	  }  
     printf("\n Erased the labels of %d proteins in pht.",nn);
     // 2nd step: Read all patterns from file ff and set their labels to the same unique value e.g. 3333
     unsigned int ppr[11]; for(int k=0; k<11;k++)  ppr[k]=0;   
     subG* s=new subG();
     read_subG(s,ppr,ff);  // reads all subgraphs in file ff and returns them in a linked list pointed to by s
	 subG* temp;
	 nn=0;
	 for (temp=s; (temp) ; temp=temp->next) {
        for (int k=0; k<temp->nol ; k++) temp->labels[k]=3333;
        nn++;
	    }
	 printf("\n The number of patterns read from file %s is %d",ff,nn);
	 // 3rd step: find the embeddings of these labelless patterns
	 for (temp=s; (temp) ; temp=temp->next) {
	    find_subG2(emptyarray,temp);
	    }


}

double cdfN(const double x)  // approximation of cdf of Normal Dist. by Abromowitz and Stegun 
{
  const double b1 =  0.319381530;
  const double b2 = -0.356563782;
  const double b3 =  1.781477937;
  const double b4 = -1.821255978;
  const double b5 =  1.330274429;
  const double p  =  0.2316419;
  const double c  =  0.39894228;

  if(x >= 0.0) {
      double t = 1.0 / ( 1.0 + p * x );
      return (1.0 - c * exp( -x * x / 2.0 ) * t *( t *( t * ( t * ( t * b5 + b4 ) + b3 ) + b2 ) + b1 ));
  }
  else {
      double t = 1.0 / ( 1.0 - p * x );
      return ( c * exp( -x * x / 2.0 ) * t *( t *( t * ( t * ( t * b5 + b4 ) + b3 ) + b2 ) + b1 ));
    }
}

// An implementation of Peter J. Acklam’s algorithm for approximating 
// Inverse cdf of standard normal distribution (by V. Natarajan)
#define  A1  (-3.969683028665376e+01)
#define  A2   2.209460984245205e+02
#define  A3  (-2.759285104469687e+02)
#define  A4   1.383577518672690e+02
#define  A5  (-3.066479806614716e+01)
#define  A6   2.506628277459239e+00

#define  B1  (-5.447609879822406e+01)
#define  B2   1.615858368580409e+02
#define  B3  (-1.556989798598866e+02)
#define  B4   6.680131188771972e+01
#define  B5  (-1.328068155288572e+01)

#define  C1  (-7.784894002430293e-03)
#define  C2  (-3.223964580411365e-01)
#define  C3  (-2.400758277161838e+00)
#define  C4  (-2.549732539343734e+00)
#define  C5   4.374664141464968e+00
#define  C6   2.938163982698783e+00

#define  D1   7.784695709041462e-03
#define  D2   3.224671290700398e-01
#define  D3   2.445134137142996e+00
#define  D4   3.754408661907416e+00

#define P_LOW   0.02425
/* P_high = 1 - p_low*/
#define P_HIGH  0.97575

long double normsinv(double p)
{
long double x;
long double q, r;
if ((0 < p )  && (p < P_LOW)){
   q = sqrt(-2*log(p));
   x = (((((C1*q+C2)*q+C3)*q+C4)*q+C5)*q+C6) / ((((D1*q+D2)*q+D3)*q+D4)*q+1);
}
else{
        if ((P_LOW <= p) && (p <= P_HIGH)){
           q = p - 0.5;
           r = q*q;
           x = (((((A1*r+A2)*r+A3)*r+A4)*r+A5)*r+A6)*q /(((((B1*r+B2)*r+B3)*r+B4)*r+B5)*r+1);
        }
        else{
                if ((P_HIGH < p)&&(p < 1)){
                   q = sqrt(-2*log(1-p));
                   x = -(((((C1*q+C2)*q+C3)*q+C4)*q+C5)*q+C6) / ((((D1*q+D2)*q+D3)*q+D4)*q+1);
                }
        }
}
return x;
}

long double bonferroni_correct_Z_score(double Zscore, int experiments){
    double i=0;
	//long double dara=normsinv((long double)((1.0-cdfN(Zscore))*experiments));
    double p=(double)1.0-cdfN(Zscore);
	//std::numeric_limits<double>::min();
	// cdfN()returns zero if Z-score is greater than 8, therefore
	if (p==0) p=(double)1.0-cdfN(8.25);   
	//printf("\np:%.20f",p);
	double pp=p*experiments;
	if (((pp)>0.0)&&((pp)<0.5))
      return -normsinv((long double)(pp));
	else if (pp<1.0){
	  return -normsinv(((pp-0.5)>(1.0-pp))?(pp-0.5):(1.0-pp));
	  }
	else {
	  for(i=0.5 ; (pp-i)>0.5 ; i=i+0.5);
	  return -normsinv((long double)(pp-i));
	  }
}

int main(int argc, char** argv)
{   /*
	double p=cdfN(5.0*10000000000/0.0000000001);
    printf("\n%.20f",733*(1.0-(p)));
    long double b=-1.0*normsinv(0.99);  //(long double)((733)*(1.0-(p)))
    printf("\n%.20Le\n",b);	
	printf("\nBonferroni corrected:%.20f\n",bonferroni_correct_Z_score(5.0*10000000000/0.0000000001,1));

	double z=0;
	for (double i=0.0; i<5.0 ; i=i+0.1)
		printf("\nZ:%.6f , corrected:%.6f . %d",i,z=bonferroni_correct_Z_score(i,1000),((z-i)>0.00001)?333:0);
	if (1.0==p) printf("\n Yes 1.0");
	*/
	//if (strcmp(argv[1],"test")==0) {test2();exit(1);}
    //test9();
	//printf("\n%d %d",hashf("P35536",19200), ((rand()*rand())%10000));
	// exit(1);
	printf("\nArguments:%d",argc);
	if (argc<4) {printf("\nSAMPLE INVOCATION:\n"
	  "goinpin7_sp2.exe human_ophid_go_annotations.txt human_ophid_PIN_annotated.txt gene_ontology_edit.obo.txt"
                       "subG3.txt 100 8 0 exclude use_embeddings toberemoved.txt Z RGM=F\n");
                 exit(1);
	            }
	char buffer[100];
	char p1[10];
	char p2[10];
	float conf;
	char label[10];
	unsigned int n=0;
	unsigned int hv;
	for(int i=0; argv[i] ; i++) {strcat(command,argv[i]); strcat(command," ");}
	toplist=new subGlist(1000); 
    // read pf file
	printf("\nReading PROTEIN file %s ...", argv[1]);
	FILE* pf=fopen(argv[1], "r");
	fgets( buffer, 100, pf ); 
	sscanf( buffer, "%10s ", p1 );
	printf("\none protein: %s \n",p1);
	while (!feof(pf)) {
	   sscanf( buffer, "%10s ", p1 );
       n++;
	   fgets( buffer, 100, pf );  
	   }
	printf("\none last protein: %s \n",p1);
	printf("Annotations count:%d\n",n);
	printf("size of int data:%d\n Asumed size is 4 . \n",sizeof(int));
	fclose(pf);   
    hts=n*5;   // hash table size
	pht=new PROINFO[hts];
	for(unsigned int i=0; i<hts;i++) {pht[i].id=NULL; pht[i].label=0; pht[i].non=0; pht[i].nl=NULL;
	                                  for(int j=0; j<5 ; j++) pht[i].otherlabels[j]=0;}
    n=0;
	pf=fopen(argv[1], "r");
	fgets( buffer, 100, pf );
	sscanf( buffer, "%10s ", p1 );
	printf("one protein: %s  \n",p1);
	printf("hash value:%d   hts: %d\n",(hashf(p1,hts)), hts); // first protein in the list is not added to the hash table
	while (!feof(pf)) {
	   sscanf( buffer, "%10s%10s ", p1, label );
	   hv=hashf(p1,hts);  // calculate hash value of p1
	   //printf("protein: % s   hash value:%d label:%s\n",p1,hv,label);
	   add_protein(pht,p1,label, hts);
	   n++;
	   fgets( buffer, 100, pf );  
	  }
	fclose(pf);
	
	//test hash table
	unsigned int xx=get_protein(pht,"P00358",hts);
	if (xx<hts) printf("label of %s is %d \n",pht[xx].id, pht[xx].label);
	// Test OK

	
	// read the pin file and update the neighbour list of each protein in the hash table.
	n=0;
	if (argv[2]) { 
	  printf("\nReading PIN file %s ...", argv[2]);
	  FILE* pinf=fopen(argv[2], "r");
	  fgets( buffer, 100, pinf ); 
	  sscanf( buffer, "%10s%10s%f ", p1, p2,&conf );
  	  while (!feof(pinf)) {
	     sscanf( buffer, "%10s%10s%f ", p1, p2,&conf );
         insert_assoc(pht,p1,p2,conf);
		 //printf("%s  %s %f %d\n",p1,p2,conf,n);// last line in the file is somehow processed twice
         n++;
	     fgets( buffer, 100, pinf ); 
	  }
	  printf("%s  %s %f \n",p1,p2,conf);
	  printf("Inserted Associations:%d\n",n);
	  fclose(pinf);
	  pinn=n;
	  }
	fflush(stdout);
    srand( (unsigned)time(NULL) );
	clean_graph();
	if (strcmp(argv[11],"P")==0) ordering_value=PVALUE;
	//printf("\nypl010w=%d",get_protein(pht,"ypl010w",hts));
	if (strcmp(argv[10],"not")) {   // if there is a list of proteins to be removed from pht
		 remove_proteins(argv[10]);  // remove them prior to discovery
	     printf("\nRemoved proteins in %s from pht",argv[10]);}
	if (argv[7]) { // the level parameter , either number or go subset name as in obo file
		if (argv[7][0]=='0') { tlevel=0; strcpy(go_subset,"");}
	  else if (atoi(argv[7])==0) {/*go_subset=new char(40);*/ strcpy(go_subset,argv[7]);tlevel=0;}
	  else { tlevel=atoi(argv[7]); strcpy(go_subset,"");}
	}
	printf("\ntlevel=%d , go_subset=%s",tlevel,go_subset);  
	if (argv[3])    
		if (!read_goterms(argv[3], 100000))
			printf("\nSomething went wrong with obo file parsing.."); 
	//printf("\nFinished with obo file parsing.."); 
	set_gtt_levels(3674);  
	if ((tlevel)||(strcmp(go_subset,""))) 
		update_labels3(tlevel);  // level zero is the root level 
	fflush(stdout);
	pht_vs_gtt();  // compatibility test
	//test8();
	//exit(1);
	//test5();
    if (argv[5]) support=atoi(argv[5]);
	if (argv[6]) maxdepth=atoi(argv[6]);
    lpt= new labp[pinn * 2];
	memset(lpt,0,pinn * 2 * sizeof(labp));
	calc_label_pair_freq();
	if (strcmp(argv[8],"exclude")==0)  exclusions();  // some edges are excluded from the solution
	else printf("\nNo exclusions are made.");
    if (strstr(argv[9],"use")) use_embeddings=1;
	//exit(1);
	int nn=((rand()*rand()*17191)%1000 );
	sprintf(maxresult_file,"maxsubs_s%d_d%d_L%d_ex%d_ue%d_%d.txt",support,maxdepth,tlevel,exclude,use_embeddings,nn);
	//n=((rand()*rand()*17191)%1000 );
    sprintf(resfile,"results_s%d_d%d_L%d_ex%d_ue%d_%d.txt",support,maxdepth,tlevel,exclude,use_embeddings,nn);
    sprintf(rankfile,"rank_s%d_d%d_L%d_ex%d_ue%d_%d.txt",support,maxdepth,tlevel,exclude,use_embeddings,nn);

	if (strstr(argv[4],"rdp:")) {  //read subGs from file argv[4], remove duplicates and write to file "rdp_xxx"
	   if (strlen(argv[4])>5){
		  remove_duplicate_patterns(argv[4]+4,false);
		  exit(1);
	      }
	   } 
	else if (strstr(argv[4],"bonfer:")) {  //read subGs from file argv[4], remove duplicates and write to file "rdp_xxx"
	   if (strlen(argv[4])>8){
		  bonferroni_correct_file(argv[4]+7);
		  exit(1);
	      }
	   } 
	else if (strstr(argv[4],"intersect:")) {  
       char* gg=strchr(argv[4],'-');
	   (*gg)='\0';
	   gg++;
       ;
	   printf("Files %s and %s have %d intersecting patterns. ",argv[4]+10,gg,find_patterns_in_both(argv[4]+10,gg));
	   exit(1);
	   } 
	else if (strstr(argv[4],"fenp:")) {  //read subGs from file argv[4], 
		                                 //find embeddings of non-labelled patterns in this file
	   if (strlen(argv[4])>5){
		  find_embeddings_of_non_labelled_patterns(argv[4]+5);
		  exit(1);
	      }
	   } 
	else
	   printf("\nThis subG has %d nonoverlaping embeddings",nn=read_subG_find(argv[4],pht));
	//exit(1);
	if (argv[12]) RGM=argv[12];
	if (ordering_value==ZSCORE)
		randomPINcount=create_randomPINS(pht,10,15000); 
	
    //random_graph_statistics();

	//printf("\nThis subG has %d nonoverlaping embeddings.",read_subG_find(argv[4],randomPIN[0]));
/*	
	subG* s=new subG();
	unsigned int ppr[11]; for(int k=0; k<11;k++)  ppr[k]=0; 
    read_subG(s,ppr,argv[4]);  // reads all subgraphs in file ff and returns them in a linked list pointed to by s
	printf("\nIts Z-score is %f .",compute_Z_score(s));
*/  
	if (strcmp(argv[13],"FIARN")==0) {
      // Lets see if a random network contains frequent patterns
	  PROINFO* pin1=pht;
	  pht=randomPIN[0];
	  randomPIN[0]=pin1;
	  printf("\n Real network abondoned: Capturing frequent patterns in one RANDOM NETWORK!!!");
	  }
	printf("\nRunning all_subgraphs() : ");
	all_subgraphs();

   if (strcmp(argv[4],"rdp:")==0) //read subGs from maxresult_file, remove duplicates and write to file rdp_maxresult_file
remove_duplicate_patterns(maxresult_file,true);


	delete pht;
	end_pgm(argv);

	return 0;


}



