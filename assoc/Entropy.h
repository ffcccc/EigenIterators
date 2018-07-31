#ifndef ENTROPY_H_
#define ENTROPY_H_ 1

#include <iostream>
#include <vector>
#include <map>
#include <string>
#include <cassert>
#include <valarray>

using namespace std;

typedef enum {X, Y} SVar;
typedef enum {A, B, C, D} TVal;
//typedef enum {TP, FP, TN, FN} TVal;

template<class A>
class Mapper{
public:
	// costruttore da array di valori da mappare
	Mapper(const valarray<A> &X_);
	// restituisce l'indice corrispondente al valore 'value'
	inline size_t getPos(const A &value) const {return XSet.find(value)->second;};
	// restituisce il valore associato alla posizione 'pos'
	A getEvent(const size_t pos) const;
	// restituisce la cardinalit della variabile X
	inline size_t size() const {return XSet.size();};
private:
	map<A, size_t> XSet;
};

template<class A>
Mapper<A>::Mapper(const valarray<A> &X_){
	XSet.clear();
	for (unsigned int j = 0; j < X_.size(); j++) {
		XSet.insert(make_pair(X_[j], XSet.size()));
	}
}

//template<class A>
//size_t Mapper<A>::getPos(const A &value) const {
//	return XSet.find(value)->second;
//}

template<class A>
A Mapper<A>::getEvent(const size_t pos) const {
	typename map<A,size_t>::const_iterator iter = XSet.begin();
	size_t cont = 0;
	while(cont++ < pos) iter++;
	return iter->first;
}

//template<class A>
//size_t Mapper<A>::size() const {
//	return XSet.size();
//}

template<class T, class A>
class ProbBox{

	class MDRCells{
	public:
		size_t size() const { return _Cell.size();};
		void resize(size_t n){
			_Cell.resize(n);
			for(unsigned int j=0; j<n; j++)
				_Cell[j].resize(4, T(0));
		};
		void run(size_t pos, bool cc){
			for(unsigned int j=0; j<_Cell.size(); j++){
				bool exposed(pos==j);
				if(exposed){
					if(cc)
						_Cell[j][0]++;	//A
					else
						_Cell[j][1]++;	//B
				}else {
					if(cc)
						_Cell[j][2]++;	//C
					else
						_Cell[j][3]++;	//D
				}
			}
		};
		std::vector< valarray<T> > _Cell;
	};

  public:
	ProbBox();  
	//ProbBox(const A *X_, const A *Y_); // per performance !
	ProbBox(const valarray<A> &X_, const valarray<A> &Y_);
	//ProbBox(const A *X_, const A *Y_, const valarray<T> &CT_); // per performance !
	ProbBox(const valarray<A> &X_, const valarray<A> &Y_, const valarray<T> &CT_);
	~ProbBox();
	
	void updateClassVar(const valarray<int> &CC_);
    // codifica var. rimane uguale, cambiano i vettori dati (eventi nei samples) e aggiorno tabelle 
	void updateContingencyTable(const valarray<A> &X_, const valarray<A> &Y_);
	// codifica var. rimane uguale, cambiano i vettori dati, ma carico la tabella dei conteggi direttamente 
	void updateContingencyTable(const valarray<T> &CT_);
	// aggiorno la codifica delle var.
	void updateVarMaps(const valarray<A> &X_, const valarray<A> &Y_);

    T prob(SVar var, A event);
    T jointProbXY(A eventX, A eventY);
    T conditionalProbYX(A eventX, A eventY);
    T entropy(SVar var);
    T jointEntropyXY();
    T conditionalEntropy();
    T mutualInformation();
	//
	inline T N() {return _N;};
	T freq(SVar var, A event);
	T expected(A eventX, A eventY);
	T observed(A eventX, A eventY);
	T chi2();
	T LR();
	T MDR();

  private:
	  static inline T invLog2()	{return T(3.32192809488736234791); }; // 1.0/log10(2.0);};		
	  static inline T Log2()	{return T(0.30102999566398119521); }; // log10(2.0);};
	
	void	doJointFreqTable(const valarray<A> &X_, const valarray<A> &Y_, MDRCells &cells);
	//void	doJointProbTable();
	size_t	getIndex(const A &eventX, const A &eventY);
	std::slice* dataSlice(SVar var, A event);
	Mapper<A>* Var[2];
	// array per contenere la tavola di contingenza delle 2 variabili (conteggio degli eventi)
	valarray<T> jFT;
	// array per contenere la tavola delle probabilit� marginali (conteggi / N)
	valarray<T> jPT;
	// numero totale di eventi della tab. contingenza
	T _N;
	// variabile di classe: caso=1 e controllo=0
	valarray<int> CC;
	MDRCells Cells;

};

template<class T, class A> inline
ProbBox<T,A>::ProbBox(){
	Var[X]=NULL;
	Var[Y]=NULL;
	_N = T(0);
	Cells.resize(0);
}

template<class T, class A> inline
ProbBox<T,A>::ProbBox(const valarray<A> &X_, const valarray<A> &Y_){
	Var[X]=NULL;
	Var[Y]=NULL;
	_N = T(0);
	Cells.resize(0);
	
	assert(X_.size() == Y_.size());
	updateVarMaps(X_, Y_);
	updateContingencyTable(X_, Y_);
}

template<class T, class A> inline
ProbBox<T,A>::ProbBox(const valarray<A> &X_, const valarray<A> &Y_, const valarray<T> &CT_){
	Var[X]=NULL;
	Var[Y]=NULL;
	_N = T(0);
	Cells.resize(0);

	assert(X_.size() == Y_.size());
	assert(CT_.size() == X_.size() * Y_.size());
	updateVarMaps(X_, Y_);
	updateContingencyTable(CT_);
}

template<class T, class A> inline
ProbBox<T,A>::~ProbBox(){
//	cout << " -pb destroy- ";
	if(Var[X]!=NULL) delete Var[X];
	if(Var[Y]!=NULL) delete Var[Y];
}

template<class T, class A> inline
void ProbBox<T,A>::updateVarMaps(const valarray<A> &X_, const valarray<A> &Y_){
	assert(X_.size() == Y_.size());
	if(Var[X]!=NULL) delete Var[X];
	if(Var[Y]!=NULL) delete Var[Y];
	Var[X] = new Mapper<A>(X_);
	Var[Y] = new Mapper<A>(Y_);
//	updateContingencyTable(X_, Y_);
}

template<class T, class A> inline
void ProbBox<T,A>::updateClassVar(const valarray<int> &CC_){
	CC = CC_;
}

template<class T, class A> inline
void ProbBox<T,A>::updateContingencyTable(const valarray<A> &X_, const valarray<A> &Y_){
	assert(X_.size() == Y_.size());
	_N = T(X_.size());
	jFT.resize(Var[X]->size() * Var[Y]->size(), T(0));
	jPT.resize(Var[X]->size() * Var[Y]->size(), T(0));
	doJointFreqTable(X_, Y_, Cells);
	//doJointProbTable();
}

template<class T, class A> inline
void ProbBox<T,A>::updateContingencyTable(const valarray<T> &CT_){
	assert(CT_.size() == Var[X]->size() * Var[Y]->size());
	jFT.resize(Var[X]->size() * Var[Y]->size(), T(0));
	jPT.resize(Var[X]->size() * Var[Y]->size(), T(0));
	jFT = CT_;
	_N = jFT.sum();
	//doJointProbTable();
	jPT = jFT / N();
}

template<class T, class A> inline
size_t ProbBox<T,A>::getIndex(const A &eventX, const A &eventY){
	size_t pos = Var[Y]->getPos(eventY) * Var[X]->size() + Var[X]->getPos(eventX);
	return pos;
}

template<class T, class A> inline
void ProbBox<T,A>::doJointFreqTable(const valarray<A> &X_, const valarray<A> &Y_, MDRCells &cells){
	assert(CC.size() == X_.size());
	cells.resize(jFT.size());
	for(unsigned int i=0; i<X_.size(); i++){
		size_t pos = getIndex(X_[i], Y_[i]);	//int(Y_[i])*Var[X]->size() + int(X_[i]);
		assert(pos >= 0);
		assert(pos <  jPT.size());
		jFT[pos]++;
		// callback !
		cells.run(pos, CC[i] != 0);
	}
	// doJointProbTable
	assert(_N == jFT.sum());
	jPT = jFT / N();
}

//template<class T, class A> inline
//void ProbBox<T,A>::doJointProbTable(){
//	assert(_N == jFT.sum());
//	jPT = jFT / N();
//}

template<class T, class A> inline
slice * ProbBox<T,A>::dataSlice(SVar var, A event){
	slice *vaSlice;
	if(var == X){
		vaSlice = new slice(Var[X]->getPos(event) , Var[Y]->size() , Var[X]->size() );
	} else {
		vaSlice = new slice(Var[X]->size()*Var[Y]->getPos(event), Var[X]->size(), 1);
	}
	return vaSlice;
}

template<class T, class A> inline
T ProbBox<T,A>::prob(SVar var, A event){
	slice *vaSlice = dataSlice(var, event);
	valarray<T> temp = jPT[*vaSlice];
	T result = temp.sum();
	delete vaSlice;
	return result;
}

template<class T, class A> inline
T ProbBox<T,A>::freq(SVar var, A event){
	slice *vaSlice = dataSlice(var, event);
	valarray<T> temp = jFT[*vaSlice];
	T result = temp.sum();
	delete vaSlice;
	return result;
}

template<class T, class A> inline
T ProbBox<T,A>::expected(A eventX, A eventY){
	T result = ( freq(X, eventX) * freq(Y, eventY) ) / N();
	return result;
}

template<class T, class A> inline
T ProbBox<T,A>::observed(A eventX, A eventY){
	unsigned int pos = getIndex(eventX, eventY);	//eventY*Var[X]->size() + eventX;
	assert(pos >= 0);
	assert(pos <  jFT.size());
	return jFT[pos];
}

template<class T, class A> inline
T ProbBox<T,A>::jointProbXY(A eventX, A eventY){
	size_t pos = getIndex(eventX, eventY);	//eventY*Var[X]->size() + eventX;
	assert(pos >= 0);
	assert(pos <  jPT.size());
	return jPT[pos];
}

template<class T, class A> inline
T ProbBox<T,A>::conditionalProbYX(A eventX, A eventY){
	return jointProbXY(eventX, eventY) / prob(X, eventX);
}

template<class T, class A> inline
T ProbBox<T,A>::entropy(SVar var){
	T result = T(0);
	T Pb = T(0);
	for(unsigned int i=0;i<Var[var]->size();i++){
		Pb = prob(var, Var[var]->getEvent(i));
		result += Pb * (log10(Pb)*invLog2());	// log base 2 in funz. di logaritmi base 10
	}
	return -result;
}

template<class T, class A> inline
T ProbBox<T,A>::jointEntropyXY(){
	T result = T(0); 
	for(unsigned int i=0;i<jPT.size();i++){
		if(jPT[i] > 0.01){
			result += jPT[i]*(log10(jPT[i])*invLog2());	// log base 2 in funz. di logaritmi base 10
		}
	}
	return -result;
}

//template<class T>
//T ProbBox<T,A>::conditionalEntropy(){}

template<class T, class A> inline
T ProbBox<T,A>::mutualInformation(){
	T a=entropy(X);
	T b=entropy(Y);
	T c=jointEntropyXY();
	T result = a+b-c;	//entropyX()+entropyY()-jointEntropyXY(); 
	return result;
}

template<class T, class A> inline
T ProbBox<T,A>::chi2(){
	T result = T(0);
	for(unsigned int x=0;x<Var[X]->size();x++){
		A evX = Var[X]->getEvent(x);
		for(unsigned int y=0;y<Var[Y]->size();y++){
			A evY = Var[Y]->getEvent(y);
			T expxy = expected(evX,evY);
			result += pow(T(observed(evX,evY)-expxy), T(2.0))/expxy;
		}
	}
		
	return result;
}

template<class T, class A> inline
T ProbBox<T,A>::LR(){
	T result = T(0);
	for(unsigned int x=0;x<Var[X]->size();x++){
		//for each( pair<A, size_t> x in Var[X]->XSet){
		A evX = Var[X]->getEvent(x);
		for(unsigned int y=0;y<Var[Y]->size();y++){
			//for each( pair<A, size_t> y in Var[Y]->XSet){
			A evY = Var[Y]->getEvent(y);
			T expxy = expected(evX, evY);
			T obsxy = observed(evX, evY);
			result += obsxy * log(obsxy/expxy);	// logaritmo naturale
		}
	}
		
	return result*2;
}

template<class T, class A> inline
T ProbBox<T,A>::MDR(){
	T result = T(0);
	valarray<string> vExpos(2);		vExpos[0]="exp";	vExpos[1]="notExp";
	valarray<string> vStatus(2);	vStatus[0]="case";	vStatus[1]="control";
	valarray<T> newvar(T(0), Cells.size());
	ProbBox<T, string> pb;
	pb.updateVarMaps(vExpos, vStatus);

	for(unsigned int c=0;c<Cells.size();c++){
		pb.updateContingencyTable(Cells._Cell[c]);
		//if((Cells[c][0]/Cells[c][2]) > (Cells[c][1]/Cells[c][3])){
		//	newvar[c]=T(1);
		//}
		result = pb.chi2();
	}
	return result;
}

#endif /*ENTROPY_H_*/
