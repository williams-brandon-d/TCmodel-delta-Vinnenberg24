/* Created by Language version: 7.7.0 */
/* NOT VECTORIZED */
#define NRN_VECTORIZED 0
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "scoplib_ansi.h"
#undef PI
#define nil 0
#include "md1redef.h"
#include "section.h"
#include "nrniv_mf.h"
#include "md2redef.h"
 
#if METHOD3
extern int _method3;
#endif

#if !NRNGPU
#undef exp
#define exp hoc_Exp
extern double hoc_Exp(double);
#endif
 
#define nrn_init _nrn_init__GABAa_S
#define _nrn_initial _nrn_initial__GABAa_S
#define nrn_cur _nrn_cur__GABAa_S
#define _nrn_current _nrn_current__GABAa_S
#define nrn_jacob _nrn_jacob__GABAa_S
#define nrn_state _nrn_state__GABAa_S
#define _net_receive _net_receive__GABAa_S 
#define release release__GABAa_S 
 
#define _threadargscomma_ /**/
#define _threadargsprotocomma_ /**/
#define _threadargs_ /**/
#define _threadargsproto_ /**/
 	/*SUPPRESS 761*/
	/*SUPPRESS 762*/
	/*SUPPRESS 763*/
	/*SUPPRESS 765*/
	 extern double *getarg();
 static double *_p; static Datum *_ppvar;
 
#define t nrn_threads->_t
#define dt nrn_threads->_dt
#define gmax _p[0]
#define i _p[1]
#define g _p[2]
#define C _p[3]
#define R _p[4]
#define R0 _p[5]
#define R1 _p[6]
#define TimeCount _p[7]
#define lastrelease _p[8]
#define _g _p[9]
#define _nd_area  *_ppvar[0]._pval
#define pre	*_ppvar[2]._pval
#define _p_pre	_ppvar[2]._pval
 
#if MAC
#if !defined(v)
#define v _mlhv
#endif
#if !defined(h)
#define h _mlhh
#endif
#endif
 
#if defined(__cplusplus)
extern "C" {
#endif
 static int hoc_nrnpointerindex =  2;
 /* external NEURON variables */
 /* declaration of user functions */
 static double _hoc_exptable();
 static double _hoc_release();
 static int _mechtype;
extern void _nrn_cacheloop_reg(int, int);
extern void hoc_register_prop_size(int, int, int);
extern void hoc_register_limits(int, HocParmLimits*);
extern void hoc_register_units(int, HocParmUnits*);
extern void nrn_promote(Prop*, int, int);
extern Memb_func* memb_func;
 
#define NMODL_TEXT 1
#if NMODL_TEXT
static const char* nmodl_file_text;
static const char* nmodl_filename;
extern void hoc_reg_nmodl_text(int, const char*);
extern void hoc_reg_nmodl_filename(int, const char*);
#endif

 extern Prop* nrn_point_prop_;
 static int _pointtype;
 static void* _hoc_create_pnt(_ho) Object* _ho; { void* create_point_process();
 return create_point_process(_pointtype, _ho);
}
 static void _hoc_destroy_pnt();
 static double _hoc_loc_pnt(_vptr) void* _vptr; {double loc_point_process();
 return loc_point_process(_pointtype, _vptr);
}
 static double _hoc_has_loc(_vptr) void* _vptr; {double has_loc_point();
 return has_loc_point(_vptr);
}
 static double _hoc_get_loc_pnt(_vptr)void* _vptr; {
 double get_loc_point_process(); return (get_loc_point_process(_vptr));
}
 extern void _nrn_setdata_reg(int, void(*)(Prop*));
 static void _setdata(Prop* _prop) {
 _p = _prop->param; _ppvar = _prop->dparam;
 }
 static void _hoc_setdata(void* _vptr) { Prop* _prop;
 _prop = ((Point_process*)_vptr)->_prop;
   _setdata(_prop);
 }
 /* connect user functions to hoc names */
 static VoidFunc hoc_intfunc[] = {
 0,0
};
 static Member_func _member_func[] = {
 "loc", _hoc_loc_pnt,
 "has_loc", _hoc_has_loc,
 "get_loc", _hoc_get_loc_pnt,
 "exptable", _hoc_exptable,
 "release", _hoc_release,
 0, 0
};
#define _f_exptable _f_exptable_GABAa_S
#define exptable exptable_GABAa_S
 extern double _f_exptable( double );
 extern double exptable( double );
 /* declare global and static user variables */
#define Alpha Alpha_GABAa_S
 double Alpha = 10.5;
#define Beta Beta_GABAa_S
 double Beta = 0.166;
#define Cdur Cdur_GABAa_S
 double Cdur = 0.3;
#define Cmax Cmax_GABAa_S
 double Cmax = 0.5;
#define Deadtime Deadtime_GABAa_S
 double Deadtime = 1;
#define Erev Erev_GABAa_S
 double Erev = -80;
#define Prethresh Prethresh_GABAa_S
 double Prethresh = 0;
#define Rtau Rtau_GABAa_S
 double Rtau = 0;
#define Rinf Rinf_GABAa_S
 double Rinf = 0;
#define usetable usetable_GABAa_S
 double usetable = 1;
 /* some parameters have upper and lower limits */
 static HocParmLimits _hoc_parm_limits[] = {
 "usetable_GABAa_S", 0, 1,
 0,0,0
};
 static HocParmUnits _hoc_parm_units[] = {
 "Cmax_GABAa_S", "mM",
 "Cdur_GABAa_S", "ms",
 "Alpha_GABAa_S", "/ms",
 "Beta_GABAa_S", "/ms",
 "Erev_GABAa_S", "mV",
 "Deadtime_GABAa_S", "ms",
 "Rtau_GABAa_S", "ms",
 "gmax", "umho",
 "i", "nA",
 "g", "umho",
 "C", "mM",
 "TimeCount", "ms",
 0,0
};
 static double delta_t = 1;
 static double v = 0;
 /* connect global user variables to hoc */
 static DoubScal hoc_scdoub[] = {
 "Cmax_GABAa_S", &Cmax_GABAa_S,
 "Cdur_GABAa_S", &Cdur_GABAa_S,
 "Alpha_GABAa_S", &Alpha_GABAa_S,
 "Beta_GABAa_S", &Beta_GABAa_S,
 "Erev_GABAa_S", &Erev_GABAa_S,
 "Prethresh_GABAa_S", &Prethresh_GABAa_S,
 "Deadtime_GABAa_S", &Deadtime_GABAa_S,
 "Rinf_GABAa_S", &Rinf_GABAa_S,
 "Rtau_GABAa_S", &Rtau_GABAa_S,
 "usetable_GABAa_S", &usetable_GABAa_S,
 0,0
};
 static DoubVec hoc_vdoub[] = {
 0,0,0
};
 static double _sav_indep;
 static void nrn_alloc(Prop*);
static void  nrn_init(_NrnThread*, _Memb_list*, int);
static void nrn_state(_NrnThread*, _Memb_list*, int);
 static void nrn_cur(_NrnThread*, _Memb_list*, int);
static void  nrn_jacob(_NrnThread*, _Memb_list*, int);
 static void _hoc_destroy_pnt(_vptr) void* _vptr; {
   destroy_point_process(_vptr);
}
 
static int _ode_count(int);
 /* connect range variables in _p that hoc is supposed to know about */
 static const char *_mechanism[] = {
 "7.7.0",
"GABAa_S",
 "gmax",
 0,
 "i",
 "g",
 "C",
 "R",
 "R0",
 "R1",
 "TimeCount",
 0,
 0,
 "pre",
 0};
 
extern Prop* need_memb(Symbol*);

static void nrn_alloc(Prop* _prop) {
	Prop *prop_ion;
	double *_p; Datum *_ppvar;
  if (nrn_point_prop_) {
	_prop->_alloc_seq = nrn_point_prop_->_alloc_seq;
	_p = nrn_point_prop_->param;
	_ppvar = nrn_point_prop_->dparam;
 }else{
 	_p = nrn_prop_data_alloc(_mechtype, 10, _prop);
 	/*initialize range parameters*/
 	gmax = 0;
  }
 	_prop->param = _p;
 	_prop->param_size = 10;
  if (!nrn_point_prop_) {
 	_ppvar = nrn_prop_datum_alloc(_mechtype, 3, _prop);
  }
 	_prop->dparam = _ppvar;
 	/*connect ionic variables to this model*/
 
}
 static void _initlists();
 extern Symbol* hoc_lookup(const char*);
extern void _nrn_thread_reg(int, int, void(*)(Datum*));
extern void _nrn_thread_table_reg(int, void(*)(double*, Datum*, Datum*, _NrnThread*, int));
extern void hoc_register_tolerance(int, HocStateTolerance*, Symbol***);
extern void _cvode_abstol( Symbol**, double*, int);

 void _gabaa_reg() {
	int _vectorized = 0;
  _initlists();
 	_pointtype = point_register_mech(_mechanism,
	 nrn_alloc,nrn_cur, nrn_jacob, nrn_state, nrn_init,
	 hoc_nrnpointerindex, 0,
	 _hoc_create_pnt, _hoc_destroy_pnt, _member_func);
 _mechtype = nrn_get_mechtype(_mechanism[1]);
     _nrn_setdata_reg(_mechtype, _setdata);
 #if NMODL_TEXT
  hoc_reg_nmodl_text(_mechtype, nmodl_file_text);
  hoc_reg_nmodl_filename(_mechtype, nmodl_filename);
#endif
  hoc_register_prop_size(_mechtype, 10, 3);
  hoc_register_dparam_semantics(_mechtype, 0, "area");
  hoc_register_dparam_semantics(_mechtype, 1, "pntproc");
  hoc_register_dparam_semantics(_mechtype, 2, "pointer");
 	hoc_register_cvode(_mechtype, _ode_count, 0, 0, 0);
 	hoc_register_var(hoc_scdoub, hoc_vdoub, hoc_intfunc);
 	ivoc_help("help ?1 GABAa_S C:/Users/brndn/OneDrive/Desktop/DLGN_NEW/gabaa.mod\n");
 hoc_register_limits(_mechtype, _hoc_parm_limits);
 hoc_register_units(_mechtype, _hoc_parm_units);
 }
 static double *_t_exptable;
static int _reset;
static char *modelname = "simple GABAa receptors";

static int error;
static int _ninits = 0;
static int _match_recurse=1;
static void _modl_cleanup(){ _match_recurse=1;}
static int release();
 static double _n_exptable(double);
 
static int  release (  ) {
   TimeCount = TimeCount - dt ;
   if ( TimeCount < - Deadtime ) {
     if ( pre > Prethresh ) {
       C = Cmax ;
       R0 = R ;
       lastrelease = t ;
       TimeCount = Cdur ;
       }
     }
   else if ( TimeCount > 0.0 ) {
     }
   else if ( C  == Cmax ) {
     R1 = R ;
     C = 0. ;
     }
   if ( C > 0.0 ) {
     R = Rinf + ( R0 - Rinf ) * exptable ( _threadargscomma_ - ( t - lastrelease ) / Rtau ) ;
     }
   else {
     R = R1 * exptable ( _threadargscomma_ - Beta * ( t - ( lastrelease + Cdur ) ) ) ;
     }
   
/*VERBATIM*/
	return 0;
  return 0; }
 
static double _hoc_release(void* _vptr) {
 double _r;
    _hoc_setdata(_vptr);
 _r = 1.;
 release (  );
 return(_r);
}
 static double _mfac_exptable, _tmin_exptable;
 static void _check_exptable();
 static void _check_exptable() {
  static int _maktable=1; int _i, _j, _ix = 0;
  double _xi, _tmax;
  if (!usetable) {return;}
  if (_maktable) { double _x, _dx; _maktable=0;
   _tmin_exptable =  - 10.0 ;
   _tmax =  10.0 ;
   _dx = (_tmax - _tmin_exptable)/2000.; _mfac_exptable = 1./_dx;
   for (_i=0, _x=_tmin_exptable; _i < 2001; _x += _dx, _i++) {
    _t_exptable[_i] = _f_exptable(_x);
   }
  }
 }

 double exptable(double _lx){ _check_exptable();
 return _n_exptable(_lx);
 }

 static double _n_exptable(double _lx){ int _i, _j;
 double _xi, _theta;
 if (!usetable) {
 return _f_exptable(_lx); 
}
 _xi = _mfac_exptable * (_lx - _tmin_exptable);
 if (isnan(_xi)) {
  return _xi; }
 if (_xi <= 0.) {
 return _t_exptable[0];
 }
 if (_xi >= 2000.) {
 return _t_exptable[2000];
 }
 _i = (int) _xi;
 return _t_exptable[_i] + (_xi - (double)_i)*(_t_exptable[_i+1] - _t_exptable[_i]);
 }

 
double _f_exptable (  double _lx ) {
   double _lexptable;
 if ( ( _lx > - 10.0 )  && ( _lx < 10.0 ) ) {
     _lexptable = exp ( _lx ) ;
     }
   else {
     _lexptable = 0. ;
     }
   
return _lexptable;
 }
 
static double _hoc_exptable(void* _vptr) {
 double _r;
    _hoc_setdata(_vptr);
  _r =  exptable (  *getarg(1) );
 return(_r);
}
 
static int _ode_count(int _type){ hoc_execerror("GABAa_S", "cannot be used with CVODE"); return 0;}

static void initmodel() {
  int _i; double _save;_ninits++;
 _save = t;
 t = 0.0;
{
 {
   R = 0.0 ;
   C = 0.0 ;
   Rinf = Cmax * Alpha / ( Cmax * Alpha + Beta ) ;
   Rtau = 1.0 / ( ( Alpha * Cmax ) + Beta ) ;
   lastrelease = - 9e9 ;
   R1 = 0.0 ;
   TimeCount = - 1.0 ;
   }
  _sav_indep = t; t = _save;

}
}

static void nrn_init(_NrnThread* _nt, _Memb_list* _ml, int _type){
Node *_nd; double _v; int* _ni; int _iml, _cntml;
#if CACHEVEC
    _ni = _ml->_nodeindices;
#endif
_cntml = _ml->_nodecount;
for (_iml = 0; _iml < _cntml; ++_iml) {
 _p = _ml->_data[_iml]; _ppvar = _ml->_pdata[_iml];
#if CACHEVEC
  if (use_cachevec) {
    _v = VEC_V(_ni[_iml]);
  }else
#endif
  {
    _nd = _ml->_nodelist[_iml];
    _v = NODEV(_nd);
  }
 v = _v;
 initmodel();
}}

static double _nrn_current(double _v){double _current=0.;v=_v;{ {
   g = gmax * R ;
   i = g * ( v - Erev ) ;
   }
 _current += i;

} return _current;
}

static void nrn_cur(_NrnThread* _nt, _Memb_list* _ml, int _type){
Node *_nd; int* _ni; double _rhs, _v; int _iml, _cntml;
#if CACHEVEC
    _ni = _ml->_nodeindices;
#endif
_cntml = _ml->_nodecount;
for (_iml = 0; _iml < _cntml; ++_iml) {
 _p = _ml->_data[_iml]; _ppvar = _ml->_pdata[_iml];
#if CACHEVEC
  if (use_cachevec) {
    _v = VEC_V(_ni[_iml]);
  }else
#endif
  {
    _nd = _ml->_nodelist[_iml];
    _v = NODEV(_nd);
  }
 _g = _nrn_current(_v + .001);
 	{ _rhs = _nrn_current(_v);
 	}
 _g = (_g - _rhs)/.001;
 _g *=  1.e2/(_nd_area);
 _rhs *= 1.e2/(_nd_area);
#if CACHEVEC
  if (use_cachevec) {
	VEC_RHS(_ni[_iml]) -= _rhs;
  }else
#endif
  {
	NODERHS(_nd) -= _rhs;
  }
 
}}

static void nrn_jacob(_NrnThread* _nt, _Memb_list* _ml, int _type){
Node *_nd; int* _ni; int _iml, _cntml;
#if CACHEVEC
    _ni = _ml->_nodeindices;
#endif
_cntml = _ml->_nodecount;
for (_iml = 0; _iml < _cntml; ++_iml) {
 _p = _ml->_data[_iml];
#if CACHEVEC
  if (use_cachevec) {
	VEC_D(_ni[_iml]) += _g;
  }else
#endif
  {
     _nd = _ml->_nodelist[_iml];
	NODED(_nd) += _g;
  }
 
}}

static void nrn_state(_NrnThread* _nt, _Memb_list* _ml, int _type){
Node *_nd; double _v = 0.0; int* _ni; int _iml, _cntml;
#if CACHEVEC
    _ni = _ml->_nodeindices;
#endif
_cntml = _ml->_nodecount;
for (_iml = 0; _iml < _cntml; ++_iml) {
 _p = _ml->_data[_iml]; _ppvar = _ml->_pdata[_iml];
 _nd = _ml->_nodelist[_iml];
#if CACHEVEC
  if (use_cachevec) {
    _v = VEC_V(_ni[_iml]);
  }else
#endif
  {
    _nd = _ml->_nodelist[_iml];
    _v = NODEV(_nd);
  }
 v=_v;
{
 { error =  release();
 if(error){fprintf(stderr,"at line 112 in file gabaa.mod:\n	SOLVE release\n"); nrn_complain(_p); abort_run(error);}
 }}}

}

static void terminal(){}

static void _initlists() {
 int _i; static int _first = 1;
  if (!_first) return;
   _t_exptable = makevector(2001*sizeof(double));
_first = 0;
}

#if NMODL_TEXT
static const char* nmodl_filename = "gabaa.mod";
static const char* nmodl_file_text = 
  "TITLE simple GABAa receptors\n"
  "\n"
  "COMMENT\n"
  "-----------------------------------------------------------------------------\n"
  "\n"
  "	Simple model for GABAa receptors\n"
  "	================================\n"
  "\n"
  "  - FIRST-ORDER KINETICS, FIT TO WHOLE-CELL RECORDINGS\n"
  "\n"
  "    Whole-cell recorded GABA-A postsynaptic currents (Otis et al, J. Physiol. \n"
  "    463: 391-407, 1993) were used to estimate the parameters of the present\n"
  "    model; the fit was performed using a simplex algorithm (see Destexhe et\n"
  "    al., J. Neurophysiol. 72: 803-818, 1994).\n"
  "\n"
  "  - SHORT PULSES OF TRANSMITTER (0.3 ms, 0.5 mM)\n"
  "\n"
  "    The simplified model was obtained from a detailed synaptic model that \n"
  "    included the release of transmitter in adjacent terminals, its lateral \n"
  "    diffusion and uptake, and its binding on postsynaptic receptors (Destexhe\n"
  "    and Sejnowski, 1995).  Short pulses of transmitter with first-order\n"
  "    kinetics were found to be the best fast alternative to represent the more\n"
  "    detailed models.\n"
  "\n"
  "  - ANALYTIC EXPRESSION\n"
  "\n"
  "    The first-order model can be solved analytically, leading to a very fast\n"
  "    mechanism for simulating synapses, since no differential equation must be\n"
  "    solved (see references below).\n"
  "\n"
  "\n"
  "\n"
  "References\n"
  "\n"
  "   Destexhe, A., Mainen, Z.F. and Sejnowski, T.J.  An efficient method for\n"
  "   computing synaptic conductances based on a kinetic model of receptor binding\n"
  "   Neural Computation 6: 10-14, 1994.  \n"
  "\n"
  "   Destexhe, A., Mainen, Z.F. and Sejnowski, T.J. Synthesis of models for\n"
  "   excitable membranes, synaptic transmission and neuromodulation using a \n"
  "   common kinetic formalism, Journal of Computational Neuroscience 1: \n"
  "   195-230, 1994.\n"
  "\n"
  "See also:\n"
  "\n"
  "   http://cns.iaf.cnrs-gif.fr\n"
  "\n"
  "Written by A. Destexhe, 1995\n"
  "27-11-2002: the pulse is implemented using a counter, which is more\n"
  "	stable numerically (thanks to Yann LeFranc)\n"
  "\n"
  "-----------------------------------------------------------------------------\n"
  "ENDCOMMENT\n"
  "\n"
  "\n"
  "\n"
  "INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}\n"
  "\n"
  "NEURON {\n"
  "	POINT_PROCESS GABAa_S\n"
  "	POINTER pre\n"
  "	RANGE C, R, R0, R1, g, gmax, TimeCount\n"
  "	NONSPECIFIC_CURRENT i\n"
  "	GLOBAL Cmax, Cdur, Alpha, Beta, Erev, Prethresh, Deadtime, Rinf, Rtau\n"
  "}\n"
  "UNITS {\n"
  "	(nA) = (nanoamp)\n"
  "	(mV) = (millivolt)\n"
  "	(umho) = (micromho)\n"
  "	(mM) = (milli/liter)\n"
  "}\n"
  "\n"
  "PARAMETER {\n"
  "	dt		(ms)\n"
  "	Cmax	= 0.5	(mM)		: max transmitter concentration\n"
  "	Cdur	= 0.3	(ms)		: transmitter duration (rising phase)\n"
  "	Alpha	= 10.5	(/ms mM)	: forward (binding) rate\n"
  "	Beta	= 0.166	(/ms)		: backward (unbinding) rate\n"
  "	Erev	= -80	(mV)		: reversal potential\n"
  "	Prethresh = 0 			: voltage level nec for release\n"
  "	Deadtime = 1	(ms)		: mimimum time between release events\n"
  "	gmax		(umho)		: maximum conductance\n"
  "}\n"
  "\n"
  "\n"
  "ASSIGNED {\n"
  "	v		(mV)		: postsynaptic voltage\n"
  "	i 		(nA)		: current = g*(v - Erev)\n"
  "	g 		(umho)		: conductance\n"
  "	C		(mM)		: transmitter concentration\n"
  "	R				: fraction of open channels\n"
  "	R0				: open channels at start of release\n"
  "	R1				: open channels at end of release\n"
  "	Rinf				: steady state channels open\n"
  "	Rtau		(ms)		: time constant of channel binding\n"
  "	pre 				: pointer to presynaptic variable\n"
  "	lastrelease	(ms)		: time of last spike\n"
  "	TimeCount	(ms)		: time counter\n"
  "}\n"
  "\n"
  "INITIAL {\n"
  "	R = 0\n"
  "	C = 0\n"
  "	Rinf = Cmax*Alpha / (Cmax*Alpha + Beta)\n"
  "	Rtau = 1 / ((Alpha * Cmax) + Beta)\n"
  "	lastrelease = -9e9\n"
  "	R1=0\n"
  "	TimeCount=-1\n"
  "}\n"
  "\n"
  "BREAKPOINT {\n"
  "	SOLVE release\n"
  "	g = gmax * R\n"
  "	i = g*(v - Erev)\n"
  "}\n"
  "\n"
  "PROCEDURE release() { \n"
  "	:will crash if user hasn't set pre with the connect statement \n"
  "\n"
  "	TimeCount = TimeCount - dt		: time since last release ended\n"
  "\n"
  "						: ready for another release?\n"
  "	if (TimeCount < -Deadtime) {\n"
  "		if (pre > Prethresh) {		: spike occured?\n"
  "			C = Cmax			: start new release\n"
  "			R0 = R\n"
  "			lastrelease = t\n"
  "			TimeCount=Cdur\n"
  "		}\n"
  "						\n"
  "        } else if (TimeCount > 0) {		: still releasing?\n"
  "	\n"
  "		: do nothing\n"
  "	\n"
  "	} else if (C == Cmax) {			: in dead time after release\n"
  "		R1 = R\n"
  "		C = 0.\n"
  "	}\n"
  "\n"
  "\n"
  "\n"
  "	if (C > 0) {				: transmitter being released?\n"
  "\n"
  "	   R = Rinf + (R0 - Rinf) * exptable (- (t - lastrelease) / Rtau)\n"
  "				\n"
  "	} else {				: no release occuring\n"
  "\n"
  "  	   R = R1 * exptable (- Beta * (t - (lastrelease + Cdur)))\n"
  "	}\n"
  "\n"
  "	VERBATIM\n"
  "	return 0;\n"
  "	ENDVERBATIM\n"
  "}\n"
  "\n"
  "FUNCTION exptable(x) { \n"
  "	TABLE  FROM -10 TO 10 WITH 2000\n"
  "\n"
  "	if ((x > -10) && (x < 10)) {\n"
  "		exptable = exp(x)\n"
  "	} else {\n"
  "		exptable = 0.\n"
  "	}\n"
  "}\n"
  ;
#endif
