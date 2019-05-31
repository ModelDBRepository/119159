//--------------------------------------------------------------------------
// Author: Thomas Nowotny
//
// Institute: Institute for Nonlinear Dynamics
//            University of California San Diego
//            La Jolla, CA 92093-0402
//
// email to:  tnowotny@ucsd.edu
//
// initial version: 2002-01-24
//
//--------------------------------------------------------------------------

//--------------------------------------------------------------------------
// adapted from LPneuronJorge with
// a) own fit for Ih
// b) GHK formalism for ICa + own fit with it
//--------------------------------------------------------------------------

#ifndef LPRNEURON_H
#define LPRNEURON_H

#include "CN_neuron.h"
#include <cmath>

#define LPR_PNO 73
#define LPR_IVARNO 16

//--------------------------------------------------------------------------
// parameters of the LPR neuron (typical set, can be modified at any time
// later on. Modifications take effect immediately in the next evaluation)

double stdLPR_p[LPR_PNO]= {
   80,                  // 0 - gNa: Na conductance [muS]
   47.1843,             // 1 - ENa: Na equi potential in mV
   0.247102,            // 2 - gCa1: first Ca conductance
   0.0272075,           // 3 - gCa2: second Ca conductance
   2.87471,             // 4 - goCa: KCa conductance
   -83.8393,            // 5 - EK: K rev. potential
   0.401598,            // 6 - gd: Kd conductance in [muS]
   2.26681,             // 7 - gA: A channel conductance
   0.0448226,           // 8 - gh: Ih channel conductance
   -11.6324,            // 9 - EIh: Ih channel rev. potential
   0.0100637,           // 10 - gleak: leak conductance neuropil
   -50.8427,            // 11 - Eleak: leak reversal potential
   0.0502388,           // 12 - Cmem: membrane capacitance neuropil
   0.0449852,           // 13 - c_iCa: factor in Ca dynamics [mM/pCb]
   0.0446877,           // 14 - k_Ca: time scale Ca removal [1/ms]
   3.44474,             // 15 - tfac: rescaling all taum/tauh
   2.70653,             // 16 - Vshift: shift of V to outside world
   0.372713,            // 17 - IDC: DC input current to set dyn. regime
   1.05799,             // 18 - tfac_fast: recaling fast taum/tauh (Na, Kd)
   0.0492504,           // 19 - gVV: conductance between neuropil & soma
   0.693825,            // 20 - Cs: capacitance of soma
   0.0214249,           // 21 - gleaks: leak conduct. soma
   11.4866,             // 22 - RTF
   1.0,                 // 23 - AreaAxon
   1.0,                 // 24 - Area
   0.045,               // 25 - k_aCa1
   0.02,                // 26 - k_bCa1
   0.0608644,           // 27 - k_aCa2
   -5.7,                // 28 - V_aCa1
   -38,                 // 29 - V_bCa1
   29.1336,             // 30 - V_aCa2
   -7.61858,            // 31 - s_aCa1
   1.2,                 // 32 - s_bCa1
   -4.43083,            // 33 - s_aCa2
   0.00155982,          // 34 - P_Ca
   15123.5,             // 35 - Ca_out
   -15,                 // 36 - V_kbCa1
   -3.8,                // 37 - s_kbCa1
   0.480858,            // 38 - k_oa
   0.0374617,           // 39 - k_ob
   0.17997,             // 40 - V_ao1
   -19.6809,            // 41 - V_ao2
   -23.3747,            // 42 - s_ao1
   -5.01798,            // 43 - s_ao2
   0.606697,            // 44 - f
   2.54667,             // 45 - c1
   0.710174,            // 46 - c2
   0.598139,            // 47 - c3
   0.018541,            // 48 - Ca_0
   0.142092,            // 49 - k_aA
   0.0475747,           // 50 - k_bA1
   0.00399543,          // 51 - c_A2
   -11.5622,            // 52 - V_aA
   -60.4254,            // 53 - V_bA
   -37.9344,            // 54 - V_kbA2
   5.14585,             // 55 - V_x
   -25.8765,            // 56 - s_aA
   5.4265,              // 57 - s_bA
   -12.6589,            // 58 - s_kbA2
   -16.7579,            // 59 - s_x
   2.10007e-05,         // 60 - c_r
   -80.5426,            // 61 - V_r
   -64.2891,            // 62 - V_kr
   7.98176,             // 63 - s_r
   -21.5582,            // 64 - s_kr
   0.3,                 // 65 - I_scale
   -93.6857,            // 66 - E_M
   0.128435,            // 67 - g_M
   7.66637e-05,         // 68 - k_M
   -20.0768,            // 69 - V_M
   -4.67792,            // 70 - s_M
   -60.0,               // 71 - V_kM
   4.0,                 // 72 - s_kM
};

double *LPR_p= stdLPR_p;

char *LPR_p_text[LPR_PNO]= {
  "0 - gNa: Na conductance [muS]",
  "1 - ENa: Na equi potential in mV",
  "2 - gCa1: first Ca conductance",
  "3 - gCa2: second Ca conductance",
  "4 - goCa: KCa conductance",
  "5 - EK: K rev. potential",
  "6 - gd: Kd conductance in [muS]",
  "7 - gA: A channel conductance",
  "8 - gh: Ih channel conductance",
  "9 - EIh: Ih channel rev. potential",
  "10 - gleak: leak conductance neuropil",
  "11 - Eleak: leak reversal potential",
  "12 - Cmem: membrane capacitance neuropil",
  "13 - c_iCa: factor in Ca dynamics [mM/pCb]",
  "14 - k_Ca: time scale Ca removal [1/ms]",
  "15 - tfac: rescaling all taum/tauh",
  "16 - Vshift: shift of V to outside world",
  "17 - IDC: DC input current to set dyn. regime",
  "18 - tfac_fast: recaling fast taum/tauh (Na, Kd)",
  "19 - gVV: conductance between neuropil & soma",
  "20 - Cs: capacitance of soma",
  "21 - gleaks: leak conduct. soma",
  "22 - RTF",
  "23 - V_m",
  "24 - s_m",
  "25 - k_aCa1",
  "26 - k_bCa1",
  "27 - k_aCa2",
  "28 - V_aCa1",
  "29 - V_bCa1",
  "30 - V_aCa2",
  "31 - s_aCa1",
  "32 - s_bCa1",
  "33 - s_aCa2",
  "34 - P_Ca",
  "35 - Ca_out",
  "36 - V_kbCa1",
  "37 - s_kbCa1",
  "38 - k_oa",
  "39 - k_ob",
  "40 - V_ao1",
  "41 - V_ao2",
  "42 - s_ao1",
  "43 - s_ao2",
  "44 - f",
  "45 - c1",
  "46 - c2",
  "47 - c3",
  "48 - Ca_0",
  "49 - k_aA",
  "50 - k_bA1",
  "51 - c_A2",
  "52 - V_aA",
  "53 - V_bA",
  "54 - V_kbA2",
  "55 - V_x",
  "56 - s_aA",
  "57 - s_bA",
  "58 - s_kbA2",
  "59 - s_x",
  "60 - c_r",
  "61 - V_r",
  "62 - V_kr",
  "63 - s_r",
  "64 - s_kr",
  "65 - I_scale",
  "66 - E_M",
  "67 - g_M",
  "68 - k_M",
  "69 - V_M",
  "70 - s_M",
  "71 - V_kM",
  "72 - s_kM"
};


//----------------------------------------------------------------------
// the following initial values reflect the steady state with no input

double stdLPR_INIVARS[LPR_IVARNO]= {
  -56.0,                       // 0 - membrane potential 
  0.0,                         // 1 - Na channel activation mNa
  0.2,                         // 2 - Na channel unblocking hNa
  0.0,                         // 3 - Ca1 channel activation mCa1
  1.0,                         // 4 - Ca1 channel unblocking hCa1
  0.0,                         // 5 - Ca2 channel activation mCa2
  0.0,                         // 6 - oCa channel activation moCa
  1.0,                         // 7 - oCa channel unblocking hoCa
  0.0,                         // 8 - d channel activation md
  0.0,                         // 9 - A channel activation mA
  1.0,                         // 10 - A channel unblocking hA1
  1.0,                         // 11 - A channel unblocking hA2
  0.0,                         // 12 - Ih channel activation
  0.04,                        // 13 - Ca concentration
  0.01,                        // 14 - M current activation var
  -56.0                        // 15 - membrane potential soma
};
double *LPR_INIVARS= stdLPR_INIVARS;

char *LPR_INIVARSTEXT[LPR_IVARNO]= {
  "0 - membrane potential",
  "1 - Na channel activation mNa",
  "2 - Na channel unblocking hNa",
  "3 - Ca1 channel activation mCa1",
  "4 - Ca1 channel unblocking hCa1",
  "5 - Ca2 channel activation mCa2",
  "6 - oCa channel activation moCa",
  "7 - oCa channel unblocking hoCa",
  "8 - d channel activation md",
  "9 - A channel activation mA",
  "10 - A channel unblocking hA1",
  "11 - A channel unblocking hA2",
  "12 - Ih channel activation",
  "13 - Ca concentration",
  "14 - M current activation var",
  "15 - membrane potential soma"
};


// the LPR neuron class itself

class LPRneuron: public neuron
{
 public:
  LPRneuron(int, double *);
  LPRneuron(int, tnvector<int>, double *);
  virtual ~LPRneuron() { }
  inline virtual double E(double *);
  virtual void currents(ostream &, double *);
  virtual void derivative(double *, double *);
};

// for fitting purposes: type == 0 means additive changes
//                       type == 1 means multiplicative changes

double LPR_p_type[LPR_PNO]= {
  1,      // 0 - gNa: Na conductance [muS] paper
  0,      // 1 - ENa: Na equi potential in mV
  1,      // 2 - gCa1: first Ca conductance
  1,      // 3 - gCa2: second Ca conductance
  1,      // 4 - goCa: KCa conductance
  0,      // 5 - EK: K rev. potential
  1,      // 6 - gd: Kd conductance in [muS]
  1,      // 7 - gA: A channel conductance
  1,      // 8 - gh: Ih channel conductance
  0,      // 9 - EIh: Ih channel rev. potential
  1,      // 10 - gleak: leak conductance neuropil
  0,      // 11 - Eleak: leak reversal potential
  1,      // 12 - Cmem: membrane capacitance neuropil
  1,      // 13 - c_iCa: factor in Ca dynamics [mM/pCb]
  1,      // 14 - k_Ca: time scale Ca removal [1/ms]
  1,      // 15 - tfac: rescaling all taum/tauh
  0,      // 16 - Vshift: shift of V to outside world
  0,      // 17 - IDC: DC input current to set dyn. regime
  1,      // 18 - tfac_fast: recaling fast taum/tauh (Na, Kd)
  1,      // 19 - gVV: conductance between neuropil & soma
  1,      // 20 - Cs: capacitance of soma
  1,      // 21 - gleaks: leak conduct. soma
  1,      // 22 - RTF
  1,      // 23 - AreaAxon
  1,      // 24 - Area
  1,      // 25 - k_aCa1
  1,      // 26 - k_bCa1
  1,      // 27 - k_aCa2
  0,      // 28 - V_aCa1
  0,      // 29 - V_bCa1
  0,      // 30 - V_aCa2
  1,      // 31 - s_aCa1
  1,      // 32 - s_bCa1
  1,      // 33 - s_aCa2
  1,      // 34 - P_Ca
  1,      // 35 - Ca_out
  0,      // 36 - V_kbCa1
  1,      // 37 - s_kbCa1
  1,      // 38 - k_oa
  1,      // 39 - k_ob
  0,      // 40 - V_ao1
  0,      // 41 - V_ao2
  1,      // 42 - s_ao1
  1,      // 43 - s_ao2
  1,      // 44 - f
  1,      // 45 - c1
  1,      // 46 - c2
  1,      // 47 - c3
  1,      // 48 - Ca_0
  1,      // 49 - k_aA
  1,      // 50 - k_bA1
  1,      // 51 - c_A2
  0,      // 52 - V_aA
  0,      // 53 - V_bA
  0,      // 54 - V_kbA2
  0,      // 55 - V_x
  1,      // 56 - s_aA
  1,      // 57 - s_bA
  1,      // 58 - s_kbA2
  1,      // 59 - s_x
  1,      // 60 - c_r
  0,      // 61 - V_r
  0,      // 62 - V_kr
  1,      // 63 - s_r
  1,      // 64 - s_kr
  1,      // 65 - I_scale
  0,      // 66 - E_M
  1,      // 67 - g_M
  1,      // 68 - k_M
  0,      // 69 - V_M
  1,      // 70 - s_M
  0,      // 71 - V_kM
  1       // 72 - s_kM
};

// for fitting purposes ... the LPR_p_sens(itivity) is basically the inverse
// of the slope of the deviation of the neuron model perturbed in a given
// parameter compared to the model with the original parameter set.
// see "test_params" and "calc_sensitivity" for more details

double LPR_p_sens[LPR_PNO]= {
  5e-10, // 0 - gNa: Na conductance [muS] paper
  1e-08, // 1 - ENa: Na equi potential in mV
  5e-10, // 2 - gCa1: first Ca conductance
  5e-10, // 3 - gCa2: second Ca conductance
  5e-10, // 4 - goCa: KCa conductance
  1e-08, // 5 - EK: K rev. potential
  5e-10, // 6 - gd: Kd conductance in [muS]
  5e-10, // 7 - gA: A channel conductance
  1e-9,  // 8 - gh: Ih channel conductance
  1-07,  // 9 - EIh: Ih channel rev. potential
  1e-09, // 10 - gleak: leak conductance neuropil
  1e-09, // 11 - Eleak: leak reversal potential
  1e-09, // 12 - Cmem: membrane capacitance neuropil
  1e-09, // 13 - c_iCa: factor in Ca dynamics [mM/pCb]
  1e-09, // 14 - k_Ca: time scale Ca removal [1/ms]
  5e-10, // 15 - tfac: rescaling all taum/tauh
  1e-8,  // 16 - Vshift: shift of V to outside world
  5e-9,  // 17 - IDC: DC input current to set dyn. regime
  5e-10, // 18 - tfac_fast: recaling fast taum/tauh (Na, Kd)
  1e-9,  // 19 - gVV: conductance between neuropil & soma
  1e-9,  // 20 - Cs: capacitance of soma
  1e-9,  // 21 - gleaks: leak conduct. soma
  1e-09, // 22 - RTF
  1e-9,  // 23 - AreaAxon
  1e-9,  // 24 - Area
  1e-9,  // 25 - k_aCa1
  1e-09, // 26 - k_bCa1
  1e-08, // 27 - k_aCa2
  1e-09, // 28 - V_aCa1
  1e-09, // 29 - V_bCa1
  1e-07, // 30 - V_aCa2
  1e-11, // 31 - s_aCa1
  1e-09, // 32 - s_bCa1
  1e-08, // 33 - s_aCa2
  5e-10, // 34 - P_Ca
  5e-10, // 35 - Ca_out
  1e-08, // 36 - V_kbCa1
  5e-10, // 37 - s_kbCa1
  1e-08, // 38 - k_oa
  1e-08, // 39 - k_ob
  5e-08, // 40 - V_ao1
  1e-08, // 41 - V_ao2
  1e-09, // 42 - s_ao1
  5e-10, // 43 - s_ao2
  1e-08, // 44 - f
  1e-09, // 45 - c1
  1e-09, // 46 - c2
  1e-09, // 47 - c3
  5e-9,  // 48 - Ca_0
  1e-08, // 49 - k_aA
  1e-07, // 50 - k_bA1
  1e-08, // 51 - c_A2
  1e-08, // 52 - V_aA
  1e-08, // 53 - V_bA
  1e-07, // 54 - V_kbA2
  1e-06, // 55 - V_x
  5e-10, // 56 - s_aA
  1e-9,  // 57 - s_bA
  1e-08, // 58 - s_kbA2
  1e-07, // 59 - s_x
  5e-08, // 60 - c_r
  1e-07, // 61 - V_r
  1e-06, // 62 - V_kr
  1e-08, // 63 - s_r
  1e-07, // 64 - s_kr
  1e-9,  // 65 - I_scale
  1e-07, // 66 - E_M
  1e-9,  // 67 - g_M
  1e-9,  // 68 - k_M
  1e-8,  // 69 - V_M
  1e-9,  // 70 - s_M
  1e-8,  // 71 - V_kM
  1e-9,  // 72 - s_kM
};


#endif
