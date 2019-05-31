#--------------------------------------------------------------------------
# Author: Sachin Talathi
#
#--------------------------------------------------------------------------

C++ = g++
RM= rm -rf
HOME= ../../../simul_lrn
HOMECC=../../../C++HH/
INCLUDE= -I$(HOME)/numlib -I$(HOME)/datalib
INCLUDE1=-I$(HOMECC)
IN=-I$(HOMECC)run_slow_compute
run_slow_compute:  Test_Slow_Syn.o rk4.o Neuron.o Synapse.o Ionchannel.o GABAAsynapse.o modelGABAAsynapse.o  modelHHneuron.o HHneuron.o modelHenrysynapse.o Henrysynapse.o Calciumchannel.o AMPAsynapse.o modelAMPAsynapse.o modelGABABsynapse.o GABABsynapse.o NMDAsynapse.o modelNMDAsynapse.o mersenne.o  Couplingcurrent.o TwoDsynapse.o modelTwoDsynapse.o
	$(C++) -g -Wall -o run_slow_compute Test_Slow_Syn.o rk4.o Neuron.o Synapse.o Ionchannel.o GABAAsynapse.o modelGABAAsynapse.o modelHHneuron.o HHneuron.o modelHenrysynapse.o Henrysynapse.o Calciumchannel.o AMPAsynapse.o modelAMPAsynapse.o modelGABABsynapse.o GABABsynapse.o NMDAsynapse.o modelNMDAsynapse.o mersenne.o  Couplingcurrent.o TwoDsynapse.o modelTwoDsynapse.o

Test_Slow_Syn.o: Test_Slow_Syn.cc 
	 $(C++)  -c  -Wall -g $(INCLUDE) Test_Slow_Syn.cc -Wno-deprecated
 
Neuron.o: $(HOMECC)Neuron.cc 
	$(C++)  -c -Wall -g $(INCLUDE) $(HOMECC)Neuron.cc -Wno-deprecated

Synapse.o: $(HOMECC)Synapse.cc 
	$(C++)  -c -Wall -g $(INCLUDE) $(HOMECC)Synapse.cc -Wno-deprecated

Ionchannel.o: $(HOMECC)/Ionchannel.cc 
	$(C++)  -c -Wall -g $(INCLUDE) $(HOMECC)/Ionchannel.cc -Wno-deprecated

HHneuron.o: $(HOMECC)HHneuron/HHneuron.cc 
	$(C++)  -c -Wall -g $(INCLUDE) $(HOMECC)HHneuron/HHneuron.cc -Wno-deprecated
	
rk4.o: $(HOMECC)rk4.cc 
	$(C++)  -c -Wall -g $(HOMECC)rk4.cc -Wno-deprecated

modelHHneuron.o: $(HOMECC)HHneuron/modelHHneuron.cc 
	$(C++)  -c -Wall -g $(HOMECC)HHneuron/modelHHneuron.cc -Wno-deprecated

GABAAsynapse.o: $(HOMECC)GABAAsynapse/GABAAsynapse.cc 
	$(C++)  -c -Wall -g $(HOMECC)GABAAsynapse/GABAAsynapse.cc -Wno-deprecated

modelGABAAsynapse.o: $(HOMECC)GABAAsynapse/modelGABAAsynapse.cc 
	$(C++)  -c -Wall -g $(HOMECC)GABAAsynapse/modelGABAAsynapse.cc -Wno-deprecated

Henrysynapse.o: $(HOMECC)Henrysynapse/Henrysynapse.cc 
	$(C++)  -c -Wall -g $(HOMECC)Henrysynapse/Henrysynapse.cc -Wno-deprecated

modelHenrysynapse.o: $(HOMECC)Henrysynapse/modelHenrysynapse.cc 
	$(C++)  -c -Wall -g $(HOMECC)Henrysynapse/modelHenrysynapse.cc -Wno-deprecated

 Calciumchannel.o: $(HOMECC)Calciumchannel.cc 
	$(C++) -c -Wall -g $(INCLUDE) $(HOMECC)Calciumchannel.cc -Wno-deprecated

AMPAsynapse.o: $(HOMECC)AMPAsynapse/AMPAsynapse.cc 
	$(C++)  -c -Wall -g $(HOMECC)AMPAsynapse/AMPAsynapse.cc -Wno-deprecated

modelAMPAsynapse.o: $(HOMECC)AMPAsynapse/modelAMPAsynapse.cc 
	$(C++)  -c -Wall -g $(HOMECC)AMPAsynapse/modelAMPAsynapse.cc -Wno-deprecated

modelGABABsynapse.o: $(HOMECC)GABABsynapse/modelGABABsynapse.cc 
	$(C++)  -c -Wall -g $(HOMECC)GABABsynapse/modelGABABsynapse.cc -Wno-deprecated

GABABsynapse.o: $(HOMECC)GABABsynapse/GABABsynapse.cc 
	$(C++)  -c -Wall -g $(HOMECC)GABABsynapse/GABABsynapse.cc -Wno-deprecated

NMDAsynapse.o: $(HOMECC)NMDAsynapse/NMDAsynapse.cc 
	$(C++)  -c -Wall -g $(HOMECC)NMDAsynapse/NMDAsynapse.cc -Wno-deprecated

modelNMDAsynapse.o: $(HOMECC)NMDAsynapse/modelNMDAsynapse.cc 
	$(C++)  -c -Wall -g $(HOMECC)NMDAsynapse/modelNMDAsynapse.cc -Wno-deprecated

mersenne.o: $(HOMECC)mersenne.cc 
	$(C++) -c -Wall -g $(HOMECC)mersenne.cc -Wno-deprecated

Couplingcurrent.o: $(HOMECC)Couplingcurrent/Couplingcurrent.cc 
	$(C++)  -c -Wall -g $(HOMECC)Couplingcurrent/Couplingcurrent.cc -Wno-deprecated

TwoDsynapse.o: $(HOMECC)TwoDsynapse/TwoDsynapse.cc 
	$(C++)  -c -Wall $(INCLUDE) $(HOMECC)TwoDsynapse/TwoDsynapse.cc -Wno-deprecated

modelTwoDsynapse.o: $(HOMECC)TwoDsynapse/modelTwoDsynapse.cc 
	$(C++)  -c -Wall $(HOMECC)TwoDsynapse/modelTwoDsynapse.cc -Wno-deprecated

 clean:
	$(RM) *.o 
