#ifndef BEHAVIOR_LEARNING_H
#define BEHAVIOR_LEARNING_H

#include <cmath>
#include <map>
#include <vector>
#include <string>
#include <iostream>
#include <sstream>
#include <fstream>
#include <iterator>
#include "nest.h"
#include "dictutils.h"
#include "nest_names.h"
#include "behavior_node.h"

namespace nest{

  
class Behavior_learning
{ 
public:
	Behavior_learning(){};
	virtual bool learn(Behavior_node *behav,DictionaryDatum &d,std::vector<double > &currents)=0;
	virtual void get(DictionaryDatum &d) const=0;
	virtual void set(const DictionaryDatum& d)=0;
	virtual ~Behavior_learning(){};
};


class Supervised_behavior: public Behavior_learning
{

public: //zamenit protected:
  std::vector<std::vector<std::vector<double> > > v_learn_set_input; //number synapse, number example, sequence spikes
  std::vector<std::vector<double> > v_learn_set_desired; //number example, sequence spikes
  int iterator_desired;
  double h;//accuracy of calculation
  int number_of_incoming_synapses;
  bool trig_learn_oper_;
  bool trig_seal_;
  long number_prim_;
  double bias_;
public:

	Supervised_behavior(const DictionaryDatum &d);
 	bool learn(Behavior_node *behav,DictionaryDatum &d,std::vector<double > &currents);
	void get(DictionaryDatum &d) const;
	void set(const DictionaryDatum& d);
	long get_number_prim();
	double get_bias();
	void read_learn_set_desired_from_file(std::string s_tmp);
	void read_learn_set_input_from_files(std::string s_directory_learn_set_input,int id_neuron);
	void reset_iterator_desired();
	bool is_desired_spike(double time,std::vector<double > &v_learn_desired); // is time desired spike == current time?
	int get_number_of_incoming_synapses();
};

class Gradient_descent_behavior: public Supervised_behavior
{ 
public:
    Gradient_descent_behavior(const DictionaryDatum &d);
    bool learn(Behavior_node *behav,DictionaryDatum &d,std::vector<double > &currents);
    void get(DictionaryDatum &d) const;
    void set(const DictionaryDatum& d);
    void fill_gradient_zero();
    double (Gradient_descent_behavior::*lambda)(double potential);
    double (Gradient_descent_behavior::*lambda_der)(double potential);
    double lambda1(double potential);
    double lambda2(double potential);
    double lambda3(double potential);
    double lambda_der1(double potential);
    double lambda_der2(double potential);
    double lambda_der3(double potential);
private:
     double k;
     double gamma;
     double Th;
     std::vector<double > gradients;
};

class STDP_learning_behavior: public Supervised_behavior
{
  public:
    STDP_learning_behavior(const DictionaryDatum &d);
    bool learn(Behavior_node *behav,DictionaryDatum &d,std::vector<double > &currents);
    void get(DictionaryDatum &d) const;
    void set(const DictionaryDatum& d);
};

/*online learning*/

class Supervised_behavior_online: public Behavior_learning
{

protected: //zamenit protected:
 
  std::vector<std::vector<double> > v_learn_set_desired; //number example, sequence spikes
  int iterator_desired;
  double h;//accuracy of calculation
  int number_of_incoming_synapses;
  bool trig_learn_oper_;
  bool trig_seal_;
  long number_prim_;
  double bias_;
public:

	Supervised_behavior_online(const DictionaryDatum &d);
 	bool learn(Behavior_node *behav,DictionaryDatum &d,std::vector<double > &currents);
	void get(DictionaryDatum &d) const;
	void set(const DictionaryDatum& d);
	long get_number_prim();
	double get_bias();
	void read_learn_set_desired_from_file(std::string s_tmp);
	void read_learn_set_input_from_files(std::string s_directory_learn_set_input,int id_neuron);//remove
	void reset_iterator_desired();
	bool is_desired_spike(double time,std::vector<double > &v_learn_desired); // is time desired spike == current time?
	int get_number_of_incoming_synapses();
};

class Gradient_descent_behavior_online: public Supervised_behavior_online
{ 
public:
    Gradient_descent_behavior_online(const DictionaryDatum &d);
    bool learn(Behavior_node *behav,DictionaryDatum &d,std::vector<double > &currents);
    void get(DictionaryDatum &d) const;
    void set(const DictionaryDatum& d);
    void fill_gradient_zero();
    double (Gradient_descent_behavior_online::*lambda)(double potential);
    double (Gradient_descent_behavior_online::*lambda_der)(double potential);
    double lambda1(double potential);
    double lambda2(double potential);
    double lambda3(double potential);
    double lambda_der1(double potential);
    double lambda_der2(double potential);
    double lambda_der3(double potential);
private:
     double k;
     double gamma;
     double Th;
     std::vector<double > gradients;
};

}


 
#endif
