#ifndef BEHAVIOR_NODE_H
#define BEHAVIOR_NODE_H

#include <map>
#include <vector>
#include <string>
#include <iostream>
#include "nest.h"
#include "dictutils.h"
#include "nest_names.h"

namespace nest{

/*struct struct_behavior_node
{
  index id_target;
  index id_source;
  double weight;
  struct_behavior_node()  {};
  struct_behavior_node(index t_id_target,index t_id_source,double t_weight):
  id_target(t_id_target),id_source(t_id_source),weight(t_weight)
  {}
};*/

struct struct_behavior_node
{
  std::vector<double_t > v_time_spike_times;
  double weight;
  struct_behavior_node()  {};
  struct_behavior_node(double t_weight):
  weight(t_weight)
  {}
};
  
class Behavior_node
{ 
public:
	bool trig_fict_spike;//for fictive spike
	Behavior_node(){trig_fict_spike=false;};
	/*чушь, но по-другому неизвестно как*/
	//virtual void set_variables()=0;
	virtual void set_variables(index id_target,index id_source,double_t t_g_max,double_t t_weight,double_t t_E_rev)=0;
	virtual void add_spike_time(index id_target,index id_source,double_t t_spike)=0;
	virtual void get_parametrs(index t_id_target,index t_id_source,double_t &t_weight);
	virtual void clear()=0;
	virtual void set(index t_id_target,index t_id_source,const DictionaryDatum& d)=0;
	virtual ~Behavior_node(){};
};

/*struct struct_iaf_neuron_sm
{
  struct_behavior_node base_val;
  double g_max;
  double E_rev;

  struct_iaf_neuron_sm(){};
  
  struct_iaf_neuron_sm(index t_id_target,index t_id_source,double t_weight,double t_g_max,double t_E_rev):
  base_val(t_id_target,t_id_source,t_weight),
  g_max(t_g_max),
  E_rev(t_E_rev)
  {}
};*/

struct struct_iaf_neuron_sm
{
  struct_behavior_node base_val;
  double g_max;
  double E_rev;
  //double g_grad; // gradient

  struct_iaf_neuron_sm(){};
  
  struct_iaf_neuron_sm(double t_weight,double t_g_max,double t_E_rev):
  base_val(t_weight),
  g_max(t_g_max),
  E_rev(t_E_rev)
  {}
};

class Iaf_neuron_sm_behavior: public Behavior_node
{ 
public:

	Iaf_neuron_sm_behavior():Behavior_node(){};
	//void set_variables(){};
	void set_variables(index id_target,index id_source,double_t t_g_max,double_t t_weight,double_t t_E_rev);
	void add_spike_time(index t_id_target,index t_id_source,double_t t_spike);
	void get_parametrs(index t_id_target,index t_id_source,double_t &t_weight);
	void clear();
	void set(index t_id_target,index t_id_source,const DictionaryDatum& d){};
	//std::vector<struct_iaf_neuron_sm> val;
	std::map<std::pair<index,index>,struct_iaf_neuron_sm> val;
};

struct struct_iaf_neuron_dif_alpha
{
  struct_behavior_node base_val;
  double tau_a;

  struct_iaf_neuron_dif_alpha(){};
  
  struct_iaf_neuron_dif_alpha(double t_weight,double t_tau_a):
  base_val(t_weight),
  tau_a(t_tau_a)
  {}
};

class Iaf_neuron_dif_alpha_behavior: public Behavior_node
{ 
public:

	Iaf_neuron_dif_alpha_behavior():Behavior_node(){};
	
	void set_variables(index id_target,index id_source,double_t t_tau_a,double_t t_weight,double_t t_E_rev);
	void add_spike_time(index t_id_target,index t_id_source,double_t t_spike);
	void get_parametrs(index t_id_target,index t_id_source,double_t &t_weight);
	void clear();
	void set(index t_id_target,index t_id_source,const DictionaryDatum& d){};
	//std::vector<struct_iaf_neuron_sm> val;
	std::map<std::pair<index,index>,struct_iaf_neuron_dif_alpha> val;
};

struct struct_iaf_neuron_dif_alpha_stdp
{
  struct_behavior_node base_val;
  double tau_a;
  double tau_plus;
  double lambda;
  double Wmax;
  double dynamic_weight;
  struct_iaf_neuron_dif_alpha_stdp(){};
  
  struct_iaf_neuron_dif_alpha_stdp(double t_weight,double t_tau_a,double t_tau_plus,double t_lambda,double t_Wmax,double t_dynamic_weight):
  base_val(t_weight),
  tau_a(t_tau_a),
  tau_plus(t_tau_plus),
  lambda(t_lambda),
  Wmax(t_Wmax),
  dynamic_weight(t_dynamic_weight)
  {}
};

class Iaf_neuron_dif_alpha_stdp_behavior: public Behavior_node
{ 
public:

	Iaf_neuron_dif_alpha_stdp_behavior():Behavior_node(){};
	
	void set_variables(index id_target,index id_source,double_t t_tau_a,double_t t_weight,double_t t_E_rev);
	void add_spike_time(index t_id_target,index t_id_source,double_t t_spike);
	void get_parametrs(index t_id_target,index t_id_source,double_t &t_weight);
	void clear();
	void set(index t_id_target,index t_id_source,const DictionaryDatum& d);
	//std::vector<struct_iaf_neuron_sm> val;
	std::map<std::pair<index,index>,struct_iaf_neuron_dif_alpha_stdp> val;
};

}


 
#endif