#include "behavior_node.h"





void nest::Behavior_node::get_parametrs(index t_id_target,index t_id_source,double_t &t_weight,const DictionaryDatum& d)
{
  std::cout<<"Behavior_node::get_parametrs fictive function\n";
  //t_weight=val[std::pair<index,index>(t_id_target,t_id_source)].base_val.weight;
}

void nest::Iaf_neuron_sm_behavior::set_variables(index t_id_target,index t_id_source,const DictionaryDatum& d)
{
	double t_weight,t_g_max,t_E_rev;
	updateValue<double>(d, "weight",t_weight);
	updateValue<double>(d, "g_max",t_g_max);
	updateValue<double>(d, "E_rev",t_E_rev);
	struct_iaf_neuron_sm tmp(t_weight,t_g_max,t_E_rev);
	val[std::pair<index,index>(t_id_target,t_id_source)]=tmp;
}

void nest::Iaf_neuron_sm_behavior::add_spike_time(index t_id_target,index t_id_source,int local_id_synapse,double_t t_spike)
{
  val[std::pair<index,index>(t_id_target,t_id_source)].base_val.v_time_spike_times.push_back(t_spike);
}

void nest::Iaf_neuron_sm_behavior::get_parametrs(index t_id_target,index t_id_source,double_t &t_weight,const DictionaryDatum& d)
{
  //std::cout<<"Iaf_neuron_sm_behavior::get_parametrs fictive function\n";
  t_weight=val[std::pair<index,index>(t_id_target,t_id_source)].base_val.weight;
}

void nest::Iaf_neuron_sm_behavior::clear()
{
	val.clear();
}

void nest::Iaf_neuron_dif_alpha_behavior::set_variables(index t_id_target,index t_id_source,const DictionaryDatum& d)
{
	double t_weight,t_tau_a;
	double t_id_syn;
	updateValue<double>(d, "weight",t_weight);
	updateValue<double>(d, "tau_a",t_tau_a);
 	updateValue<double>(d,"local_id_synapse",t_id_syn);
	struct_iaf_neuron_dif_alpha tmp(t_weight,t_tau_a);
	std::pair<index,index> t_pair_target_source;
	t_pair_target_source=std::make_pair(t_id_target,t_id_source);
	std::pair<int,std::pair<index,index> > t_pair_syn_tar_sour;	t_pair_syn_tar_sour=std::make_pair(int(t_id_syn),t_pair_target_source);
	val[t_pair_syn_tar_sour]=tmp;
	//tmp
}

void nest::Iaf_neuron_dif_alpha_behavior::add_spike_time(index t_id_target,index t_id_source,int local_id_synapse,double_t t_spike)
{
  std::pair<index,index> t_pair_target_source;
  t_pair_target_source=std::make_pair(t_id_target,t_id_source);
  std::pair<int,std::pair<index,index> > t_pair_syn_tar_sour;
  t_pair_syn_tar_sour=std::make_pair(local_id_synapse,t_pair_target_source);
  val[t_pair_syn_tar_sour].base_val.v_time_spike_times.push_back(t_spike);
}

void nest::Iaf_neuron_dif_alpha_behavior::get_parametrs(index t_id_target,index t_id_source,double_t &t_weight,const DictionaryDatum& d)
{
  //std::cout<<"Iaf_neuron_sm_behavior::get_parametrs fictive function\n";
  double t_id_syn;
  updateValue<double>(d,"local_id_synapse",t_id_syn);
  std::pair<index,index> t_pair_target_source;
  t_pair_target_source=std::make_pair(t_id_target,t_id_source);
  std::pair<int,std::pair<index,index> > t_pair_syn_tar_sour;	t_pair_syn_tar_sour=std::make_pair(int(t_id_syn),t_pair_target_source);
  t_weight=val[t_pair_syn_tar_sour].base_val.weight;
}

void nest::Iaf_neuron_dif_alpha_behavior::clear()
{
	val.clear();
}

/*Iaf_neuron_dif_alpha_stdp_behavior*/

void nest::Iaf_neuron_dif_alpha_stdp_behavior::set_variables(index t_id_target,index t_id_source,const DictionaryDatum& d)
{
	//struct_iaf_neuron_dif_alpha_stdp tmp(t_weight,t_tau_a,20.0,1.0,100.0); //no_good
	double t_weight,t_weight_dynamic;
	updateValue<double>(d, "weight",t_weight);
	updateValue<double>(d, "dynamic_weight",t_weight_dynamic);
	val[std::pair<index,index>(t_id_target,t_id_source)].base_val.weight=t_weight;
	//val[std::pair<index,index>(t_id_target,t_id_source)].dynamic_weight=t_weight_dynamic;// wrong
	val[std::pair<index,index>(t_id_target,t_id_source)].dynamic_weights.push_back(t_weight_dynamic);
	//tmp
}

void nest::Iaf_neuron_dif_alpha_stdp_behavior::add_spike_time(index t_id_target,index t_id_source,int local_id_synapse,double_t t_spike)
{
  val[std::pair<index,index>(t_id_target,t_id_source)].base_val.v_time_spike_times.push_back(t_spike);
}

void nest::Iaf_neuron_dif_alpha_stdp_behavior::get_parametrs(index t_id_target,index t_id_source,double_t &t_weight,const DictionaryDatum& d)
{
  //std::cout<<"Iaf_neuron_sm_behavior::get_parametrs fictive function\n";
  t_weight=val[std::pair<index,index>(t_id_target,t_id_source)].base_val.weight;
}

void nest::Iaf_neuron_dif_alpha_stdp_behavior::clear()
{
	val.clear();
}

void nest::Iaf_neuron_dif_alpha_stdp_behavior::set(index t_id_target,index t_id_source,const DictionaryDatum& d)
{
  double weight_,tau_a_,lambda_,tau_plus_,Wmax_;
  updateValue<double>(d, "weight",weight_);
  updateValue<double>(d, "tau_a",tau_a_);
  updateValue<double>(d, "lambda",lambda_);
  updateValue<double>(d, "tau_plus",tau_plus_);
  updateValue<double>(d, "Wmax",Wmax_);
  struct_iaf_neuron_dif_alpha_stdp tmp(weight_,tau_a_,tau_plus_,lambda_,Wmax_);//no_good
  val[std::pair<index,index>(t_id_target,t_id_source)]=tmp;
}