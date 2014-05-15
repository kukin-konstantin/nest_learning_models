#include "behavior_learning.h"




nest::Supervised_behavior::Supervised_behavior(const DictionaryDatum &d)
:Behavior_learning(),iterator_desired(0)
{
  updateValue<double>(d, names::h, h);
  //updateValue<bool>(d, names::trig_learn_oper, trig_learn_oper_);
  trig_learn_oper_=false;
  //updateValue<bool>(d, names::trig_seal, trig_seal_);
  trig_seal_=false;
  //updateValue<int>(d, names::number_prim, number_prim_);
  number_prim_=0;
  //updateValue<double>(d, names::bias, bias_);
  bias_=0;
}

void nest::Supervised_behavior::read_learn_set_desired_from_file(std::string s_tmp)
{
  std::vector<double> tmp;
  std::istringstream iss;
  std::ifstream f(s_tmp.c_str());
  std::string line;
  if (f.is_open())
  {
    while ( f.good() )
    {
	getline (f,line);
	iss.str(line);
	std::copy (std::istream_iterator<double>(iss), std::istream_iterator<double>(), std::back_inserter(tmp));
	//std::cout<<"read_learn inside"<<"\n";
	v_learn_set_desired.push_back(tmp);
	tmp.clear();
	iss.clear();
    }
  }
  //std::cout<<"read_learn "<<"\n";
  f.clear();
  f.close();
}


void nest::Supervised_behavior::read_learn_set_input_from_files(std::string s_directory_learn_set_input,int id_neuron)
{
  int_t n_incom;
  std::string s_tmp;
  std::stringstream r_number_neuron;
  r_number_neuron << id_neuron;
  s_tmp=s_directory_learn_set_input+r_number_neuron.str()+"_number.txt";
  std::ifstream f_number(s_tmp.c_str());
  if (f_number.good())
  {
    f_number>>n_incom;
    number_of_incoming_synapses=n_incom;
    f_number.clear();
    f_number.close();
    v_learn_set_input.resize(n_incom);
    std::cout<<"n_incom "<<n_incom<<"\n";
    for (int i=0;i!=n_incom;i++)
    {
      std::stringstream r_number_synapse;
      r_number_synapse << i+1;
      s_tmp=s_directory_learn_set_input+r_number_neuron.str()+"_k_"+r_number_synapse.str()+".txt";
      std::cout<<"s_tmp "<<s_tmp<<"\n";
      std::vector<double> tmp;
      std::istringstream iss;
      std::ifstream f(s_tmp.c_str());
      std::string line;
      if (f.is_open())
      {
	while ( f.good() )
	{
	  getline (f,line);
	  iss.str(line);
	  std::copy (std::istream_iterator<double>(iss), std::istream_iterator<double>(), std::back_inserter(tmp));
	  //std::cout<<"read_learn inside"<<"\n";
	  v_learn_set_input[i].push_back(tmp);
	  tmp.clear();
	  iss.clear();
	}
      }
      //std::cout<<"read_learn "<<"\n";
      f.clear();
      f.close();
    }
  }
  else
  {
    std::cout<<"bad s_tmp "<<s_tmp<<"\n";
  }
}

void nest::Supervised_behavior::reset_iterator_desired()
{
  iterator_desired=0;
}

int nest::Supervised_behavior::get_number_of_incoming_synapses()
{
  return number_of_incoming_synapses;
}

bool nest::Supervised_behavior::is_desired_spike(double time,std::vector<double > &v_learn_desired)
{
  bool is_spike;
  //std::cout<<"iterator_desired "<<iterator_desired<<"\t";
  //std::cout<<"time_desired "<<v_learn_desired[iterator_desired]<<"\t";
  //std::cout<<time<<"\t";
  double time_spike_desired=v_learn_desired[iterator_desired];
  if ((fabs(time_spike_desired-time)<=(h/2.0))&&(time_spike_desired>0))
  {
    is_spike=true;
    //std::cout<<"popalo";
    if (v_learn_desired.size()!=iterator_desired+1)
      iterator_desired++;
  }
  else
  {
    is_spike=false;
  }
  //std::cout<<"exit is desired\n";
  return is_spike;
}

bool nest::Supervised_behavior::learn(Behavior_node *behav,DictionaryDatum &d,std::vector<double > &currents)
{
  return false;
}

void nest::Supervised_behavior::get(DictionaryDatum &d) const
{
  def<bool>(d, names::trig_learn_oper,trig_learn_oper_);
  def<bool>(d, names::trig_seal,trig_seal_);
  def<long>(d, names::number_prim,number_prim_);
  def<double>(d, names::bias,bias_);
  def<double>(d, names::h,h);
}

void nest::Supervised_behavior::set(const DictionaryDatum& d)
{
  updateValue<bool>(d, names::trig_learn_oper,trig_learn_oper_);
  updateValue<bool>(d, names::trig_seal,trig_seal_);
  updateValue<long>(d, names::number_prim,number_prim_);
  updateValue<double>(d, names::bias,bias_);
  updateValue<double>(d, names::h,h);
}

long nest::Supervised_behavior::get_number_prim()
{
  return number_prim_;
}

double nest::Supervised_behavior::get_bias()
{
  return bias_;
}

nest::Gradient_descent_behavior::Gradient_descent_behavior(const DictionaryDatum &d)
:Supervised_behavior(d)
{
    //change number of lambda
    //lambda=&Gradient_descent_behavior::lambda1;
    //lambda_der=&Gradient_descent_behavior::lambda_der1;
    lambda=&Gradient_descent_behavior::lambda3;
    lambda_der=&Gradient_descent_behavior::lambda_der3;
    updateValue<double>(d, names::Th, Th);
    updateValue<double>(d, names::k, k);
    updateValue<double>(d, names::gamma, gamma);
    std::cout<<" gamma in gradient "<<gamma<<"\n";
}

double nest::Gradient_descent_behavior::lambda1(double potential)
{
  return std::exp((potential-Th)/k);
}

double nest::Gradient_descent_behavior::lambda2(double potential)
{
  if (potential>=0.0)
  {
    return potential;
  }
  else
  {
    return 0.0;
  }
}

double nest::Gradient_descent_behavior::lambda3(double potential)
{
  return 1.0/(1.0+std::exp(-(potential-Th)/k));
}

double nest::Gradient_descent_behavior::lambda_der1(double potential)
{
  return (1.0/k)*lambda1(potential);
}

double nest::Gradient_descent_behavior::lambda_der2(double potential)
{
  if (potential>=0.0)
  {
    return 1.0;
  }
  else
  {
    return 0.0;
  }
}

double nest::Gradient_descent_behavior::lambda_der3(double potential)
{
  double sigm=lambda3(potential);
  return (1.0/k)*sigm*(1-sigm);
}

bool nest::Gradient_descent_behavior::learn(Behavior_node *behav,DictionaryDatum &d,std::vector<double > &currents)
{
  //Iaf_neuron_sm_behavior* behav_loc= dynamic_cast<Iaf_neuron_sm_behavior *>(behav);
  //std::cout<<"typeid behav="<<typeid(behav).name()<<"\t"<<behav<<"\n";
  Iaf_neuron_dif_alpha_behavior* behav_loc= dynamic_cast<Iaf_neuron_dif_alpha_behavior *>(behav);
  //tr/*Iaf_neuron_dif_alpha_behavior* behav_loc= new(Iaf_neuron_dif_alpha_behavior);
  //tr/* *behav_loc= *dynamic_cast<Iaf_neuron_dif_alpha_behavior *>(behav);
  //std::cout<<"typeif behav_loc="<<typeid(behav_loc).name()<<"\t"<<behav_loc<<"\n";
  //std::cout<<"typeif behav_second="<<typeid(behav).name()<<"\t"<<behav<<"\n";
  if (gradients.empty())
  {
    gradients.resize(currents.size());
    std::fill(gradients.begin(),gradients.end(),0.0);
  }
  double time;
  double u;
  updateValue<double>(d, names::time, time);
  updateValue<double>((*d), names::potential,u);
  double lambd_val=(this->*lambda)(u);
  double lambd_val_der=(this->*lambda_der)(u);
  double L=1.0-std::exp(-lambd_val*h);
  //*std::cout<<"u="<<u<<"\t";
  //*std::cout<<"lambd_val="<<lambd_val<<"\t";
  //*std::cout<<"lambd_val_der="<<lambd_val_der<<"\t";
  //*std::cout<<"std::exp "<<std::exp(-lambd_val*h)<<"\t";
  //*std::cout<<"L="<<L<<"\t";
  double learn_time=time-bias_;
  //std::cout<<"learn time s"<<learn_time<<"\n";
  //*std::cout<<"number_prim_ "<<number_prim_<<"\n";
  std::vector<double > v_learn_desired=v_learn_set_desired[number_prim_];
  bool is_spike=is_desired_spike(learn_time,v_learn_desired);
  //std::cout<<"is_spike "<<is_spike<<"\n";
  if (is_spike) //has spike
  {
    //std::cout<<"is spike "<<"\n";
    //std::map<std::pair<index,index>,struct_iaf_neuron_sm >::iterator it=behav_loc->val.begin();
     std::map<std::pair<int,std::pair<index,index> >,struct_iaf_neuron_dif_alpha >::iterator it=behav_loc->val.begin();
    int sh=0;
     
    while (it!=behav_loc->val.end())
    {
      if (L<1e-7)//magical number
      {
	L=1e-7;
	gradients[sh]-=(lambd_val_der/L)*currents[sh];
      }
      else
      {
	gradients[sh]-=((1.0-L)/L)*lambd_val_der*currents[sh];
      }
      //*std::cout<< gradients[sh]<<"("<<currents[sh]<<")"<<"("<<(*it).second.base_val.weight;
      //std::cout<<"bef_weig "<<(*it).second.base_val.weight<<"\t";
      //std::cout<<"s "<<gamma<<"\t"<<h<<"\t";
      (*it).second.base_val.weight-=gamma*gradients[sh]*h;
      //std::cout<<"after_weig "<<(*it).second.base_val.weight<<"\t";
      //*std::cout<<"->"<<(*it).second.base_val.weight<<")"<<"\t";
      gradients[sh]=0;
      it++;
      sh++;
    }
     
  }
  else
  {
    //*std::cout<<"no spike ";
    for (int i=0;i!=currents.size();i++)
    {
      gradients[i]+=lambd_val_der*currents[i];
      //*std::cout<< gradients[i]<<"("<<currents[i]<<")"<<"\t";
    }
  }
  //*std::cout<<"\n";
  //std::cout<<"is_spike "<<is_spike<<"\n";
  //tr/*delete behav_loc;
  //std::cout<<"behav="<<behav<<"\n";
  //std::cout<<"behav_loc="<<behav_loc<<"\n";
  return is_spike;
}

void nest::Gradient_descent_behavior::get(DictionaryDatum &d) const
{
  Supervised_behavior::get(d);
  def<double>(d, names::k,k);
  def<double>(d, names::gamma,gamma);
  def<double>(d, names::Th,Th);
  def<double>(d, names::bias,bias_); //?
}

void nest::Gradient_descent_behavior::set(const DictionaryDatum& d)
{
  Supervised_behavior::set(d);
  updateValue<double>(d, names::k,k);
  updateValue<double>(d, names::gamma,gamma);
  updateValue<double>(d, names::Th,Th);
  updateValue<double>(d, names::bias,bias_); //?
}

void nest::Gradient_descent_behavior::fill_gradient_zero()
{
  if (!gradients.empty())
  {
    std::fill(gradients.begin(),gradients.end(),0.0);
  }
}

nest::STDP_learning_behavior::STDP_learning_behavior(const DictionaryDatum &d)
:Supervised_behavior(d)
{

}

bool nest::STDP_learning_behavior::learn(Behavior_node *behav,DictionaryDatum &d,std::vector<double > &currents)
{
  //Iaf_neuron_sm_behavior* behav_loc= dynamic_cast<Iaf_neuron_sm_behavior *>(behav);
  //Iaf_neuron_dif_alpha_behavior* behav_loc= dynamic_cast<Iaf_neuron_dif_alpha_behavior *>(behav);
  Iaf_neuron_dif_alpha_behavior* behav_loc= dynamic_cast<Iaf_neuron_dif_alpha_behavior *>(behav);
  double time;
  updateValue<double>(d, names::time, time);
  //std::cout<<"time "<<time<<"\n";
  double learn_time=time-bias_;
  std::vector<double > v_learn_desired=v_learn_set_desired[number_prim_];
  bool is_spike=is_desired_spike(learn_time,v_learn_desired);
  delete behav_loc;
  return is_spike;
}

void nest::STDP_learning_behavior::get(DictionaryDatum &d) const
{
  Supervised_behavior::get(d);
  //def<double>(d, names::bias,bias_);//?
}

void nest::STDP_learning_behavior::set(const DictionaryDatum& d)
{
  Supervised_behavior::set(d);
  //updateValue<double>(d, names::bias,bias_);//?
}

/************************************************************************/
/*online learning*/

nest::Supervised_behavior_online::Supervised_behavior_online(const DictionaryDatum &d)
:Behavior_learning(),iterator_desired(0)
{
  updateValue<double>(d, names::h, h);
  //updateValue<bool>(d, names::trig_learn_oper, trig_learn_oper_);
  trig_learn_oper_=false;
  //updateValue<bool>(d, names::trig_seal, trig_seal_);
  trig_seal_=false;
  //updateValue<int>(d, names::number_prim, number_prim_);
  number_prim_=0;
  //updateValue<double>(d, names::bias, bias_);
  bias_=0;
}

void nest::Supervised_behavior_online::read_learn_set_desired_from_file(std::string s_tmp)
{
  std::vector<double> tmp;
  std::istringstream iss;
  std::ifstream f(s_tmp.c_str());
  std::string line;
  if (f.is_open())
  {
    while ( f.good() )
    {
	getline (f,line);
	iss.str(line);
	std::copy (std::istream_iterator<double>(iss), std::istream_iterator<double>(), std::back_inserter(tmp));
	//std::cout<<"read_learn inside"<<"\n";
	v_learn_set_desired.push_back(tmp);
	tmp.clear();
	iss.clear();
    }
  }
  //std::cout<<"read_learn "<<"\n";
  f.clear();
  f.close();
}


void nest::Supervised_behavior_online::read_learn_set_input_from_files(std::string s_directory_learn_set_input,int id_neuron)
{
  number_of_incoming_synapses=66;
}

void nest::Supervised_behavior_online::reset_iterator_desired()
{
  iterator_desired=0;
}

int nest::Supervised_behavior_online::get_number_of_incoming_synapses()
{
  return number_of_incoming_synapses;
}

bool nest::Supervised_behavior_online::is_desired_spike(double time,std::vector<double > &v_learn_desired)
{
  bool is_spike;
  //std::cout<<"iterator_desired "<<iterator_desired<<"\t";
  //std::cout<<"time_desired "<<v_learn_desired[iterator_desired]<<"\t";
  //std::cout<<time<<"\t";
  double time_spike_desired=v_learn_desired[iterator_desired];
  if ((fabs(time_spike_desired-time)<=(h/2.0))&&(time_spike_desired>0))
  {
    is_spike=true;
    //std::cout<<"popalo";
    if (v_learn_desired.size()!=iterator_desired+1)
      iterator_desired++;
  }
  else
  {
    is_spike=false;
  }
  //std::cout<<"exit is desired\n";
  return is_spike;
}

bool nest::Supervised_behavior_online::learn(Behavior_node *behav,DictionaryDatum &d,std::vector<double > &currents)
{
  return false;
}

void nest::Supervised_behavior_online::get(DictionaryDatum &d) const
{
  def<bool>(d, names::trig_learn_oper,trig_learn_oper_);
  def<bool>(d, names::trig_seal,trig_seal_);
  def<long>(d, names::number_prim,number_prim_);
  def<double>(d, names::bias,bias_);
  def<double>(d, names::h,h);
}

void nest::Supervised_behavior_online::set(const DictionaryDatum& d)
{
  updateValue<bool>(d, names::trig_learn_oper,trig_learn_oper_);
  updateValue<bool>(d, names::trig_seal,trig_seal_);
  updateValue<long>(d, names::number_prim,number_prim_);
  updateValue<double>(d, names::bias,bias_);
  updateValue<double>(d, names::h,h);
}

long nest::Supervised_behavior_online::get_number_prim()
{
  return number_prim_;
}

double nest::Supervised_behavior_online::get_bias()
{
  return bias_;
}