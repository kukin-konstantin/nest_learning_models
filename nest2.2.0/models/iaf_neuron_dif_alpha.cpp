/*
 *  iaf_neuron_dif_alpha.cpp
 *
 *
 */

#include "exceptions.h"
#include "iaf_neuron_dif_alpha.h"
#include "network.h"
#include "dict.h"
#include "integerdatum.h"
#include "doubledatum.h"
#include "dictutils.h"
#include "numerics.h"
#include "universal_data_logger_impl.h"

#include <limits>


/* ---------------------------------------------------------------- 
 * Recordables map
 * ---------------------------------------------------------------- */



int nest::iaf_neuron_dif_alpha::id_gl=-1; //nebodimo dlya togochtobi numeracia neuronov nacinalas s edinici
nest::RecordablesMap<nest::iaf_neuron_dif_alpha> nest::iaf_neuron_dif_alpha::recordablesMap_;


namespace nest
{
  // Override the create() method with one call to RecordablesMap::insert_() 
  // for each quantity to be recorded.
  template <>
  void RecordablesMap<iaf_neuron_dif_alpha>::create()
  {
    // use standard names whereever you can for consistency!
    insert_(names::V_m, &iaf_neuron_dif_alpha::get_V_m_);
  }
}

/* ---------------------------------------------------------------- 
 * Default constructors defining default parameters and state
 * ---------------------------------------------------------------- */
    
nest::iaf_neuron_dif_alpha::Parameters_::Parameters_()
  : C_      (200.0    ),  // pF
    Tau_    ( 10.0    ),  // ms
    TauR_   (  2.0    ),  // ms
    E_L_     (-70.0    ),  // mV
    V_th_  (-54.0),  // mV
    I_e_    (  0.0    )   // pA
    //n_in_(0)
    
{}

nest::iaf_neuron_dif_alpha::State_::State_()
  : 
  v_(-70.0),
    I_(0.0),  
    r_ (0)
{}

nest::iaf_neuron_dif_alpha::~iaf_neuron_dif_alpha()
{
  /*std::cout<<"destruct learn_set_input "<<"\n";
   for (int i=0;i!=behav_learning_loc->v_learn_set_input.size();i++)
   {
     for (int j=0;j!=behav_learning_loc->v_learn_set_input[i].size();j++)
     {
       for (int k=0;k!=behav_learning_loc->v_learn_set_input[i][j].size();k++)
       {
	  std::cout<<behav_learning_loc->v_learn_set_input[i][j][k]<<"\t";
       }
        std::cout<<"\n";
     }
     std::cout<<"\n\n";
   }*/
  delete behav_loc;
  delete behav_learning_loc;
  std::cout<<"iaf_neuron_dif_alpha destruct counter="<<iaf_neuron_dif_alpha::id_gl<<"\n";
  iaf_neuron_dif_alpha::id_gl--;
  //std::cout<<"behav"<<behav<<"\n";
  //delete behav;
   //std::cout<<"iaf destruct "<<"\n";
   /*std::map<std::pair<index,index>,Synapse_variables >::iterator it=m_incoming_synapses.begin();
   double sum=0;
	while (it!=m_incoming_synapses.end())
	{
	    double sum_in=0;
	    double t_spike;
	    std::cout<<"index1 "<<(*it).first.first<<"\t";
	    std::cout<<"index2 "<<(*it).first.second<<"\t";
	    std::cout<<"g_max "<<(*it).second.g_max<<"\t";
	    std::cout<<"E_rev "<<(*it).second.E_rev<<"\t";
	    std::cout<<"weight "<<(*it).second.weight<<"\n";
	    for (int i=0;i!=(*it).second.v_time_spike_times.size();i++)
	    {
		t_spike=(*it).second.v_time_spike_times[i];
		std::cout<<"t_spike "<<t_spike<<"\n";
	    }
	    
	    it++;
	}*/
   /*std::cout<<"destruct learn_set_desired "<<"\n";
   for (int i=0;i!=v_learn_set_desired.size();i++)
   {
     for (int j=0;j!=v_learn_set_desired[i].size();j++)
     {
	std::cout<<v_learn_set_desired[i][j]<<"\t";
     }
     std::cout<<"\n";
   }*/
}



/* ---------------------------------------------------------------- 
 * Parameter and state extractions and manipulation functions
 * ---------------------------------------------------------------- */

void nest::iaf_neuron_dif_alpha::Parameters_::get(DictionaryDatum &d) const
{
  
  def<double>(d, names::E_L, E_L_);   // Resting potential
  def<double>(d, names::I_e, I_e_);
  def<double>(d, names::V_th, V_th_); // threshold value
  def<double>(d, names::C_m, C_);
  def<double>(d, names::tau_m, Tau_);
  def<double>(d, names::t_ref, TauR_);
  //dobavit peremennie obuchenia
}

void nest::iaf_neuron_dif_alpha::Parameters_::set(const DictionaryDatum& d)
{
  
  updateValue<double>(d, names::E_L, E_L_);
  updateValue<double>(d, names::V_th, V_th_); 
  updateValue<double>(d, names::I_e, I_e_);
  updateValue<double>(d, names::C_m, C_);
  updateValue<double>(d, names::tau_m, Tau_);
  updateValue<double>(d, names::t_ref, TauR_);
  //dobavit peremennie obuchenia
  if ( C_ <= 0 )
    throw BadProperty("Capacitance must be strictly positive.");
    
  if ( Tau_ <= 0 || TauR_ <= 0 )
    throw BadProperty("All time constants must be strictly positive.");
}

void nest::iaf_neuron_dif_alpha::State_::get(DictionaryDatum &d, const Parameters_& p) const
{
  def<double>(d, names::V_m, v_); // Membrane potential
}

void nest::iaf_neuron_dif_alpha::State_::set(const DictionaryDatum& d, const Parameters_& p)
{
  updateValue<double>(d, names::V_m, v_);
}


nest::iaf_neuron_dif_alpha::Buffers_::Buffers_(iaf_neuron_dif_alpha &n)
  : logger_(n)
{}

nest::iaf_neuron_dif_alpha::Buffers_::Buffers_(const Buffers_ &, iaf_neuron_dif_alpha &n)
  : logger_(n)
{}


/* ---------------------------------------------------------------- 
 * Default and copy constructor for node
 * ---------------------------------------------------------------- */

nest::iaf_neuron_dif_alpha::iaf_neuron_dif_alpha()
  : Archiving_Node(), 
    P_(), 
    S_(),
    B_(*this)
{
  recordablesMap_.create();
  behav=new Iaf_neuron_dif_alpha_behavior;
  //
  DictionaryDatum *d = new DictionaryDatum(new Dictionary);
  //behav_learning=new STDP_learning_behavior((*d));
  behav_learning=new Gradient_descent_behavior((*d));
  delete d;
  std::cout<<"non-copy iaf_neuron_dif_alpha::id_gl "<<iaf_neuron_dif_alpha::id_gl<<"\n"; 
  behav_loc= dynamic_cast<Iaf_neuron_dif_alpha_behavior *>(behav);
  //behav_learning_loc=dynamic_cast<STDP_learning_behavior *>(behav_learning);
  behav_learning_loc=dynamic_cast<Gradient_descent_behavior *>(behav_learning);
  //vnesti suda Behavior_learning
  //nest::iaf_neuron_dif_alpha::sh_gl=0;
}

nest::iaf_neuron_dif_alpha::iaf_neuron_dif_alpha(const iaf_neuron_dif_alpha& n)
  : Archiving_Node(n), 
    P_(n.P_), 
    S_(n.S_),
    B_(n.B_, *this)
{
  behav=new Iaf_neuron_dif_alpha_behavior;
  DictionaryDatum *d = new DictionaryDatum(new Dictionary);
  //behav_learning=new STDP_learning_behavior((*d));
  behav_learning=new Gradient_descent_behavior((*d));
  delete d;
  //std::cout<<"iaf born by copy "<<"behav too"<<"\n"; 
  behav_loc= dynamic_cast<Iaf_neuron_dif_alpha_behavior *>(behav);
  //behav_learning_loc=dynamic_cast<STDP_learning_behavior *>(behav_learning);
  behav_learning_loc=dynamic_cast<Gradient_descent_behavior *>(behav_learning);
  /*learn function*/
  //vnesti suda Behavior_learning
  std::stringstream r;
  iaf_neuron_dif_alpha::id_gl++;
  id_neuron=iaf_neuron_dif_alpha::id_gl;
  std::cout<<"copy iaf_neuron_dif_alpha::id_gl "<<iaf_neuron_dif_alpha::id_gl<<"\n";
  r << iaf_neuron_dif_alpha::id_gl;
  //r<<++nest::iaf_neuron_dif_alpha::sh_gl;
  //id_gl=nest::iaf_neuron_dif_alpha::sh_gl;
  std::string s_directory_learn_set_desired="/home/neuron/animat_work/learn_set/desired/neuron_"+r.str()+".txt";
  //std::string s_directory_learn_set_input="/home/neuron/animat_work/learn_set/input/neuron_";
  //std::cout<<"s_tmp "<<s_tmp<<"\n";
  behav_learning_loc->read_learn_set_desired_from_file(s_directory_learn_set_desired);//perepisat vizov
   //behav_learning_loc->read_learn_set_input_from_files(s_directory_learn_set_input,iaf_neuron_dif_alpha::id_gl); //ano
  /*learn function*/
}

void nest::iaf_neuron_dif_alpha::set_learning_variables(DictionaryDatum &d,double t,std::vector<double > &currents)//prepare before learning
{
	/*std::vector<double > currents2(currents);
	//
	double learn_time=t-behav_learning_loc->get_bias();
	std::cout<<"learn time "<<learn_time<<"\n";
	int number_prim=behav_learning_loc->get_number_prim();
	std::map<std::pair<index,index>,struct_iaf_neuron_dif_alpha >::iterator it=behav_loc->val.begin();
	int number_of_synap=0;
	while (it!=behav_loc->val.end())
	{
	    double sum_inter=0;
	    double t_spike;
	    double t_tau_a=(*it).second.tau_a;
	    std::cout<<"number_prim "<<number_prim<<"\n";
	    for (int i=0;i!=behav_learning_loc->v_learn_set_input[number_of_synap][number_prim].size();i++)
	    {
		t_spike=behav_learning_loc->v_learn_set_input[number_of_synap][number_prim][i];
		if (t_spike>0.0)
		{
		  sum_inter+=alpha_function(learn_time-t_spike,t_tau_a);
		  std::cout<<"1 w learn_time-t_spike="<<learn_time-t_spike<<"\n";
		  //sum_inter+=heaviside(learn_time-t_spike);
		  //sum_inter+=gauss(learn_time-t_spike,P_.tau_a_,7.2);
		}
	    }
	    currents.push_back(sum_inter);// compare with current calculating another ways
	    it++;
	    number_of_synap++;
	}*/
	
	// it's neccessary to do loop for s
	std::map<std::pair<int,std::pair<index,index> >,struct_iaf_neuron_dif_alpha >::iterator it=behav_loc->val.begin();
	//*std::cout<<"size current="<<behav_loc->val.size()<<"\n";//change
	while (it!=behav_loc->val.end())
	{
	    double sum_inter=0;
	    double t_spike;
	    double t_tau_a=(*it).second.tau_a;
	    for (int i=0;i!=(*it).second.base_val.v_time_spike_times.size();i++)
	    {
		t_spike=(*it).second.base_val.v_time_spike_times[i];
		sum_inter+=alpha_function(t-t_spike,t_tau_a);
		//std::cout<<"2 w t-t_spike="<<t-t_spike<<"\n";
	    }
	    currents.push_back(sum_inter);// compare with current calculating another ways
	    it++;
	}
}


/* ---------------------------------------------------------------- 
 * Node initialization functions
 * ---------------------------------------------------------------- */

void nest::iaf_neuron_dif_alpha::init_state_(const Node& proto)
{
  const iaf_neuron_dif_alpha& pr = downcast<iaf_neuron_dif_alpha>(proto);
  S_ = pr.S_;
}



void nest::iaf_neuron_dif_alpha::init_buffers_()
{
  B_.spikes_.clear();    // includes resize
  B_.currents_.clear();  // include resize
  B_.logger_.reset(); // includes resize
  Archiving_Node::clear_history();
}



void nest::iaf_neuron_dif_alpha::calibrate()
{
  B_.logger_.init();
  V_.rng_ = net_->get_rng(get_thread());
  V_.RefractoryCounts_ = Time(Time::ms(P_.TauR_)).get_steps();
  assert(V_.RefractoryCounts_ >= 0); 
}

/* ---------------------------------------------------------------- 
 * Update and spike handling functions
 * ---------------------------------------------------------------- */

double nest::iaf_neuron_dif_alpha::normalization(double_t v)
{
  double k=0.1;//coff learning
  double tresh=1.0;
  double u=(v-P_.E_L_)/(P_.V_th_);
  if (u<0)
    u=0;
  double z=(u-tresh)/k;
  return z;
}

bool nest::iaf_neuron_dif_alpha::emission_procedure(double_t v)
{
  double Th=1.0;
  double k=0.1;
  const double_t h = Time::get_resolution().get_ms();
  double rand=V_.rng_->drand();
  double lambd_val=std::exp((v-Th)/k);
  //double z=normalization(v);
  double z=v;
  //double sigm=1.0/(1.0+std::exp(-z));
  //double sigm=z;
  double sigm=1.0-std::exp(-lambd_val*h);
  if (rand<=sigm)
  {
    return true;
  }
  else
  {
    return false;
  }
}


void nest::iaf_neuron_dif_alpha::update(Time const & origin, const long_t from, const long_t to)
{
  //bool ok=true;
  assert(to >= 0 && (delay) from < Scheduler::get_min_delay());
  assert(from < to);
  const double_t h = Time::get_resolution().get_ms();
  
  DictionaryDatum *d = new DictionaryDatum(new Dictionary);
  //behav_learning_loc->get((*d));
  double_t v_old;
  get_status_user((*d));
  //get_status((*d));
  double_t t=origin.get_ms();
  bool trig_seal_;
  bool trig_learn_oper_;
  updateValue<bool>((*d), names::trig_seal,trig_seal_);// d->trig_seal_
  updateValue<bool>((*d), names::trig_learn_oper,trig_learn_oper_);
  if (trig_seal_)
  {
    //*std::cout<<"trig_seal\n";
    //std::cout<<"grad "<<"\n";
    trig_seal_=0;
    def<bool>((*d), names::trig_seal, trig_seal_);
    behav_learning_loc->reset_iterator_desired();
    behav_learning_loc->set((*d));
    clear_spike_history();
    S_.I_=0;
  }
  if (trig_learn_oper_)
  {
    //std::cout<<"learn"<<"\n";
    //std::vector<double > v_learn_exam=v_learn_set_desired[P_.number_prim_];
    //std::cout<<v_learn_exam.size()<<"\n";
    //std::cout<<"l I current="<<(P_.I_e_/P_.C_) + (S_.I_/P_.C_)<<"\n";
    for ( long_t lag = from ; lag < to ; ++lag )
    {
	if ( S_.r_ == 0 ) // neuron not refractory
	{
	  v_old = S_.v_;
	  //S_.v_ +=h*( ((P_.E_L_-v_old)/P_.Tau_) + (P_.I_e_/P_.C_) + (S_.I_/P_.C_) )+B_.spikes_.get_value(lag);
	  S_.v_ +=h*( ((P_.E_L_-v_old)/P_.Tau_) + (P_.I_e_/P_.C_) + (S_.I_/P_.C_) );
	  //std::cout<<"learn S_.I_="<<S_.I_<<"\t"<<"S_.v_="<<S_.v_<<"\n";
	}
	else // neuron is absolute refractory
	  --S_.r_;
  
    
	
	//*std::cout<<"index "<<this->get_gid()<<" to "<<to<<"\t"<<"lag "<<lag<<"\t"<<"h "<<h<<"\t time"<<t+h*lag<<"\t"<<" S_.I_="<<S_.I_<<"\n";
	//behav_learning_loc->get_learning_vector(P_.number_prim_)
	//calculate currents
	//add value to t_parametrs
	def<double>((*d), names::time,t+h*lag);// t+h*lag->d
	//std::cout<<"S_.v_="<<S_.v_<<"\t";
	//std::cout<<"norm="<<normalization(S_.v_)<<"\t";
	def<double>((*d), names::potential,S_.v_);
	//std::cout<<"op"<<"\n";
	std::vector<double > currents;	
	set_learning_variables((*d),t+h*lag,currents);
	//std::cout<<"opop"<<"\n";
	
	if (behav_learning_loc->learn(behav_loc,(*d),currents))
	{
		//*std::cout<<"spike_time"<<"\n";
		//clean spike history
		S_.r_ = V_.RefractoryCounts_;
		S_.v_= P_.E_L_;   
		set_spiketime(Time::step(origin.get_steps()+lag));
		//std::cout<<"learn time spike from neuron "<<Time::step(origin.get_steps()+lag)<<"\n";
		SpikeEvent se;
		//se.set_weight(-1.0);
		network()->send(*this, se, lag);
		//clear_spike_history();
		//ok=false;
	}
	else if (S_.v_ >= P_.V_th_)
	{
		//std::cout<<"spike_time"<<"\n";
		//clean spike history
		S_.r_ = V_.RefractoryCounts_;
		S_.v_= P_.E_L_;   
		set_spiketime(Time::step(origin.get_steps()+lag)); // change -1
		SpikeEvent se;
		//se.set_weight(-1.0);
		network()->send(*this, se, lag);
 		//clear_spike_history();
		//ok=false;
	}
	//B_.spikes_.get_value(lag);// Need for clean buffer. It's necessary for working log-device wrapper.
	S_.I_=B_.currents_.get_value(lag);
	S_.I_ =S_.I_+get_sum_I(t+h*lag); //two time calclulate
	//std::cout<<"S_.I_="<<S_.I_<<"\n";;
	// voltage logging
	B_.logger_.record_data(origin.get_steps()+lag);
      }
      //std::cout<<"opop3"<<"\n";
  }
  else /*function mode */
  {
    //std::cout<<"function"<<"\n";
    //std::cout<<"I current="<<(P_.I_e_/P_.C_) + (S_.I_/P_.C_)<<"\t";
    for ( long_t lag = from ; lag < to ; ++lag )
    {
	if ( S_.r_ == 0 ) // neuron not refractory
	{
	  v_old = S_.v_;
	  //S_.v_ +=h*( ((P_.E_L_-v_old)/P_.Tau_) + (P_.I_e_/P_.C_) + (S_.I_/P_.C_) )+B_.spikes_.get_value(lag);
	  S_.v_ +=h*( ((P_.E_L_-v_old)/P_.Tau_) + (P_.I_e_/P_.C_) + (S_.I_/P_.C_) );
	  //std::cout<<"fun S_.I_="<<S_.I_<<"\t"<<"S_.v_="<<S_.v_<<"\n";
	  std::cout<<"fun S_.I_="<<S_.I_<<"\n";
	}
	else // neuron is absolute refractory
	  --S_.r_;
  
    
	// threshold crossing
	std::cout<<"index "<<this->get_gid()<<" to "<<to<<"\t"<<"lag "<<lag<<"\t"<<"h "<<h<<"\t";
	std::cout<<S_.v_<<"\t"<<"\n";
	//if (emission_procedure(S_.v_))
	if (S_.v_ >= P_.V_th_)
	{
		//std::cout<<"spike_time"<<"\n";
		//clean spike history
		S_.r_ = V_.RefractoryCounts_;
		S_.v_= P_.E_L_;   
		set_spiketime(Time::step(origin.get_steps()+lag)); // change -1
		//std::cout<<"fun time spike from neuron "<<Time::step(origin.get_steps()+lag)<<"\n";// change -1
		SpikeEvent se;
		//se.set_weight(-1.0);
		network()->send(*this, se, lag);
 		//clear_spike_history();
		//ok=false;
	}
     
	//B_.spikes_.get_value(lag);// Need for clean buffer. It's necessary for working log-device wrapper.
	S_.I_=B_.currents_.get_value(lag);
	
	S_.I_ =S_.I_+get_sum_I(t+h*lag); // надо делать инач
	//std::cout<<"fun S_.I_="<<S_.I_<<"\n";
	// voltage logging
	B_.logger_.record_data(origin.get_steps()+lag);
      }
      //std::cout<<"\n";
      /*function mode */
  }
  delete d;
}                           
                     
void nest::iaf_neuron_dif_alpha::handle(SpikeEvent & e)
{
  assert(e.get_delay() > 0);

  B_.spikes_.add_value(e.get_rel_delivery_steps(network()->get_slice_origin()),
		       e.get_weight() * e.get_multiplicity() );
  //std::cout<<"spike handle "<<e.get_stamp().get_ms()<<"\n";
}

void nest::iaf_neuron_dif_alpha::handle(CurrentEvent& e)
{
  assert(e.get_delay() > 0);

  const double_t c=e.get_current();
  const double_t w=e.get_weight();

  B_.currents_.add_value(e.get_rel_delivery_steps(network()->get_slice_origin()), 
			 w * c);
}

void nest::iaf_neuron_dif_alpha::handle(DataLoggingRequest& e)
{
  B_.logger_.handle(e);
}
