
/*
 *  stdp_poisson.cpp
 *
 *  This file is part of NEST.
 *
 *  Copyright (C) 2004 The NEST Initiative
 *
 *  NEST is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  NEST is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with NEST.  If not, see <http://www.gnu.org/licenses/>.
 *
 */
#include "network.h"
#include "dictdatum.h"
#include "connector_model.h"
#include "common_synapse_properties.h"
#include "stdp_poisson.h"
#include "event.h"


namespace nest
{

  STDPPoisson::STDPPoisson() :
    ConnectionHetWD(),
    U_(0.5),
    D_(1100.0),
    F_(50.0),
    u_(0.5),
    R(1.0),
    tau_plus_(20.0),
    lambda_(0.01),
    alpha_(1.0),
    mu_plus_(1.0),
    mu_minus_(1.0),
    Wmax_(100.0),
    Kplus_(0.0)
  { }


  STDPPoisson::STDPPoisson(const STDPPoisson &rhs) :
    ConnectionHetWD(rhs)
  {
    U_= rhs.U_;
    D_= rhs.D_;
    F_= rhs.F_;
    u_= rhs.u_;
    R= rhs.R;
    tau_plus_ = rhs.tau_plus_;
    lambda_ = rhs.lambda_;
    alpha_ = rhs.alpha_;
    mu_plus_ = rhs.mu_plus_;
    mu_minus_ = rhs.mu_minus_;
    Wmax_ = rhs.Wmax_;
    Kplus_ = rhs.Kplus_;
  }

  void STDPPoisson::get_status(DictionaryDatum & d) const
  {
    ConnectionHetWD::get_status(d);
    def<double_t>(d, "U", U_);
    def<double_t>(d, "D", D_);
    def<double_t>(d, "F", F_);
    def<double_t>(d, "u", u_);
    def<double_t>(d, "R", R);
    def<double_t>(d, "tau_plus", tau_plus_);
    def<double_t>(d, "lambda", lambda_);
    def<double_t>(d, "alpha", alpha_);
    def<double_t>(d, "mu_plus", mu_plus_);
    def<double_t>(d, "mu_minus", mu_minus_);
    def<double_t>(d, "Wmax", Wmax_);
  }

  void STDPPoisson::set_status(const DictionaryDatum & d, ConnectorModel &cm)
  {
    ConnectionHetWD::set_status(d, cm);
    updateValue<double_t>(d, "U", U_);
    updateValue<double_t>(d, "D", D_);
    updateValue<double_t>(d, "F", F_);
    updateValue<double_t>(d, "u", u_);
    updateValue<double_t>(d, "R", R);
    updateValue<double_t>(d, "tau_plus", tau_plus_);
    updateValue<double_t>(d, "lambda", lambda_);
    updateValue<double_t>(d, "alpha", alpha_);
    updateValue<double_t>(d, "mu_plus", mu_plus_);
    updateValue<double_t>(d, "mu_minus", mu_minus_);
    updateValue<double_t>(d, "Wmax", Wmax_);
  }

   /**
   * Set properties of this connection from position p in the properties
   * array given in dictionary.
   */
  void STDPPoisson::set_status(const DictionaryDatum & d, index p, ConnectorModel &cm)
  {
    ConnectionHetWD::set_status(d, p, cm);
    set_property<double_t>(d, "Us", p, U_);
    set_property<double_t>(d, "Ds", p, D_);
    set_property<double_t>(d, "Fs", p, F_);
    set_property<double_t>(d, "us", p, u_);
    set_property<double_t>(d, "Rs", p, R);
    set_property<double_t>(d, "tau_pluss", p, tau_plus_);
    set_property<double_t>(d, "lambdas", p, lambda_);
    set_property<double_t>(d, "alphas", p, alpha_);
    set_property<double_t>(d, "mu_pluss", p, mu_plus_);
    set_property<double_t>(d, "mu_minuss", p, mu_minus_);
    set_property<double_t>(d, "Wmaxs", p, Wmax_);
  }

  void STDPPoisson::initialize_property_arrays(DictionaryDatum & d) const
  {
    ConnectionHetWD::initialize_property_arrays(d);
    initialize_property_array(d, "Us");
    initialize_property_array(d, "Ds");
    initialize_property_array(d, "Fs");
    initialize_property_array(d, "us");
    initialize_property_array(d, "Rs");
    initialize_property_array(d, "tau_pluss"); 
    initialize_property_array(d, "lambdas"); 
    initialize_property_array(d, "alphas"); 
    initialize_property_array(d, "mu_pluss"); 
    initialize_property_array(d, "mu_minuss");
    initialize_property_array(d, "Wmaxs");
  }

  /**
   * Append properties of this connection to the given dictionary. If the
   * dictionary is empty, new arrays are created first.
   */
  void STDPPoisson::append_properties(DictionaryDatum & d) const
  {
    ConnectionHetWD::append_properties(d);
    append_property<double_t>(d, "Us", U_);
    append_property<double_t>(d, "Ds", D_);
    append_property<double_t>(d, "Fs", F_);
    append_property<double_t>(d, "us", u_);
    append_property<double_t>(d, "Rs", R);
    append_property<double_t>(d, "tau_pluss", tau_plus_); 
    append_property<double_t>(d, "lambdas", lambda_); 
    append_property<double_t>(d, "alphas", alpha_); 
    append_property<double_t>(d, "mu_pluss", mu_plus_); 
    append_property<double_t>(d, "mu_minuss", mu_minus_);
    append_property<double_t>(d, "Wmaxs", Wmax_);
  }

} // of namespace nest
