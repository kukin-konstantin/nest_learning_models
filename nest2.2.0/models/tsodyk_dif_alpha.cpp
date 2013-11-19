
/*
 *  tsodyk_dif_alpha.cpp
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
#include "tsodyk_dif_alpha.h"
#include "event.h"


namespace nest
{

  TsodykDifAlpha::TsodykDifAlpha() :
    ConnectionHetWD(),
    U_(0.5),
    D_(1100.0),
    F_(50.0),
    u_(0.5),
    R(1.0),
    tau_a_(10.0),
    t_count(0)
  { }


  TsodykDifAlpha::TsodykDifAlpha(const TsodykDifAlpha &rhs) :
    ConnectionHetWD(rhs)
  {
    U_= rhs.U_;
    D_= rhs.D_;
    F_= rhs.F_;
    u_= rhs.u_;
    R= rhs.R;
    tau_a_=rhs.tau_a_;
    t_count=rhs.t_count;
  }

  void TsodykDifAlpha::get_status(DictionaryDatum & d) const
  {
    ConnectionHetWD::get_status(d);
    def<double_t>(d, "U", U_);
    def<double_t>(d, "D", D_);
    def<double_t>(d, "F", F_);
    def<double_t>(d, "u", u_);
    def<double_t>(d, "R", R);
    def<double_t>(d, names::tau_a, tau_a_);
  }

  void TsodykDifAlpha::set_status(const DictionaryDatum & d, ConnectorModel &cm)
  {
    ConnectionHetWD::set_status(d, cm);
    updateValue<double_t>(d, "U", U_);
    updateValue<double_t>(d, "D", D_);
    updateValue<double_t>(d, "F", F_);
    updateValue<double_t>(d, "u", u_);
    updateValue<double_t>(d, "R", R);
    updateValue<double_t>(d, names::tau_a, tau_a_);
  }

   /**
   * Set properties of this connection from position p in the properties
   * array given in dictionary.
   */
  void TsodykDifAlpha::set_status(const DictionaryDatum & d, index p, ConnectorModel &cm)
  {
    ConnectionHetWD::set_status(d, p, cm);
    set_property<double_t>(d, "Us", p, U_);
    set_property<double_t>(d, "Ds", p, D_);
    set_property<double_t>(d, "Fs", p, F_);
    set_property<double_t>(d, "us", p, u_);
    set_property<double_t>(d, "Rs", p, R);
    set_property<double_t>(d, "tau_as", p, tau_a_);
  }

  void TsodykDifAlpha::initialize_property_arrays(DictionaryDatum & d) const
  {
    ConnectionHetWD::initialize_property_arrays(d);
    initialize_property_array(d, "Us");
    initialize_property_array(d, "Ds");
    initialize_property_array(d, "Fs");
    initialize_property_array(d, "us");
    initialize_property_array(d, "Rs");
    initialize_property_array(d, "tau_as");
  }

  /**
   * Append properties of this connection to the given dictionary. If the
   * dictionary is empty, new arrays are created first.
   */
  void TsodykDifAlpha::append_properties(DictionaryDatum & d) const
  {
    ConnectionHetWD::append_properties(d);
    append_property<double_t>(d, "Us", U_);
    append_property<double_t>(d, "Ds", D_);
    append_property<double_t>(d, "Fs", F_);
    append_property<double_t>(d, "us", u_);
    append_property<double_t>(d, "Rs", R);
    append_property<double_t>(d, "tau_as", tau_a_);
  }

} // of namespace nest
