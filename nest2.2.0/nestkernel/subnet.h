/*
 *  subnet.h
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

#ifndef SUBNET_H
#define SUBNET_H

#include <vector>
#include <string>
#include "node.h"
#include "dictdatum.h"
#include "multirange.h"

/* BeginDocumentation

Name: subnet - Root node for subnetworks.

Description:
A network node of type subnet serves as a root node for subnetworks

Parameters:
Parameters that can be accessed via the GetStatus and SetStatus functions:

customdict (dictionarytype) -
   A user-defined dictionary, which may be used to store additional
   data.
label (stringtype) -
  A user-defined string, which may be used to give a symbolic name to
  the node.
number_of_children (integertype) -
  The number of direct children of the subnet

  SeeAlso: modeldict, Node
*/

namespace nest{

  using std::vector;

  class Node;

  /**
   * Base class for all subnet nodes.
   * This class can be used
   * - to group other Nodes into "sub-networks"
   * - to construct Node classes which are composed of multiple 
   *   subnodes.
   */
  class Subnet: public Node
  {
  public:

    Subnet();   

    Subnet(const Subnet &);

    virtual ~Subnet(){}
   
    void set_status(const DictionaryDatum&);
    void get_status(DictionaryDatum&) const;

    bool has_proxies() const;
          
    size_t global_size() const;  //!< Returns total number of children.
    size_t local_size() const; //!< Returns number of childern in local process.
    bool   global_empty() const; //!< returns true if subnet is empty *globally*
    bool   local_empty() const; //!< returns true if subnet has no local nodes

    void   reserve(size_t);

    /**
     * Add a local node to the subnet.
     * This function adds a node to the subnet and returns its local id.
     * The node is appended to the subnet child-list. 
     */ 
    index add_node(Node *);

    /**
     * Add a remote node to the subnet.
     * This function increments the next local id to be assigned.
     */ 
    index add_remote_node(index gid, index mid);

    /**
     * Return iterator to the first local child node.
     */
    vector<Node*>::iterator local_begin();

    /**
     * Return iterator to the end of the local child-list.
     */
    vector<Node*>::iterator local_end();

    /**
     * Return const iterator to the first local child node.
     */
    vector<Node*>::const_iterator local_begin() const;

    /**
     * Return const iterator to the end of the local child-list.
     */
    vector<Node*>::const_iterator local_end() const;

    /**
     * Return pointer to Node at given LID if it is local.
     * @note Defined for dense subnets only (all children local)
     */
    Node* at_lid(index) const;

    /**
     * Return the subnets's user label.
     * Each subnet can be given a user-defined string as a label, which
     * may be used to give a symbolic name to the node. From the SLI
     * level, the GetGlobalNodes command may be used to find a subnet's
     * GID from its label.
     */
    std::string get_label() const;
    
    /**
     * Set the subnet's user label.
     * Each subnet can be given a user-defined string as a label, which
     * may be used to give a symbolic name to the node. From the SLI
     * level, the GetGlobalNodes command may be used to find a subnet's
     * GID from its label.
     */
    void set_label(std::string const);

    /**
     * Set the subnet's custom dictionary.
     * Each subnet can be given a user-defined dictionary, which
     * may be used to store additional data. From the SLI
     * level, the SetStatus command may be used to set a subnet's
     * custom dictionary.
     */
    DictionaryDatum get_customdict() const;
    
    /**
     * Return pointer to the subnet's custom dictionary.
     * Each subnet contains a user-definable dictionary, which
     * may be used to store additional data. From the SLI
     * level, the SetStatus command may be used to set a subnet's
     * custom dictionary.
     */
    void set_customdict(DictionaryDatum const dict);
    
    std::string print_network(int , int, std::string = "");

    bool get_children_on_same_vp() const;
    void set_children_on_same_vp(bool);

    thread get_children_vp() const;
    void set_children_vp(thread);
    
    bool allow_entry() const;

    bool is_homogeneous() const; 

  protected:
    void init_node_(const Node&) {}
    void init_state_(const Node&) {}
    void init_buffers_() {}

    void calibrate() {}
    void update(Time const &, const long_t, const long_t) {}
    
    /**
     * Pointer to child nodes.
     * This vector contains the pointers to the child nodes.
     * Since deletion of Nodes is possible, entries in this
     * vector may be NULL. Note that all code must handle
     * this case gracefully.
     */
    vector<Node *> nodes_;

    /**
     * GIDs of global child nodes.
     * This Multirange contains the GIDs of all child nodes on all
     * processes.
     */
    Multirange gids_;

    /**
     * flag indicating if all children of this subnet have to
     * be created on the same virtual process or not. Use with
     * care. This may lead to severe performance problems!
     */
    bool children_on_same_vp_;
    thread children_vp_;
    
  private:
    void get_dimensions_(std::vector<int>&) const;

    std::string     label_;      //!< user-defined label for this node.
    DictionaryDatum customdict_; //!< user-defined dictionary for this node.
    // note that DictionaryDatum is a pointer and must be initialized in the constructor.
    bool homogeneous_;           //!< flag which indicates if the subnet contains different kinds of models.
    index last_mid_;             //!< model index of last child
  };

  /**
   * Add a local node to the subnet.
   */
  inline
  index Subnet::add_node(Node *n)
  {
    const index lid = gids_.size();
    const index mid = n->get_model_id();
    if ((homogeneous_) && (lid > 0))
      if (mid != last_mid_)
	homogeneous_ = false;
    n->set_lid_(lid);
    n->set_subnet_index_(nodes_.size());
    nodes_.push_back(n);
    n->set_parent_(this);
    gids_.push_back(n->get_gid());
    last_mid_ = mid;
    return lid;
  }
  /**
   * Add a remote node to the subnet.
   */
  inline
  index Subnet::add_remote_node(index gid, index mid)
  {
    const index lid = gids_.size();
    if((homogeneous_) && (lid > 0))
      if (mid != last_mid_)
	homogeneous_ = false;
    last_mid_ = mid;
    gids_.push_back(gid);
    return lid;
  }
  
  inline
  vector<Node*>::iterator Subnet::local_begin()
  {
    return nodes_.begin();
  }

  inline
  vector<Node*>::iterator Subnet::local_end()
  {
    return nodes_.end();
  }

  inline
  vector<Node*>::const_iterator Subnet::local_begin() const
  {
    return nodes_.begin();
  }

  inline
  vector<Node*>::const_iterator Subnet::local_end() const
  {
    return nodes_.end();
  }

  inline
  bool Subnet::local_empty() const
  {
    return nodes_.empty();
  }

  inline
  bool Subnet::global_empty() const
  {
    return gids_.empty();
  }

  inline
  size_t Subnet::global_size() const
  {
    return gids_.size();
  }

  inline
  size_t Subnet::local_size() const
  {
    return nodes_.size();
  }

  inline
  Node* Subnet::at_lid(index lid) const
  {
    // defined for "dense" subnets only
    assert(local_size() == global_size());

    if ( lid >= nodes_.size() )
      throw UnknownNode();

    return nodes_[lid];
  }

  inline 
  void Subnet::reserve(size_t n)
  {
    nodes_.reserve(n);
  }

  inline 
  std::string Subnet::get_label() const
  {
    return label_;
  }

  inline 
  DictionaryDatum Subnet::get_customdict() const
  {
    return customdict_;
  }
  
  inline 
  void Subnet::set_customdict(DictionaryDatum const d)
  {
    customdict_=d;
  }
  
  inline
  bool Subnet::has_proxies() const
  {
    return false;
  }

  inline
  bool Subnet::get_children_on_same_vp() const
  {
    return children_on_same_vp_;
  }

  inline
  void Subnet::set_children_on_same_vp(bool children_on_same_vp)
  {
    children_on_same_vp_ = children_on_same_vp;
  }

  inline
  thread Subnet::get_children_vp() const
  {
    return children_vp_;
  }

  inline
  void Subnet::set_children_vp(thread children_vp)
  {
    children_vp_ = children_vp;
  }

  inline
  bool Subnet::is_homogeneous() const
  {
    return homogeneous_;
  }
  
} // namespace

#endif
