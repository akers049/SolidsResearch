#ifndef dealii__mapping_q_eulerian_hp_h
#define dealii__mapping_q_eulerian_hp_h

#include <deal.II/base/config.h>
#include <deal.II/base/smartpointer.h>
#include <deal.II/base/thread_management.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/mapping_q.h>
#include <deal.II/hp/dof_handler.h>
#include <deal.II/hp/fe_values.h>
#include <deal.II/fe/fe_series.h>

#include <deal.II/base/utilities.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/lac/vector.h>
#include <deal.II/lac/block_vector.h>
#include <deal.II/lac/la_vector.h>
#include <deal.II/lac/parallel_vector.h>
#include <deal.II/lac/parallel_block_vector.h>
#include <deal.II/lac/trilinos_vector.h>
#include <deal.II/lac/trilinos_block_vector.h>
#include <deal.II/lac/trilinos_parallel_block_vector.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/fe/fe.h>
#include <deal.II/fe/fe_tools.h>
#include <deal.II/hp/dof_handler.h>
#include <deal.II/hp/fe_values.h>
#include <deal.II/fe/fe_series.h>


DEAL_II_NAMESPACE_OPEN


template <typename> class Vector;


/*!@addtogroup mapping */
/*@{*/

/**
 * This class is an extension of the MappingQ1Eulerian class to higher order
 * $Q_p$ mappings.  It is useful when one wants to calculate shape function
 * information on a domain that is deforming as the computation proceeds.
 *
 * <h3>Usage</h3>
 *
 * The constructor of this class takes three arguments: the polynomial degree
 * of the desire Qp mapping, a reference to the vector that defines the
 * mapping from the initial configuration to the current configuration, and a
 * reference to the DoFHandler. The most common case is to use the solution
 * vector for the problem under consideration as the shift vector. The key
 * requirement is that the number of components of the given vector field be
 * equal to (or possibly greater than) the number of space dimensions. If
 * there are more components than space dimensions (for example, if one is
 * working with a coupled problem where there are additional solution
 * variables), the first <tt>dim</tt> components are assumed to represent the
 * displacement field, and the remaining components are ignored.  If this
 * assumption does not hold one may need to set up a separate DoFHandler on
 * the triangulation and associate the desired shift vector to it.
 *
 * Typically, the DoFHandler operates on a finite element that is constructed
 * as a system element (FESystem) from continuous FE_Q() objects. An example
 * is shown below:
 * @code
 *    FESystem<dim> fe(FE_Q<dim>(2), dim, FE_Q<dim>(1), 1);
 *    DoFHandler<dim> dof_handler(triangulation);
 *    dof_handler.distribute_dofs(fe);
 *    Vector<double> displacement_field(dof_handler.n_dofs());
 *    // ... compute displacement field somehow...
 *    MappingQEulerian<dim> q2_mapping(2, dof_handler, displacement_field);
 * @endcode
 *
 * In this example, our element consists of <tt>(dim+1)</tt> components. Only
 * the first <tt>dim</tt> components will be used, however, to define the Q2
 * mapping.  The remaining components are ignored.
 *
 * Note that it is essential to call the distribute_dofs(...) function before
 * constructing a mapping object.
 *
 * Also note that since the vector of shift values and the dof handler are
 * only associated to this object at construction time, you have to make sure
 * that whenever you use this object, the given objects still represent valid
 * data.
 *
 * To enable the use of the MappingQ1Eulerian class also in the context of
 * parallel codes using the PETSc or Trilinos wrapper classes, the type
 * of the vector can be specified as template parameter <tt>VectorType</tt>.
 *
 * @author Joshua White, 2008
 */
template <int dim, typename VectorType = Vector<double>, int spacedim=dim >
class MappingQEulerian_hp : public MappingQ<dim, spacedim>
{
public:
  /**
   * Constructor.
   *
   * @param[in] degree The polynomial degree of the desired $Q_p$ mapping.
   * @param[in] euler_dof_handler A DoFHandler object that defines a finite
   * element space. This space needs to have at least dim components and the
   * first dim components of the space will be considered displacements
   * relative to the original positions of the cells of the triangulation.
   * @param[in] euler_vector A finite element function in the space defined by
   * the second argument. The first dim components of this function will be
   * interpreted as the displacement we use in defining the mapping, relative
   * to the location of cells of the underlying triangulation.
   */
  MappingQEulerian_hp (const unsigned int        degree,
                    const hp::DoFHandler<dim,spacedim> &euler_dof_handler,
                    const VectorType               &euler_vector);
  virtual ~MappingQEulerian_hp(){};

  /**
   * @deprecated Use the constructor with the reverse order of second and
   * third argument.
   */
  // MappingQEulerian_hp (const unsigned int        degree,
  //                   const VectorType               &euler_vector,
  //                   const hp::DoFHandler<dim,spacedim> &euler_dof_handler) DEAL_II_DEPRECATED;

  /**
   * Return the mapped vertices of the cell. For the current class, this
   * function does not use the support points from the geometry of the current
   * cell but instead evaluates an externally given displacement field in
   * addition to the geometry of the cell.
   */
  virtual
  std_cxx11::array<Point<spacedim>, GeometryInfo<dim>::vertices_per_cell>
  get_vertices (const typename Triangulation<dim,spacedim>::cell_iterator &cell) const;

  /**
   * Return a pointer to a copy of the present object. The caller of this copy
   * then assumes ownership of it.
   */
  virtual
  Mapping<dim,spacedim> *clone () const;

  /**
   * Always returns @p false because MappingQ1Eulerian does not in general
   * preserve vertex locations (unless the translation vector happens to
   * provide for zero displacements at vertex locations).
   */
  bool preserves_vertex_locations () const;

  /**
   * Exception
   */
  DeclException0 (ExcInactiveCell);

protected:
  /**
   * Compute mapping-related information for a cell. See the documentation of
   * Mapping::fill_fe_values() for a discussion of purpose, arguments, and
   * return value of this function.
   *
   * This function overrides the function in the base class since we cannot
   * use any cell similarity for this class.
   */
  virtual
  CellSimilarity::Similarity
  fill_fe_values (const typename Triangulation<dim,spacedim>::cell_iterator &cell,
                  const CellSimilarity::Similarity                           cell_similarity,
                  const Quadrature<dim>                                     &quadrature,
                  const typename Mapping<dim,spacedim>::InternalDataBase    &internal_data,
                  internal::FEValues::MappingRelatedData<dim,spacedim>      &output_data) const;

  /**
   * Reference to the vector of shifts.
   */
  SmartPointer<const VectorType, MappingQEulerian_hp<dim,VectorType,spacedim> > euler_vector;

  /**
   * Pointer to the DoFHandler to which the mapping vector is associated.
   */
  SmartPointer<const hp::DoFHandler<dim,spacedim>,MappingQEulerian_hp<dim,VectorType,spacedim> > euler_dof_handler;


private:

  /**
   * A class derived from MappingQGeneric that provides the generic mapping
   * with support points on boundary objects so that the corresponding Q3
   * mapping ends up being C1.
   */
  class MappingQEulerianGeneric_hp : public MappingQGeneric<dim,spacedim>
  {
  public:

    /**
     * Constructor.
     */
    MappingQEulerianGeneric_hp (const unsigned int        degree,
                             const MappingQEulerian_hp<dim,VectorType,spacedim> &mapping_q_eulerian);
    ~MappingQEulerianGeneric_hp(){delete fe_values;};

    /**
     * Return the mapped vertices of the cell. For the current class, this
     * function does not use the support points from the geometry of the
     * current cell but instead evaluates an externally given displacement
     * field in addition to the geometry of the cell.
     */
    virtual
    std_cxx11::array<Point<spacedim>, GeometryInfo<dim>::vertices_per_cell>
    get_vertices (const typename Triangulation<dim,spacedim>::cell_iterator &cell) const;

    /**
     * Compute the positions of the support points in the current
     * configuration. See the documentation of
     * MappingQGeneric::compute_mapping_support_points() for more information.
     */
    virtual
    std::vector<Point<spacedim> >
    compute_mapping_support_points(const typename Triangulation<dim,spacedim>::cell_iterator &cell) const;

  private:
    /**
     * Reference to the surrounding object off of which we live.
     */
    const MappingQEulerian_hp<dim,VectorType,spacedim> &mapping_q_eulerian;


    /**
     * Special quadrature rule used to define the support points in the
     * reference configuration.
     */

    class SupportQuadrature : public Quadrature<dim>
    {
    public:
      /**
       * Constructor, with an argument defining the desired polynomial degree.
       */

      SupportQuadrature (const unsigned int map_degree);

    };

    /**
     * A member variable holding the quadrature points in the right order.
     */
    hp::QCollection<dim> support_quadrature_collection;

    /**
     * FEValues object used to query the the given finite element field at the
     * support points in the reference configuration.
     *
     * The variable is marked as mutable since we have to call
     * FEValues::reinit from compute_mapping_support_points, a function that
     * is 'const'.
     */
    mutable hp::FEValues<dim,spacedim>* fe_values = NULL;

    /**
     * A variable to guard access to the fe_values variable.
     */
    mutable Threads::Mutex fe_values_mutex;
  };

};


template <int dim, typename VectorType, int spacedim>
inline
bool
MappingQEulerian_hp<dim,VectorType,spacedim>::preserves_vertex_locations () const
{
  return false;
}


// ---------------------------------------------------------------------
//
// Copyright (C) 2001 - 2016 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE at
// the top level of the deal.II distribution.
//
// ---------------------------------------------------------------------



// .... MAPPING Q EULERIAN CONSTRUCTOR

template <int dim, class VectorType, int spacedim>
MappingQEulerian_hp<dim, VectorType, spacedim>::MappingQEulerianGeneric_hp::
MappingQEulerianGeneric_hp (const unsigned int                                 degree,
                         const MappingQEulerian_hp<dim, VectorType, spacedim> &mapping_q_eulerian)
  :
  MappingQGeneric<dim,spacedim>(degree),
  mapping_q_eulerian (mapping_q_eulerian)
{
  for(unsigned int i = 0; i < mapping_q_eulerian.euler_dof_handler->get_fe().size(); i++)
    support_quadrature_collection.push_back(SupportQuadrature(degree));

  fe_values = new hp::FEValues<dim, spacedim> (mapping_q_eulerian.euler_dof_handler->get_fe(),
                                          support_quadrature_collection,
                                          update_values | update_quadrature_points);
}

// template <int dim, class VectorType, int spacedim>
// MappingQEulerian_hp<dim, VectorType, spacedim>::
// MappingQEulerian_hp (const unsigned int                degree,
//                   const VectorType               &euler_vector,
//                   const hp::DoFHandler<dim,spacedim> &euler_dof_handler)
//   :
//   MappingQ<dim,spacedim>(degree, true),
//   euler_vector(&euler_vector),
//   euler_dof_handler(&euler_dof_handler)
// {
//   // reset the q1 mapping we use for interior cells (and previously
//   // set by the MappingQ constructor) to a MappingQ1Eulerian with the
//   // current vector
//   this->q1_mapping.reset (new MappingQEulerianGeneric_hp(1, *this));

//   // also reset the qp mapping pointer with our own class
//   this->qp_mapping.reset (new MappingQEulerianGeneric_hp(degree,*this));
// }



template <int dim, class VectorType, int spacedim>
MappingQEulerian_hp<dim, VectorType, spacedim>::
MappingQEulerian_hp (const unsigned int              degree,
                  const hp::DoFHandler<dim,spacedim> &euler_dof_handler,
                  const VectorType               &euler_vector)
  :
  MappingQ<dim,spacedim>(degree, true),
  euler_vector(&euler_vector),
  euler_dof_handler(&euler_dof_handler)
{
  // reset the q1 mapping we use for interior cells (and previously
  // set by the MappingQ constructor) to a MappingQ1Eulerian with the
  // current vector
  this->q1_mapping.reset (new MappingQEulerianGeneric_hp(1, *this));

  // also reset the qp mapping pointer with our own class
  this->qp_mapping.reset (new MappingQEulerianGeneric_hp(degree,*this));
}



template <int dim, class VectorType, int spacedim>
Mapping<dim,spacedim> *
MappingQEulerian_hp<dim, VectorType, spacedim>::clone () const
{
  return new MappingQEulerian_hp<dim,VectorType,spacedim>(this->get_degree(),
                                                       *euler_dof_handler,
                                                       *euler_vector);
}



// .... SUPPORT QUADRATURE CONSTRUCTOR

template <int dim, class VectorType, int spacedim>
MappingQEulerian_hp<dim,VectorType,spacedim>::MappingQEulerianGeneric_hp::
SupportQuadrature::
SupportQuadrature (const unsigned int map_degree)
  :
  Quadrature<dim>(Utilities::fixed_power<dim>(map_degree+1))
{
  // first we determine the support points on the unit cell in lexicographic
  // order, which are (in accordance with MappingQ) the support points of
  // QGaussLobatto.
  const QGaussLobatto<dim> q_iterated(map_degree+1);
  const unsigned int n_q_points = q_iterated.size();

  // we then need to define a renumbering vector that allows us to go from a
  // lexicographic numbering scheme to a hierarchic one.  this fragment is
  // taking almost verbatim from the MappingQ class.
  std::vector<unsigned int> renumber(n_q_points);
  std::vector<unsigned int> dpo(dim+1, 1U);
  for (unsigned int i=1; i<dpo.size(); ++i)
    dpo[i]=dpo[i-1]*(map_degree-1);

  FETools::lexicographic_to_hierarchic_numbering (FiniteElementData<dim> (dpo, 1, map_degree),
                                                  renumber);

  // finally we assign the quadrature points in the required order.
  for (unsigned int q=0; q<n_q_points; ++q)
    this->quadrature_points[renumber[q]] = q_iterated.point(q);
}



// .... COMPUTE MAPPING SUPPORT POINTS

template <int dim, class VectorType, int spacedim>
std_cxx11::array<Point<spacedim>, GeometryInfo<dim>::vertices_per_cell>
MappingQEulerian_hp<dim, VectorType, spacedim>::
get_vertices
(const typename Triangulation<dim,spacedim>::cell_iterator &cell) const
{
  // get the vertices as the first 2^dim mapping support points
  const std::vector<Point<spacedim> > a
    = dynamic_cast<const MappingQEulerianGeneric_hp &>(*this->qp_mapping).compute_mapping_support_points(cell);

  std_cxx11::array<Point<spacedim>, GeometryInfo<dim>::vertices_per_cell> vertex_locations;
  std::copy (a.begin(),
             a.begin()+GeometryInfo<dim>::vertices_per_cell,
             vertex_locations.begin());

  return vertex_locations;
}



template <int dim, class VectorType, int spacedim>
std_cxx11::array<Point<spacedim>, GeometryInfo<dim>::vertices_per_cell>
MappingQEulerian_hp<dim, VectorType, spacedim>::MappingQEulerianGeneric_hp::
get_vertices
(const typename Triangulation<dim,spacedim>::cell_iterator &cell) const
{
  return mapping_q_eulerian.get_vertices (cell);
}




template <int dim, class VectorType, int spacedim>
std::vector<Point<spacedim> >
MappingQEulerian_hp<dim, VectorType, spacedim>::MappingQEulerianGeneric_hp::
compute_mapping_support_points (const typename Triangulation<dim,spacedim>::cell_iterator &cell) const
{
  // first, basic assertion with respect to vector size,

  const types::global_dof_index n_dofs  = mapping_q_eulerian.euler_dof_handler->n_dofs();
  const types::global_dof_index vector_size = mapping_q_eulerian.euler_vector->size();
  (void)n_dofs;
  (void)vector_size;

  AssertDimension(vector_size,n_dofs);

  // we then transform our tria iterator into a dof iterator so we can access
  // data not associated with triangulations
  typename hp::DoFHandler<dim,spacedim>::cell_iterator dof_cell(*cell,
                                                            mapping_q_eulerian.euler_dof_handler);

  Assert (dof_cell->active() == true, ExcInactiveCell());

  // our quadrature rule is chosen so that each quadrature point corresponds
  // to a support point in the undeformed configuration.  we can then query
  // the given displacement field at these points to determine the shift
  // vector that maps the support points to the deformed configuration.

  // we assume that the given field contains dim displacement components, but
  // that there may be other solution components as well (e.g. pressures).
  // this class therefore assumes that the first dim components represent the
  // actual shift vector we need, and simply ignores any components after
  // that.  this implies that the user should order components appropriately,
  // or create a separate dof handler for the displacements.

  // fill shift vector for each support point using an fe_values object. make
  // sure that the fe_values variable isn't used simultaneously from different
  // threads
  Threads::Mutex::ScopedLock lock(fe_values_mutex);
  fe_values->reinit(dof_cell);

  const FEValues<dim> &next_fe_values = fe_values->get_present_fe_values();
  const unsigned int n_support_pts = next_fe_values.n_quadrature_points;
  const unsigned int n_components  = mapping_q_eulerian.euler_dof_handler->get_fe().n_components();

  Assert (n_components >= spacedim, ExcDimensionMismatch(n_components, spacedim) );

  std::vector<Vector<typename VectorType::value_type> >
  shift_vector(n_support_pts,
               Vector<typename VectorType::value_type>(n_components));


  next_fe_values.get_function_values(*mapping_q_eulerian.euler_vector, shift_vector);

  // and finally compute the positions of the support points in the deformed
  // configuration.
  std::vector<Point<spacedim> > a(n_support_pts);
  for (unsigned int q=0; q<n_support_pts; ++q)
    {
      a[q] = next_fe_values.quadrature_point(q);
      for (unsigned int d=0; d<spacedim; ++d)
        a[q](d) += shift_vector[q](d);
    }

  return a;
}



template<int dim, class VectorType, int spacedim>
CellSimilarity::Similarity
MappingQEulerian_hp<dim,VectorType,spacedim>::
fill_fe_values (const typename Triangulation<dim,spacedim>::cell_iterator &cell,
                const CellSimilarity::Similarity                           ,
                const Quadrature<dim>                                     &quadrature,
                const typename Mapping<dim,spacedim>::InternalDataBase    &internal_data,
                internal::FEValues::MappingRelatedData<dim,spacedim>      &output_data) const
{
  // call the function of the base class, but ignoring
  // any potentially detected cell similarity between
  // the current and the previous cell
  MappingQ<dim,spacedim>::fill_fe_values (cell,
                                          CellSimilarity::invalid_next_cell,
                                          quadrature,
                                          internal_data,
                                          output_data);
  // also return the updated flag since any detected
  // similarity wasn't based on the mapped field, but
  // the original vertices which are meaningless
  return CellSimilarity::invalid_next_cell;
}

DEAL_II_NAMESPACE_CLOSE


#endif // dealii__mapping_q_eulerian_h
