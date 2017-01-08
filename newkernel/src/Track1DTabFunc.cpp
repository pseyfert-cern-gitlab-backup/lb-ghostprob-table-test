//============================================================================
/** @file Track1DTabFunc.cpp
 *
 *  Implementation file for class : Track::TabulatedFunction1D
 *
 *  @author Paul Seyfert   Paul.Seyfert@cern.ch
 *  @date   2016-12-20
 */
//============================================================================

// STL
#include <algorithm>

// GaudiKernel
#include "GaudiKernel/GaudiException.h"

// local
#include "TrackKernel/Track1DTabFunc.h"

using namespace Track;

//============================================================================

// Destructor
TabulatedFunction1D::~TabulatedFunction1D( ) { }

//============================================================================

// Constructor
TabulatedFunction1D::TabulatedFunction1D( std::initializer_list<float> x ) :
  m_xedges(x) {
    if (m_xedges.size()<2) {
      throw GaudiException( "TabulatedFunction1D() : must be initialized with more than one bin edge",
          "*Track::TabulatedFunction1D*", StatusCode::FAILURE );
    }
    m_width = 1.f/(m_xedges.size() - 1);
    if (!std::is_sorted(m_xedges.begin(),m_xedges.end())) {
      throw GaudiException( "TabulatedFunction1D() : must be initialized with sorted braced initializer list",
          "*Track::TabulatedFunction1D*", StatusCode::FAILURE );
    }
    m_begin = m_xedges.begin();
    m_end = m_xedges.end();
  }

//============================================================================

// evaluation function
float 
TabulatedFunction1D::value( const float x ) const
{
  // out of range check
  if (x<=(*(m_begin))) return 0.f;
  if (x>=(*(m_end-1))) return 1.f;

  // iterator to the first element that is not smaller than x
  // may be end() - if x is larger than the last element
  // may be begin() - if x is smaller than the first element
  auto up = std::lower_bound (m_begin, m_end, x);
  // iterator to the last element that is smaller than x
  // (may be out of range)
  auto low = up - 1;

  // y-value for the lower edge of the x-bin we're in
  float edge = m_width * (low - m_begin);

  // by what fraction did we enter the x bin
  float add  = m_width * (x - (*low))/((*up) - (*low));

  return edge + add;
}

//============================================================================
