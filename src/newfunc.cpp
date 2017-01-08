//============================================================================
/** @file Track1DTabFunc.cpp
 *
 *  Implementation file for class : Track::NewFunc
 *
 *  @author Paul Seyfert   Paul.Seyfert@cern.ch
 *  @date   2016-12-20
 */
//============================================================================

// STL
#include <algorithm>

// local
#include "newfunc.h"

using namespace Track;

//============================================================================

// Destructor
NewFunc::~NewFunc( ) { }

//============================================================================

// Constructor
NewFunc::NewFunc( std::initializer_list<float> x ) :
  m_xedges(x) {
    if (m_xedges.size()<2) {
      throw 11;
    }
    m_width = 1.f/(m_xedges.size() - 1);
    if (!std::is_sorted(m_xedges.begin(),m_xedges.end())) {
      throw 11;
    }
    m_begin = m_xedges.begin();
    m_end = m_xedges.end();
  }

//============================================================================

// evaluation function
float 
NewFunc::value( const float x ) const
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
