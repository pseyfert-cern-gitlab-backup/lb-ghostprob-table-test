//============================================================================
/** @file Track1DTabFunc.h
 *
 *  Header file for utility class : Track::TabulatedFunction1D
 *
 *  @author Paul Seyfert      Paul.Seyfert@cern.ch
 *  @date   2016-12-20
 */
//============================================================================

#ifndef TRACKKERNEL_TRACK1DTABFUNC_H
#define TRACKKERNEL_TRACK1DTABFUNC_H 1

#include <vector>

namespace Track
{

  //============================================================================
  /** @class Track::TabulatedFunction1D Track1DTabFunc.h
   *
   *  A class describing a tabulated function with equidistant y-binning from 0 to 1.
   *
   *  @author Paul Seyfert      Paul.Seyfert@cern.ch
   *  @date   2016-12-20
   */
  //============================================================================

  class TabulatedFunction1D
  {

  public:

    /** Default Constructor with optional interpolator type argument
     *
     *  @param x         braced initializer list for the x bins (N bins = N+1 x values)
     */
    TabulatedFunction1D( std::initializer_list<float> x ) ;

    /// Destructor
    virtual ~TabulatedFunction1D( );

  public:

    /** Computes the function value (y) for the given parameter (x) value
     *  with linear interpolation between bin edges
     *  does out of range check
     *
     *  @param x The parameter value
     *
     *  @return The value of the function at the given parameter value
     */
    float value( const float x ) const ;

  private: // data
    
    // edges of x bins
    const std::vector<float>                 m_xedges;

    // width of y bins
    float                                    m_width;

    // never changes after initialization
    std::vector<float>::const_iterator       m_begin;
    std::vector<float>::const_iterator       m_end;
  };

}

#endif // TRACKKERNEL_TRACK1DTABFUNC_H
