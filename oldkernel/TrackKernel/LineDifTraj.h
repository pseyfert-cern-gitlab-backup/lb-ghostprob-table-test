// $Id: LineDifTraj.h,v 1.1 2009-07-08 13:33:45 wouter Exp $
#ifndef TRACKKERNEL_LineDifTraj_H
#define TRACKKERNEL_LineDifTraj_H 1

// Include files
#include <memory>
#include <ostream>
#include "Kernel/DifTraj.h"

namespace LHCb
{

    class LineDifTraj: public LHCb::DifTraj<4> {

    public:

      /// Destructor
      ~LineDifTraj() {}

      // clone thyself...
      std::unique_ptr<Trajectory> clone() const override;

      /// Constructor from the middle point and a direction vector
      LineDifTraj( const Point& middle,
                   const Vector& dir,
                   const Range& range );

      /// Constructor from a begin and an end point
      LineDifTraj( const Point& begPoint,
                   const Point& endPoint );

      /// Point on the trajectory at arclength from the starting point
      Point position( double arclength ) const override;

      /// First derivative of the trajectory at arclength from the starting point
      Vector direction( double arclength=0 ) const override;

      /// Second derivative of the trajectory at arclength from the starting point
      Vector curvature( double arclength=0 ) const override;

      /// Create a parabolic approximation to the trajectory
      /// at arclength from the starting point
      void expansion( double arclength,
                      Point& p,
                      Vector& dp,
                      Vector& ddp ) const override;

      /// Determine the distance in arclenghts to the
      /// closest point on the trajectory to a given point
      double muEstimate( const Point& point ) const override;

      /// Number of arclengths until deviation of the trajectory from the
      /// expansion reaches the given tolerance.
      double distTo1stError( double arclength,
                             double tolerance,
                             int pathDirection = +1 ) const override;

      /// Number of arclengths until deviation of the trajectory from the
      /// expansion reaches the given tolerance.
      double distTo2ndError( double arclength,
                             double tolerance,
                             int pathDirection = +1 ) const override;

      using LHCb::DifTraj<4>::arclength;
      /// Distance, along the Trajectory, between position(mu1) and
      /// position(mu2). Trivial because LineDifTraj is parameterized in
      /// arclength.
      double arclength(double mu1, double mu2) const override { return mu2 - mu1 ; }

      typedef LHCb::DifTraj<4>::Parameters Parameters;
      typedef LHCb::DifTraj<4>::Derivative Derivative;
      Parameters parameters() const override;
      LineDifTraj& operator+=(const Parameters&) override;
      Derivative derivative(double arclen) const override;

      std::ostream& print(std::ostream& s=std::cout) const;

    private:
      Vector m_dir;
      Point  m_pos;

    }; // class LineDifTraj

}
#endif /// LineDifTraj
