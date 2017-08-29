#include <iostream>
#include <fstream>
#include <ostream>
#include <istream>

using namespace std;

//#define SENZA

class TPoint;
typedef class TPoint *PTPoint;
typedef class TPoint &RTPoint;

class TPoint
{
  public:
   
    int PID;
    double x,y,z;
    double Error;

    TPoint( double xi=0, double yi=0, double zi = 0 )
      : PID(-1), x(xi), y(yi), z(zi), Error(0.0){};
  
#ifndef SENZA    
    friend ostream& operator<< ( ostream&, RTPoint p );
    friend istream& operator>> ( istream&, RTPoint p );    
#endif
};

#ifndef SENZA    

ostream& operator<< ( ostream& os, RTPoint p )
{
  os << p.x << " " << p.y << " " << p.z;
  return(os);
}

istream& operator>> ( istream& is, RTPoint p )
{
  is >> p.x >> p.y >> p.z;
  return(is);
}

#endif
