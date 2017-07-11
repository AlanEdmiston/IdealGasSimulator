using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace IdealGasSimulator
{
    class Molecule2
    {
        public Dvector position, velocity, acceleration, angularVelocity, AngularMomentum;
        public double mass, bondLength, bondSpringConstant, angularFrequency, amplitude, theta;
        public Atom atom1, atom2;
        public double[,] interiaTensor = new double[3, 3];
    }
}
