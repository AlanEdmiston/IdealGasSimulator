using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace IdealGasSimulator
{
    class Atom
    {
        public Dvector position, velocity, acceleration;
        public double mass, radius;
        public bool hasCollided, isInMolecule;
        public int catalogue;
    }
}
