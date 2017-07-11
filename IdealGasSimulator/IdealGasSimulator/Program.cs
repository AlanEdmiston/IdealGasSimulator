using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Data.SqlClient;
using System.Threading;
using System.Data.SqlClient;

namespace IdealGasSimulator
{
    class Program
    {
        int atomsNo = new int();
        int molecules2No = new int();
        double timeIncrement, time, runTime;
        Dvector size = new Dvector();
        List<Atom> atoms;
        List<Molecule2> molecules2s;
        string str;
        Random rand = new Random();

        static void Main(string[] args)
        {
            Program p = new Program();
            p.Process();
        }

        void Process()
        {

            for (time = 0; time < runTime; time += timeIncrement)
            {
                foreach (Atom atom in atoms)
                {
                    System.Threading.ThreadPool.QueueUserWorkItem(new WaitCallback(AtomComputerDelegate), atom);
                }
                foreach (Molecule2 molecule in molecules2s)
                {
                    System.Threading.ThreadPool.QueueUserWorkItem(new WaitCallback(Molecule2ComputerDelegate), molecule);
                }
            }

        }

        void AtomComputer(Atom atom)
        {
            if (atom.isInMolecule == false)
            {
                atom.position = AddVector(atom.position, ScalarMultiple(timeIncrement, atom.velocity));
            }
            if (atom.hasCollided == false)
            {
                if (atom.position.X > size.X || atom.position.X < 0)
                {
                    atom.velocity.X = -atom.velocity.X;
                    atom.hasCollided = true;
                }
                if (atom.position.X > size.X || atom.position.X < 0)
                {
                    atom.velocity.X = -atom.velocity.X;
                    atom.hasCollided = true;
                }
                if (atom.position.Z > size.Z || atom.position.Z < 0)
                {
                    atom.velocity.Z = -atom.velocity.Z;
                    atom.hasCollided = true;
                }

                foreach (Atom atom2 in atoms)
                {
                    if (atom != atom2 && AbsVect(AddVector(atom.position, ScalarMultiple(-1, atom2.position))) < atom.radius + atom2.radius)
                    {
                        Collision(atom, atom2);

                        atom.hasCollided = true;
                    }
                }
            }
            else
            {
                atom.hasCollided = false;
            }

        }

        void Molecule2Computer(Molecule2 molecule)
        {
            molecule.position = AddVector(molecule.position, ScalarMultiple(timeIncrement, molecule.velocity));
            double reducedMass;
            reducedMass = molecule.atom1.mass * molecule.atom2.mass / molecule.mass;
            Dvector r = ScalarMultiple(molecule.atom1.mass / molecule.mass , AddVector(molecule.atom1.position, ScalarMultiple(-1, molecule.atom2.position)));

            molecule.interiaTensor[0, 0] = (Math.Pow(molecule.position.Y - molecule.atom1.position.Y, 2) + Math.Pow(molecule.position.Z - molecule.atom1.position.Z, 2)) * molecule.atom1.mass
                + (Math.Pow(molecule.position.Y - molecule.atom2.position.Y, 2) + Math.Pow(molecule.position.Z - molecule.atom2.position.Z, 2)) * molecule.atom2.mass;
            molecule.interiaTensor[1, 1] = (Math.Pow(molecule.position.X - molecule.atom1.position.X, 2) + Math.Pow(molecule.position.Z - molecule.atom1.position.Z, 2)) * molecule.atom1.mass
                + (Math.Pow(molecule.position.X - molecule.atom2.position.X, 2) + Math.Pow(molecule.position.Z - molecule.atom2.position.Z, 2)) * molecule.atom2.mass;
            molecule.interiaTensor[2, 2] = (Math.Pow(molecule.position.X - molecule.atom1.position.X, 2) + Math.Pow(molecule.position.Y - molecule.atom1.position.Y, 2)) * molecule.atom1.mass
                + (Math.Pow(molecule.position.X - molecule.atom2.position.X, 2) + Math.Pow(molecule.position.Y - molecule.atom2.position.Y, 2)) * molecule.atom2.mass;

            molecule.interiaTensor[1, 0] = (molecule.position.X - molecule.atom1.position.X) * (molecule.atom1.position.Y - molecule.position.Y) * molecule.atom1.mass + (molecule.position.X - molecule.atom2.position.X) * (molecule.atom2.position.Y - molecule.position.Y) * molecule.atom2.mass;
            molecule.interiaTensor[2, 0] = (molecule.position.X - molecule.atom1.position.X) * (molecule.atom1.position.Z - molecule.position.Z) * molecule.atom1.mass + (molecule.position.X - molecule.atom2.position.X) * (molecule.atom2.position.Z - molecule.position.Z) * molecule.atom2.mass;
            molecule.interiaTensor[2, 1] = (molecule.position.Y - molecule.atom1.position.Y) * (molecule.atom1.position.Z - molecule.position.Z) * molecule.atom1.mass + (molecule.position.Y - molecule.atom2.position.Y) * (molecule.atom2.position.Z - molecule.position.Z) * molecule.atom2.mass;

            molecule.interiaTensor[0, 1] = molecule.interiaTensor[1, 0];
            molecule.interiaTensor[0, 2] = molecule.interiaTensor[2, 0];

            molecule.interiaTensor[1, 2] = molecule.interiaTensor[2, 1];

            if (molecule.atom1.position.X > size.X || molecule.atom1.position.X < 0)
            {
                molecule.atom1.velocity.X = -molecule.atom1.velocity.X;
                molecule.velocity.X += molecule.atom1.velocity.X * 2 * (molecule.atom1.mass) / molecule.mass;

                molecule.amplitude = Math.Abs(2 * molecule.atom1.velocity.X / Math.Sqrt(molecule.bondSpringConstant * reducedMass));
                molecule.angularFrequency = Math.Sqrt(molecule.bondSpringConstant * reducedMass);
                molecule.theta = Math.PI / 2;

                molecule.angularVelocity = Crossproduct(ScalarMultiple(2, molecule.atom1.velocity), AddVector(molecule.atom1.position, ScalarMultiple(-1, molecule.atom2.position)));
            }
            if (molecule.atom1.position.Y > size.Y || molecule.atom1.position.Y < 0)
            {
                molecule.atom1.velocity.Y = -molecule.atom1.velocity.Y;
                molecule.velocity.Y += molecule.atom1.velocity.Y * 2 * (molecule.atom1.mass) / molecule.mass;

                molecule.amplitude = Math.Abs(2 * molecule.atom1.velocity.Y / Math.Sqrt(molecule.bondSpringConstant * reducedMass));
                molecule.angularFrequency = Math.Sqrt(molecule.bondSpringConstant * reducedMass);
                molecule.theta = Math.PI / 2;

                molecule.angularVelocity = Crossproduct(ScalarMultiple(2, molecule.atom1.velocity), AddVector(molecule.atom1.position, ScalarMultiple(-1, molecule.atom2.position)));
            }
            if (molecule.atom1.position.Z > size.Z || molecule.atom1.position.Z < 0)
            {
                molecule.atom1.velocity.Z = -molecule.atom1.velocity.Z;
                molecule.velocity.Z += molecule.atom1.velocity.Z * 2 * (molecule.atom1.mass) / molecule.mass;

                molecule.amplitude = Math.Abs(2 * molecule.atom1.velocity.Z / Math.Sqrt(molecule.bondSpringConstant * reducedMass));
                molecule.angularFrequency = Math.Sqrt(molecule.bondSpringConstant * reducedMass);
                molecule.theta = Math.PI / 2;

                molecule.angularVelocity = Crossproduct(ScalarMultiple(2, molecule.atom1.velocity), AddVector(molecule.atom1.position, ScalarMultiple(-1, molecule.atom2.position)));
            }

            if (molecule.atom2.position.X > size.X || molecule.atom2.position.X < 0)
            {
                molecule.atom2.velocity.X = -molecule.atom2.velocity.X;
                molecule.velocity.X += molecule.atom2.velocity.X * 2 * (molecule.atom2.mass) / molecule.mass;

                molecule.amplitude = Math.Abs(2 * molecule.atom2.velocity.X / Math.Sqrt(molecule.bondSpringConstant * reducedMass));
                molecule.angularFrequency = Math.Sqrt(molecule.bondSpringConstant * reducedMass);
                molecule.theta = Math.PI / 2;

                molecule.angularVelocity = Crossproduct(ScalarMultiple(2, molecule.atom1.velocity), AddVector(molecule.atom1.position, ScalarMultiple(-1, molecule.atom2.position)));
            }
            if (molecule.atom2.position.Y > size.Y || molecule.atom2.position.Y < 0)
            {
                molecule.atom2.velocity.Y = -molecule.atom2.velocity.Y;
                molecule.velocity.Y += molecule.atom2.velocity.Y * 2 * (molecule.atom2.mass) / molecule.mass;

                molecule.amplitude = Math.Abs(2 * molecule.atom2.velocity.Y / Math.Sqrt(molecule.bondSpringConstant * reducedMass));
                molecule.angularFrequency = Math.Sqrt(molecule.bondSpringConstant * reducedMass);
                molecule.theta = Math.PI / 2;

                molecule.angularVelocity = Crossproduct(ScalarMultiple(2, molecule.atom1.velocity), AddVector(molecule.atom1.position, ScalarMultiple(-1, molecule.atom2.position)));
            }
            if (molecule.atom2.position.Z > size.Z || molecule.atom2.position.Z < 0)
            {
                molecule.atom2.velocity.Z = -molecule.atom2.velocity.Z;
                molecule.velocity.Z += molecule.atom2.velocity.Z * 2 * (molecule.atom2.mass) / molecule.mass;

                molecule.amplitude = Math.Abs(2 * molecule.atom2.velocity.Z / Math.Sqrt(molecule.bondSpringConstant * reducedMass));
                molecule.angularFrequency = Math.Sqrt(molecule.bondSpringConstant * reducedMass);
                molecule.theta = Math.PI / 2;

                molecule.angularVelocity = Crossproduct(ScalarMultiple(2, molecule.atom1.velocity), AddVector(molecule.atom1.position, ScalarMultiple(-1, molecule.atom2.position)));
            }

            foreach (Atom atom in atoms)
            {
                if(molecule.atom1 != atom && AbsVect(AddVector(atom.position, ScalarMultiple(-1, molecule.atom1.position))) < atom.radius + molecule.atom1.radius)
                {
                    Collision(molecule.atom1, atom);
                    molecule.angularVelocity = Crossproduct(AddVector(molecule.atom1.velocity, ScalarMultiple(-1, molecule.atom2.velocity)), AddVector(molecule.atom1.position, ScalarMultiple(-1, molecule.atom2.position)));
                }

                if (molecule.atom2 != atom && AbsVect(AddVector(atom.position, ScalarMultiple(-1, molecule.atom2.position))) < atom.radius + molecule.atom2.radius)
                {
                    Collision(molecule.atom2, atom);
                    molecule.angularVelocity = Crossproduct(AddVector(molecule.atom1.velocity, ScalarMultiple(-1, molecule.atom2.velocity)), AddVector(molecule.atom1.position, ScalarMultiple(-1, molecule.atom2.position)));
                }
            }

            molecule.angularVelocity = Crossproduct(AddVector(molecule.atom2.velocity, ScalarMultiple(-1, molecule.atom1.velocity)), AddVector(molecule.atom2.velocity, ScalarMultiple(-1, molecule.atom1.velocity)));
            molecule.AngularMomentum.X = molecule.interiaTensor[0, 0] * molecule.angularVelocity.X + molecule.interiaTensor[0, 1] * molecule.angularVelocity.Y + molecule.interiaTensor[0, 2] * molecule.angularVelocity.Z;
            molecule.AngularMomentum.Y = molecule.interiaTensor[1, 0] * molecule.angularVelocity.X + molecule.interiaTensor[1, 1] * molecule.angularVelocity.Y + molecule.interiaTensor[1, 2] * molecule.angularVelocity.Z;
            molecule.AngularMomentum.Z = molecule.interiaTensor[2, 0] * molecule.angularVelocity.X + molecule.interiaTensor[2, 1] * molecule.angularVelocity.Y + molecule.interiaTensor[2, 2] * molecule.angularVelocity.Z;
        }

        void Collision(Atom atom, Atom atom2)
        {
            double radialV1, radialV2;
            Dvector position = new Dvector();
            Dvector rad1UnitV = new Dvector();
            position = AddVector(atom2.position, ScalarMultiple(-1, atom.position));
            radialV1 = DotProduct(atom.velocity, position) / AbsVect(position);
            radialV2 = DotProduct(atom2.velocity, ScalarMultiple(-1, position)) / AbsVect(position);
            rad1UnitV = ScalarMultiple(1 / AbsVect(position), AddVector(atom2.position, ScalarMultiple(-1, atom.position)));
            Dvector momentum = new Dvector();
            momentum = AddVector(ScalarMultiple(atom.mass, atom.velocity), ScalarMultiple(atom2.mass, atom2.velocity));
            double deltaV1 = (radialV2 - 2 * radialV1) * (atom.mass + 1) - radialV2 * atom2.mass;
            double deltaV2 = (radialV2 - 2 * radialV1) * atom.mass - radialV2 * (atom2.mass + 1);

            atom.velocity = AddVector(atom.velocity, ScalarMultiple(-deltaV1, rad1UnitV));
        }

        void AtomComputerDelegate(object atomDelegate)
        {
            Atom atom = (Atom)atomDelegate;
            AtomComputer(atom);
        }
        void Molecule2ComputerDelegate(object Molecule2Delegate)
        {
            Molecule2 molecule = (Molecule2)Molecule2Delegate;
            Molecule2Computer(molecule);
        }

        Atom AtomCreator(int catalogue, double energy, double mass, double radius)
        {
            Atom atom = new Atom();
            atom.catalogue = catalogue;
            atom.mass = mass;
            atom.radius = radius;
            atom.position.X = rand.Next((int)radius, (int)(size.X - radius));
            atom.position.Y = rand.Next((int)radius, (int)(size.Y - radius));
            atom.position.Z = rand.Next((int)radius, (int)(size.Z - radius));

            atom.velocity.X = Math.Sqrt(2 * energy / (mass * 3)) * rand.Next(-100, 101) / 100;
            atom.velocity.X = Math.Sqrt(2 * energy / (mass * 3)) * rand.Next(-100, 101) / 100;
            atom.velocity.X = Math.Sqrt(2 * energy / (mass * 3)) * rand.Next(-100, 101) / 100;

            return atom;
        }

        Molecule2 molecule2Creator(Atom atom1, Atom atom2, double bondLength)
        {
            Molecule2 molecule = new Molecule2();
            molecule.bondLength = bondLength;
            atom2.position = atom1.position;
            atom2.position.X += bondLength;
            if (atom2.position.X > size.X)
            {
                atom2.position.X -= 2 * bondLength;
            }
            return molecule;
        }

        #region //vector algebra methods
        double AbsVect(Dvector vector)
        {
            double absvect = Math.Pow(vector.X, 2) + Math.Pow(vector.Y, 2) + Math.Pow(vector.Z, 2) + Math.Pow(vector.W, 2);
            return absvect;
        }
        Dvector AddVector(Dvector vector1, Dvector vector2)
        {
            Dvector SumVector = new Dvector();
            SumVector.X = vector1.X + vector2.X;
            SumVector.Y = vector1.Y + vector2.Y;
            SumVector.Z = vector1.Z + vector2.Z;
            SumVector.W = vector1.W + vector2.W;
            return SumVector;
        }
        Dvector Crossproduct(Dvector vector1, Dvector vector2)
        {
            Dvector crossProduct = new Dvector();
            crossProduct.X = vector1.Y * vector2.Z - vector2.Y * vector1.Z;
            crossProduct.Y = vector1.Z * vector2.X - vector2.Z * vector1.X;
            crossProduct.Z = vector1.X * vector2.Y - vector2.X * vector1.Y;
            return crossProduct;
        }
        double DotProduct(Dvector vector1, Dvector vector2)
        {
            double dotProduct;
            dotProduct = vector1.X * vector2.X + vector1.Y * vector2.Y + vector1.Z * vector2.Z + vector1.W * vector2.W;
            return dotProduct;
        }
        Dvector ScalarMultiple(double scalar, Dvector vector)
        {
            Dvector scalarMultiple = new Dvector();
            scalarMultiple.X = scalar * vector.X;
            scalarMultiple.Y = scalar * vector.Y;
            scalarMultiple.Z = scalar * vector.Z;
            scalarMultiple.W = scalar * vector.W;
            return scalarMultiple;
        }
        #endregion
    }
}
