using System;
using System.Collections;
using System.Collections.Generic;
using System.Text;
using System.Drawing;
using System.Timers;
using System.Diagnostics;

namespace TSP {
    class ProblemAndSolver {
        private class TSPSolution {
            /// <summary>
            /// we use the representation [cityB,cityA,cityC] 
            /// to mean that cityB is the first city in the solution, cityA is the second, cityC is the third 
            /// and the edge from cityC to cityB is the final edge in the path.  
            /// you are, of course, free to use a different representation if it would be more convenient or efficient 
            /// for your node data structure and search algorithm. 
            /// </summary>
            public ArrayList
                Route;

            public TSPSolution(ArrayList iroute) {
                Route = new ArrayList(iroute);
            }


            /// <summary>
            ///  compute the cost of the current route.  does not check that the route is complete, btw.
            /// assumes that the route passes from the last city back to the first city. 
            /// </summary>
            /// <returns></returns>
            public double costOfRoute() {
                // go through each edge in the route and add up the cost. 
                int x;
                City here;
                double cost = 0D;

                for (x = 0; x < Route.Count - 1; x++) {
                    here = Route[x] as City;
                    cost += here.costToGetTo(Route[x + 1] as City);
                }
                // go from the last city to the first. 
                here = Route[Route.Count - 1] as City;
                cost += here.costToGetTo(Route[0] as City);
                return cost;
            }
        }

        #region private members

        private const int DEFAULT_SIZE = 25;

        private const int CITY_ICON_SIZE = 5;

        /// <summary>
        /// the cities in the current problem.
        /// </summary>
        private City[] Cities;
        /// <summary>
        /// a route through the current problem, useful as a temporary variable. 
        /// </summary>
        private ArrayList Route;
        /// <summary>
        /// best solution so far. 
        /// </summary>
        private TSPSolution bssf;

        /// <summary>
        /// how to color various things. 
        /// </summary>
        private Brush cityBrushStartStyle;
        private Brush cityBrushStyle;
        private Pen routePenStyle;


        /// <summary>
        /// keep track of the seed value so that the same sequence of problems can be 
        /// regenerated next time the generator is run. 
        /// </summary>
        private int _seed;
        /// <summary>
        /// number of cities to include in a problem. 
        /// </summary>
        private int _size;

        /// <summary>
        /// random number generator. 
        /// </summary>
        private Random rnd;
        #endregion

        #region public members.
        public int Size {
            get { return _size; }
        }

        public int Seed {
            get { return _seed; }
        }
        #endregion

        public const int DEFAULT_SEED = -1;

        #region Constructors
        public ProblemAndSolver() {
            initialize(DEFAULT_SEED, DEFAULT_SIZE);
        }

        public ProblemAndSolver(int seed) {
            initialize(seed, DEFAULT_SIZE);
        }

        public ProblemAndSolver(int seed, int size) {
            initialize(seed, size);
        }
        #endregion

        #region Private Methods

        /// <summary>
        /// reset the problem instance. 
        /// </summary>
        private void resetData() {
            Cities = new City[_size];
            Route = new ArrayList(_size);
            bssf = null;

            for (int i = 0; i < _size; i++)
                Cities[i] = new City(rnd.NextDouble(), rnd.NextDouble());

            cityBrushStyle = new SolidBrush(Color.Black);
            cityBrushStartStyle = new SolidBrush(Color.Red);
            routePenStyle = new Pen(Color.LightGray, 1);
            routePenStyle.DashStyle = System.Drawing.Drawing2D.DashStyle.Solid;
        }

        private void initialize(int seed, int size) {
            this._seed = seed;
            this._size = size;
            if (seed != DEFAULT_SEED)
                this.rnd = new Random(seed);
            else
                this.rnd = new Random();
            this.resetData();
        }

        #endregion

        #region Public Methods

        /// <summary>
        /// make a new problem with the given size.
        /// </summary>
        /// <param name="size">number of cities</param>
        public void GenerateProblem(int size) {
            this._size = size;
            resetData();
        }

        /// <summary>
        /// return a copy of the cities in this problem. 
        /// </summary>
        /// <returns>array of cities</returns>
        public City[] GetCities() {
            City[] retCities = new City[Cities.Length];
            Array.Copy(Cities, retCities, Cities.Length);
            return retCities;
        }

        /// <summary>
        /// draw the cities in the problem.  if the bssf member is defined, then
        /// draw that too. 
        /// </summary>
        /// <param name="g">where to draw the stuff</param>
        public void Draw(Graphics g) {
            float width = g.VisibleClipBounds.Width - 45F;
            float height = g.VisibleClipBounds.Height - 15F;
            Font labelFont = new Font("Arial", 10);

            g.DrawString("n(c) means this node is the nth node in the current solution and incurs cost c to travel to the next node.", labelFont, cityBrushStartStyle, new PointF(0F, 0F));

            // Draw lines
            if (bssf != null) {
                // make a list of points. 
                Point[] ps = new Point[bssf.Route.Count];
                int index = 0;
                foreach (City c in bssf.Route) {
                    if (index < bssf.Route.Count - 1)
                        g.DrawString(" " + index + "(" + c.costToGetTo(bssf.Route[index + 1] as City) + ")", labelFont, cityBrushStartStyle, new PointF((float)c.X * width + 3F, (float)c.Y * height));
                    else
                        g.DrawString(" " + index + "(" + c.costToGetTo(bssf.Route[0] as City) + ")", labelFont, cityBrushStartStyle, new PointF((float)c.X * width + 3F, (float)c.Y * height));
                    ps[index++] = new Point((int)(c.X * width) + CITY_ICON_SIZE / 2, (int)(c.Y * height) + CITY_ICON_SIZE / 2);
                }

                if (ps.Length > 0) {
                    g.DrawLines(routePenStyle, ps);
                    g.FillEllipse(cityBrushStartStyle, (float)Cities[0].X * width - 1, (float)Cities[0].Y * height - 1, CITY_ICON_SIZE + 2, CITY_ICON_SIZE + 2);
                }

                // draw the last line. 
                g.DrawLine(routePenStyle, ps[0], ps[ps.Length - 1]);
            }

            // Draw city dots
            foreach (City c in Cities) {
                g.FillEllipse(cityBrushStyle, (float)c.X * width, (float)c.Y * height, CITY_ICON_SIZE, CITY_ICON_SIZE);
            }

        }

        /// <summary>
        ///  return the cost of the best solution so far. 
        /// </summary>
        /// <returns></returns>
        public double costOfBssf() {
            if (bssf != null)
                return (bssf.costOfRoute());
            else
                return -1D;
        }
        #endregion


        #region My Methods
        /// <summary>
        ///  solve the problem.  This is the entry point for the solver when the run button is clicked
        /// right now it just picks a simple solution. 
        /// </summary>
        public void solveProblem() {

            //GET AN INITIAL BSSF

            //Nearest Neighbor greedy algorithm to get the baseline function
            Route = new ArrayList();
            double minCost;
            int tempCityIndex = 0;
            int currentCityIndex = 0;
            ArrayList bssfIndexes = new ArrayList();
            bssfIndexes.Add(0);
            Route.Add(Cities[currentCityIndex]);
            while (Route.Count < Cities.Length) {
                minCost = double.PositiveInfinity;
                for (int i = 0; i < Cities.Length; i++) {
                    if (i != currentCityIndex && !Route.Contains(Cities[i])) {
                        double cost = Cities[currentCityIndex].costToGetTo(Cities[i]);
                        if (cost < minCost) {
                            minCost = cost;
                            tempCityIndex = i;
                        }
                    }
                }
                currentCityIndex = tempCityIndex;
                Route.Add(Cities[currentCityIndex]);
                bssfIndexes.Add(currentCityIndex);

            }
            bssf = new TSPSolution(Route);
            double bssfCost = bssf.costOfRoute();

            // END OF INITIAL BSSF

            //create matrix
            //Get an average cost of all the paths
            double totalCosts = 0;
            double[,] matrix = new double[Cities.Length, Cities.Length];
            for (int i = 0; i < Cities.Length; i++) {
                matrix[i, i] = double.PositiveInfinity;
                for (int j = i + 1; j < Cities.Length; j++) {
                    double cost = Cities[i].costToGetTo(Cities[j]);
                    matrix[i, j] = cost;
                    matrix[j, i] = cost;
                    totalCosts += cost;
                }
            }
            double averageCost = totalCosts / (Cities.Length * Cities.Length * 1.0) / 2.0;

            //Another try to take depth into account
            double[] minPerRow = new double[Cities.Length];
            for (int i = 0; i < Cities.Length; i++) {
                double bestSoFar = Double.PositiveInfinity;
                for (int j = 0; j < Cities.Length; j++) {
                    if (i != j) {
                        double cost = Cities[i].costToGetTo(Cities[j]);
                        if (cost < bestSoFar) {
                            bestSoFar = cost;
                        }
                    }
                }
                minPerRow[i] = bestSoFar;
            }

            //Get the Children possible from the initial state (all of them)
            ArrayList childrenIndexes = new ArrayList();
            for (int i = 1; i < Cities.Length; i++) {
                childrenIndexes.Add(i);
            }
            //Turn that into the initial state
            State s = new State(0, childrenIndexes, (double[,])matrix.Clone(), new ArrayList(), 0, 0, Cities, 0);
            //Get the Bound on the newly Made initial State
            double bound = s.lowerBound;
            //Put the state on the new Queue
            PriorityQueueDemo.PriorityQueue<double, State> q = new PriorityQueueDemo.PriorityQueue<double, State>();
            q.Enqueue(bound, s);

            //Kick off the StopWatch
            Stopwatch watch = new Stopwatch();
            watch.Start();

            int maxPqSize = 0;

            //Branch and Bound It Up
            //Cut off if we go past a minute
            while (!q.IsEmpty && bssfCost > bound && watch.Elapsed.Minutes < 1) {
                if (q.Count > maxPqSize)
                    maxPqSize = q.Count;

                State u = q.Dequeue().Value;

                //If its possible to be a better state
                if (u.lowerBound < bssfCost) {
                    ArrayList childStates = u.generateChildrenStates();
                    foreach (State w in childStates) {
                        if (w.lowerBound < bssfCost) {
                            if (w.solutionFound && w.currentCost < bssfCost) {
                                bssfCost = w.currentCost;
                                bssfIndexes = w.pathSoFar;
                            }
                            else {
                                double pVal = w.lowerBound;
                                foreach (int i in w.childrenIndexes) {
                                    pVal += minPerRow[i];
                                }
                                q.Enqueue(pVal, w);
                                //Other Attempts at getting a Better Priority Value based on Depth
                                //q.Enqueue(w.lowerBound, w); //36
                                //q.Enqueue(w.lowerBound-(w.depthIntoSolution*averageCost), w); // 29
                            }
                        }

                    }
                }
            }

            Route = new ArrayList();
            foreach (int i in bssfIndexes) {
                Route.Add(Cities[i]);
            }

            bssf = new TSPSolution(Route);

            System.Diagnostics.Debug.WriteLine(maxPqSize);
            // update the cost of the tour. 
            Program.MainForm.tbCostOfTour.Text = " " + bssf.costOfRoute();
            //Output the Time it took
            Program.MainForm.tbElapsedTime.Text = watch.Elapsed.Minutes + ":" + watch.Elapsed.Seconds;
            // do a refresh. 
            Program.MainForm.Invalidate();
        }

        private void printMatrix(double[,] matrix) {
            for (int i = 0; i <= matrix.GetUpperBound(0); i++) {
                for (int j = 0; j <= matrix.GetUpperBound(1); j++) {
                    System.Diagnostics.Debug.Write(matrix[i, j] + " ");
                }
                System.Diagnostics.Debug.WriteLine("");
            }
        }


        #endregion


        #region State Object

        private class State {
            public double lowerBound;
            public ArrayList childrenIndexes;
            public double[,] costArray;
            public ArrayList pathSoFar;
            public double currentCost;
            public int currentIndex;
            public City[] cities;
            public int depthIntoSolution;
            public bool solutionFound;

            public State(double b, ArrayList c, double[,] costs, ArrayList path, double startCost, int index, City[] cit, int d) {
                lowerBound = b;
                childrenIndexes = c;
                costArray = costs;
                pathSoFar = path;
                solutionFound = false;
                currentCost = startCost;
                currentIndex = index;
                cities = cit;
                depthIntoSolution = d;

                //If we have no children, we are done. 
                if (childrenIndexes.Count == 0) {
                    solutionFound = true;
                }

                //generate the lower bound
                this.lowerBound = b + reducedCostMatrix();

                if (solutionFound) {
                    currentCost += cit[0].costToGetTo(cit[index]);
                    pathSoFar.Add(index);
                }
            }

            private double reducedCostMatrix() {
                double[,] reducedMatrix = costArray;
                double costSoFar = 0;
                double noChangeRows = 0;

                //reduce rows
                ArrayList rowsToReduce = (ArrayList)childrenIndexes.Clone();
                rowsToReduce.Add(currentIndex);
                foreach (int i in rowsToReduce) {
                    double minCost = double.PositiveInfinity;
                    foreach (int j in childrenIndexes) {
                        if (reducedMatrix[i, j] < minCost) {
                            minCost = reducedMatrix[i, j];
                        }
                    }

                    if (!double.IsPositiveInfinity(minCost)) {
                        foreach (int j in childrenIndexes) {
                            reducedMatrix[i, j] -= minCost;
                        }
                        costSoFar += minCost;
                    }
                    else {
                        noChangeRows++;
                    }
                }

                //reduce columns
                foreach (int i in childrenIndexes) {
                    double minCost = double.PositiveInfinity;
                    foreach (int j in rowsToReduce) {
                        if (reducedMatrix[j, i] < minCost) {
                            minCost = reducedMatrix[j, i];
                        }
                    }
                    if (!double.IsPositiveInfinity(minCost)) {
                        foreach (int j in rowsToReduce) {
                            reducedMatrix[j, i] -= minCost;
                        }
                    }
                    else {
                        noChangeRows++;
                    }
                }

                if (noChangeRows >= reducedMatrix.GetUpperBound(0)) {
                    solutionFound = true;
                    costSoFar += cities[1].costToGetTo(cities[currentIndex]);
                }

                if (solutionFound) {
                    return 0;
                }

                return costSoFar;
            }

            public ArrayList generateChildrenStates() {
                ArrayList childStates = new ArrayList();
                //For every possible child index, create a state
                foreach (int i in childrenIndexes) {
                    //Make New Children
                    ArrayList newChildren = (ArrayList)childrenIndexes.Clone();
                    newChildren.Remove(i);
                    // Copy and update the Path
                    ArrayList newPath = (ArrayList)pathSoFar.Clone();
                    newPath.Add(currentIndex);
                    //Add new City to the Path
                    double cost = cities[currentIndex].costToGetTo(cities[i]);

                    //make new matrix
                    double[,] newCost = (double[,])costArray.Clone();
                    for (int j = 0; j <= newCost.GetUpperBound(0); j++) {
                        newCost[j, currentIndex] = double.PositiveInfinity;
                    }

                    //Add it to the ChildStates
                    childStates.Add(new State(lowerBound + costArray[currentIndex, i], newChildren, newCost, newPath, currentCost + cost, i, cities, depthIntoSolution + 1));
                }
                //Return the childStates
                return childStates;
            }
        }
        #endregion
    }

}
