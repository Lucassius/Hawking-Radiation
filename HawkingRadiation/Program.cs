using System;
using System.Collections;
using System.Collections.Generic;


// Implementation of a Genetic Algorithm in C#
// Rastrigin Function: f(x) = A * n + ∑[i=1 to n] (x[i]^2 - A * cos(2 * π * x[i]))

/*The initial parameters used for 
 *starting the genetic algorithm were crossover= 80%, mutation = 5 %,
 *population size= 100, generation size = 2000 and chromosome size = 2.
*/


public delegate double GAFunction(double[] values);

// A Genetic Algorithm class
public class GA
{
    public double MutationRate;
    public double CrossoverRate;
    public int ChromosomeLength;
    public int PopulationSize;
    public int GenerationSize;
    public double TotalFitness;
    public bool Elitism;
    private ArrayList CurrentGenerationList;
    private ArrayList NextGenerationList;
    private ArrayList FitnessList;
    static private Random rand = new Random();
    static private GAFunction getFitness;
    public GAFunction FitnessFunction
    {
        get { return getFitness; }
        set { getFitness = value; }
    }

    // Constructor with user specified crossover rate,
    // mutation rate, population size, generation size
    // and chromosome length.
    public GA(double XoverRate, double mutRate, int popSize, int genSize, int ChromLength)
    {
        Elitism = false;
        MutationRate = mutRate;
        CrossoverRate = XoverRate;
        PopulationSize = popSize;
        GenerationSize = genSize;
        ChromosomeLength = ChromLength;
    }

    // Method which launches the GA into execution mode.
    public void LaunchGA()
    {
        // Create the arrays to hold the fitness,
        // current and next generation lists.
        FitnessList = new ArrayList();
        CurrentGenerationList = new ArrayList(GenerationSize);
        NextGenerationList = new ArrayList(GenerationSize);

        // Initialize the mutation rate.
        Chromosome.ChromosomeMutationRate = MutationRate;

        // Create the initial chromosome population by repeatedly
        // calling the user supplied fitness function.
        for (int i = 0; i < PopulationSize; i++)
        {
            Chromosome g = new Chromosome(ChromosomeLength, true);
            CurrentGenerationList.Add(g);
        }

        // Rank the initial chromosome population.
        RankPopulation();

        // Loop through the entire generation size creating
        // and evaluating generations of new chromosomes.
        for (int i = 0; i < GenerationSize; i++)
        {
            CreateNextGeneration();
            RankPopulation();
        }
    }

    // After ranking all the chromosomes by fitness, use a
    // "roulette wheel" selection method that allocates a large
    // probability of being selected to those chromosomes with the
    // highest fitness.
    private int RouletteSelection()
    {
        double randomFitness = rand.NextDouble() * TotalFitness;
        int idx = -1;
        int mid;
        int first = 0;
        int last = PopulationSize - 1;
        mid = (last - first) / 2;
        while (idx == -1 && first <= last)
        {
            if (randomFitness < (double)FitnessList[mid])
            {
                last = mid;
            }
            else if (randomFitness > (double)FitnessList[mid])
            {
                first = mid;
            }
            mid = (first + last) / 2;
            if ((last - first) == 1) idx = last;
        }
        return idx;
    }

    // Rank population and then sort it in order of fitness.
    private void RankPopulation()
    {
        TotalFitness = 0;
        for (int i = 0; i < PopulationSize; i++)
        {
            Chromosome g = ((Chromosome)CurrentGenerationList[i]);
            g.ChromosomeFitness = FitnessFunction(g.ChromosomeGenes);
            TotalFitness += g.ChromosomeFitness;
        }

        CurrentGenerationList.Sort(new ChromosomeComparer());
        double fitness = 0.0;
        FitnessList.Clear();
        for (int i = 0; i < PopulationSize; i++)
        {
            fitness += ((Chromosome)CurrentGenerationList[i]).ChromosomeFitness;
            FitnessList.Add((double)fitness);
        }
    }

    // Create a new generation of chromosomes.
    private void CreateNextGeneration()
    {
        NextGenerationList.Clear();
        Chromosome g = null;
        if (Elitism)
            g = (Chromosome)CurrentGenerationList[PopulationSize - 1];

        for (int i = 0; i < PopulationSize; i += 2)
        {
            int pidx1 = RouletteSelection();
            int pidx2 = RouletteSelection();
            Chromosome parent1, parent2, child1, child2;
            parent1 = ((Chromosome)CurrentGenerationList[pidx1]);
            parent2 = ((Chromosome)CurrentGenerationList[pidx2]);
            if (rand.NextDouble() < CrossoverRate)
            {
                parent1.Crossover(ref parent2, out child1, out child2);
            }
            else
            {
                child1 = parent1;
                child2 = parent2;
            }
            child1.Mutate();
            child2.Mutate();
            NextGenerationList.Add(child1);
            NextGenerationList.Add(child2);
        }

        if (Elitism && g != null)
            NextGenerationList[0] = g;

        CurrentGenerationList.Clear();
        for (int i = 0; i < PopulationSize; i++)
        {
            CurrentGenerationList.Add(NextGenerationList[i]);
        }
    }

    // Extract the best values based on fitness from the current generation.
    public void GetBestValues(out double[] values, out double fitness)
    {
        Chromosome g = ((Chromosome)CurrentGenerationList[PopulationSize - 1]);
        values = new double[g.ChromosomeLength];
        g.ExtractChromosomeValues(ref values);
        fitness = (double)g.ChromosomeFitness;
    }


    public class Chromosome
    {
        public double[] ChromosomeGenes;
        public int ChromosomeLength;
        public double ChromosomeFitness;
        public static double ChromosomeMutationRate;
        private static Random rand = new Random();

        // Chromosome class constructor.
        public Chromosome(int length, bool createGenes)
        {
            ChromosomeLength = length;
            ChromosomeGenes = new double[length];
            if (createGenes)
            {
                for (int i = 0; i < ChromosomeLength; i++)
                {
                    ChromosomeGenes[i] = rand.NextDouble();
                }
            }
        }

        // Creates two offspring children using a single crossover point.
        public void Crossover(ref Chromosome Chromosome2, out Chromosome child1, out Chromosome child2)
        {
            int position = (int)(rand.NextDouble() * (double)ChromosomeLength);
            child1 = new Chromosome(ChromosomeLength, false);
            child2 = new Chromosome(ChromosomeLength, false);

            for (int i = 0; i < ChromosomeLength; i++)
            {
                if (i < position)
                {
                    child1.ChromosomeGenes[i] = ChromosomeGenes[i];
                    child2.ChromosomeGenes[i] = Chromosome2.ChromosomeGenes[i];
                }
                else
                {
                    child1.ChromosomeGenes[i] = Chromosome2.ChromosomeGenes[i];
                    child2.ChromosomeGenes[i] = ChromosomeGenes[i];
                }
            }
        }

        // Mutates the chromosome genes by randomly switching them around.
        public void Mutate()
        {
            for (int position = 0; position < ChromosomeLength; position++)
            {
                if (rand.NextDouble() < ChromosomeMutationRate)
                {
                    ChromosomeGenes[position] = (ChromosomeGenes[position] + rand.NextDouble()) / 2.0;
                }
            }
        }

        // Extracts the chromosome values.
        public void ExtractChromosomeValues(ref double[] values)
        {
            for (int i = 0; i < ChromosomeLength; i++)
            {
                values[i] = ChromosomeGenes[i];
            }
        }
    }

    // Compares two chromosomes by their fitness values.
    public sealed class ChromosomeComparer : IComparer
    {
        public int Compare(object x, object y)
        {
            if (!(x is Chromosome) || !(y is Chromosome))
            {
                throw new ArgumentException("Not of type Chromosome");
            }

            if (((Chromosome)x).ChromosomeFitness > ((Chromosome)y).ChromosomeFitness)
            {
                return 1;
            }
            else if (((Chromosome)x).ChromosomeFitness == ((Chromosome)y).ChromosomeFitness)
            {
                return 0;
            }
            else
            {
                return -1;
            }
        }
    }

    public static double GenAlgTestFcn(double[] values)
    {
        if(values.Length != 2)
        {
            throw new Exception("Should only have 1 argument");
        }

        double M = values[0];
        double temperature = (1 / (8 * Math.PI * M));

        return temperature;
    }


    static void GeneticAlgorithmTest()
    {
        GA ga = new GA(0.8, 0.05, 100, 2000, 2);
        ga.FitnessFunction = new GAFunction(GenAlgTestFcn);
        ga.Elitism = true;
        ga.LaunchGA();

        double[] values;
        double fitness;
        ga.GetBestValues(out values, out fitness);

        Console.WriteLine("Calculated value is:\nM = {0}\n", values[0]);
        Console.WriteLine("Temperature at M = {0} is: T = {1}", values[0], fitness);
        Console.WriteLine("\nPress ENTER to terminate the program");
        Console.ReadLine();
    }

    static void Main(string[] args)
    {
        GeneticAlgorithmTest();

    }


}