import simplex.Simplex;

import java.util.Arrays;
import java.util.Locale;
import java.util.Scanner;

public class Main {
    public static void main(String[] args) {
        try {
            Scanner input = new Scanner(System.in);
            input.useLocale(Locale.UK);
            System.out.printf("Input format:\n" +
                    "N - number of variables in the Objective Function\n" +
                    "[C] Nx1 - vertical vector of OF coefficients\n" +
                    "n - number of constraints\n" +
                    "[A] nxN - matrix of CF coefficients\n"+
                    "[x] Nx1 - vertical vector of initial solution\n"+
                    "[b] nx1 - vertical vector of RHS numbers\n" +
                    "e - approximation accuracy\n");
            int N = input.nextInt(); //Number of variables in the objective function
            Matrix C = new Matrix(N, 1); //Vector of coefficients of the OF
            C.readMatrix(input);
            int n = input.nextInt(); //Number of constraints
            Matrix A = new Matrix(n, N); //Matrix of coefficients of the CF
            A.readMatrix(input);
            Matrix D = new Matrix(N, N); //Matrix of initial solution
            for (int i = 0; i < N; i++) {
                D.setValue(i,i,input.nextDouble());
            }
            Matrix b = new Matrix(n, 1); //Vector of RHS numbers
            b.readMatrix(input);
            double eps = input.nextDouble();//Approximation accuracy
            double alpha1 = 0.5;
            double alpha2 = 0.9;
            /*
            System.out.println("-----");
            System.out.println("[A]");
            A.printMatrix();
            System.out.println("[D]");
            D.printMatrix();
            System.out.println("[b]");
            b.printMatrix();
            System.out.println("[C]");
            C.printMatrix(); */
            Matrix x1 = computeIterations(C,A,D,alpha1,eps);
            Matrix x2 = computeIterations(C,A,D,alpha2,eps);
            System.out.println("Final iteration of decision variables vector x for alpha = 0.5:");
            x1.printMatrix();
            double solution1 = 0;
            for (int i = 0; i < N; i++) {
                solution1 += C.getValue(i,0) * x1.getValue(i,0);
            }
            System.out.println("Value of the objective function: " + solution1);

            System.out.println("Final iteration of decision variables vector x for alpha = 0.9:");
            x2.printMatrix();
            double solution2 = 0;
            for (int i = 0; i < N; i++) {
                solution2 += C.getValue(i,0) * x2.getValue(i,0);
            }
            System.out.println("Value of the objective function: " + solution2);


            Simplex simplex = new Simplex(
                    Arrays.stream(C.matrix)
                            .flatMap(row -> Arrays.stream(row)
                                    .mapToObj(value -> (float) value))
                            .toList(),
                    Arrays.stream(A.matrix)
                            .map(row -> Arrays.stream(row)
                                    .mapToObj(value -> (float) value) // Преобразование double в Float
                                    .toList()
                            )
                            .toList(),
                    Arrays.stream(b.matrix)
                            .flatMap(row -> Arrays.stream(row)
                                    .mapToObj(value -> (float) value))
                            .toList());
            simplex.apply(5);
            System.out.println("Simplex result: " + simplex.getAnswer());
        } catch (IllegalArgumentException e) {
            System.out.println("\n" + e.getMessage());
        }

    }
    public static Matrix computeIterations(Matrix C, Matrix A, Matrix D, double alpha, double eps) {
        //System.out.printf("Iterations of interior point, alpha = %f", alpha);
        double delta = 1;
        int N = D.getRows();
        Matrix x = interiorPointMethod(C, A, D, alpha);
        int iteration = 0;
        while (delta > eps) {
            //System.out.println("Iteration: " + iteration);
            //x.printMatrix();
            Matrix newD = new Matrix(N, N);
            for (int i = 0; i < N; i++) {
                newD.setValue(i,i,x.getValue(i,0));
            }
            Matrix newx = interiorPointMethod(C, A, newD, alpha);
            delta = Math.abs(computeNorm(newx) - computeNorm(x));
            x = newx;
            iteration++;
        }
        return x;
    }
    public static double computeNorm(Matrix V) {
        if (V.getCols() > 1) {
            throw new IllegalArgumentException("Matrix V must be vector");
        }
        double norm = 0;
        for (int i = 0; i < V.getRows(); i++) {
            norm += Math.pow(V.getValue(i,0),2);
        }
        return Math.sqrt(norm);
    }
    public static Matrix interiorPointMethod(Matrix C, Matrix A, Matrix D, double alpha) {
        Matrix tildaA = A.multiply(D);
        //System.out.println("A~ :");
        //tildaA.printMatrix();
        Matrix tildaC = D.multiply(C);
        //System.out.println("C~ :");
        //tildaC.printMatrix();
        int N = A.getCols();
        int n = A.getRows();
        IdentityMatrix I = new IdentityMatrix(N);
        Matrix S1 = tildaA.multiply(tildaA.transpose()); // A~ * A~T
        //System.out.println("(A~ * A~T):");
        //S1.printMatrix();
        AugmentedMatrix AM = new AugmentedMatrix(n, S1);
        findInverseMatrix(AM, n);
        Matrix invS1 = AM.diagonalisation(); // (A~ * A~T)-1
        //System.out.println("(A~ * A~T)-1:");
        //invS1.printMatrix();
        Matrix S2 = tildaA.transpose().multiply(invS1); // A~T * (A~ * A~T)-1
        //System.out.println("A~T * (A~ * A~T)-1:");
        //S2.printMatrix();
        Matrix S3 = S2.multiply(tildaA); // A~T * (A~ * A~T)-1 * A~
        //System.out.println("A~T * (A~ * A~T)-1 * A~:");
        //S3.printMatrix();
        Matrix P = I.subtract(S3);
        //System.out.println("P:");
        //P.printMatrix();

        Matrix Cp = P.multiply(tildaC);
        //System.out.println("Cp:");
        //Cp.printMatrix();

        Matrix vector1s = new Matrix(N, 1);
        double v = 0;
        for (int i = 0; i < N; i++) {
            vector1s.setValue(i,0,1);
            v = Math.max(Math.abs(Cp.getValue(i,0)), v);
        }
        double alphaVcoeff = alpha / v; //v = 0?
        for (int i = 0; i < N; i++) {
            double newValue = Cp.getValue(i,0) * alphaVcoeff;
            Cp.setValue(i,0,newValue);
        } // a/v * Cp

        Matrix tildaX = vector1s.add(Cp);

        Matrix x = D.multiply(tildaX);
        return x;
    }

    public static double findDeterminant(int n, Matrix A) throws IllegalArgumentException {
        if (A.getRows() != A.getCols()) {
            throw new IllegalArgumentException("Unable to find determinant of a non-square matrix");
        }
        double determinant = 1;
        for (int i = 0; i < n-1; i++) {
            Matrix T = A.transpose();
            double maxPivot = T.getValue(i,i);
            int swapRow = i;
            for (int j = i; j < n; j++) {
                if (Math.abs(T.getValue(i,j)) > Math.abs(maxPivot)) {
                    swapRow = j;
                    maxPivot = T.getValue(i,j);
                }
            }
            if (swapRow > i) {
                PermutationMatrix P = new PermutationMatrix(n, i, swapRow);
                A = P.multiply(A);
                determinant *= -1;
            }
            for (int j = i+1; j < n; j++) {
                if (A.getValue(j, i) != 0) {
                    EliminationMatrix E = new EliminationMatrix(n, j, i, A);
                    A = E.multiply(A);
                }
            }
        }
        for (int i = 0; i < n; i++) {
            determinant *= A.getValue(i,i);
        }
        //System.out.println(determinant);
        return determinant;
    }

    public static Matrix findInverseMatrix(AugmentedMatrix AM, int n) throws IllegalArgumentException {
        //Matrix A = AM.L;
        //Matrix I = AM.R;
        int step = 0;
        if (findDeterminant(n, AM.L) == 0) {
            throw new IllegalArgumentException("Matrix A is singular");
        }
        for (int i = 0; i < n-1; i++) {
            Matrix T = AM.L.transpose();
            double maxPivot = T.getValue(i,i);
            int swapRow = i;
            for (int j = i; j < n; j++) {
                if (Math.abs(T.getValue(i,j)) > Math.abs(maxPivot)) {
                    swapRow = j;
                    maxPivot = T.getValue(i,j);
                }
            }
            if (swapRow > i) {
                PermutationMatrix P = new PermutationMatrix(n, i, swapRow);
                AM.L = P.multiply(AM.L);
                AM.R = P.multiply(AM.R);
                //System.out.println("step #" + ++step + ": permutation");
                //AM.display();

            }

            for (int j = i+1; j < n; j++) {
                if (AM.L.getValue(j, i) != 0) {
                    EliminationMatrix E = new EliminationMatrix(n, j, i, AM.L);
                    AM.L = E.multiply(AM.L);
                    AM.R = E.multiply(AM.R);
                    //System.out.println("step #" + ++step + ": elimination");
                    //AM.display();
                }
            }
        }
        for (int i = n-1; i > 0; i--) {
            for (int j = i - 1; j >= 0; j--) {
                if (AM.L.getValue(j, i) != 0) {
                    EliminationMatrix E = new EliminationMatrix(n, j, i, AM.L);
                    AM.L = E.multiply(AM.L);
                    AM.R = E.multiply(AM.R);
                    //System.out.println("step #" + ++step + ": elimination");
                    //AM.display();
                }
            }
        }
        return AM.R;
    }
}

class Matrix {
    protected double[][] matrix;
    protected int rows;
    protected int cols;

    public Matrix(int rows, int cols) {
        matrix = new double[rows][cols];
        this.rows = rows;
        this.cols = cols;
    }

    public void readMatrix(Scanner scanner) {
        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < cols; j++) {
                matrix[i][j] = scanner.nextDouble();
            }
        }
        scanner.nextLine();
    }
    public void printMatrix() {
        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < cols; j++) {
                System.out.printf("%.5f ", matrix[i][j]);
            }
            System.out.println();
        }
    }

    public Matrix add(Matrix other) throws IllegalArgumentException {
        if (this.rows != other.rows || this.cols != other.cols) {
            throw new IllegalArgumentException("Matrices have different dimensions");
        }
        Matrix result = new Matrix(this.rows, this.cols);
        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < cols; j++) {
                result.matrix[i][j] = matrix[i][j] + other.matrix[i][j];
            }
        }
        return result;
    }
    public Matrix subtract(Matrix other) throws IllegalArgumentException {
        if (this.rows != other.rows || this.cols != other.cols) {
            throw new IllegalArgumentException("Matrices have different dimensions");
        }
        Matrix result = new Matrix(this.rows, this.cols);
        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < cols; j++) {
                result.matrix[i][j] = matrix[i][j] - other.matrix[i][j];
            }
        }
        return result;
    }
    public Matrix multiply(Matrix other) throws IllegalArgumentException {
        if (this.cols != other.rows) {
            throw new IllegalArgumentException("Matrices cannot be multiplied");
        }
        Matrix result = new Matrix(this.rows, other.cols);
        for (int i = 0; i < rows; ++i) {
            for (int j = 0; j < other.cols; ++j) {
                for (int k = 0; k < cols; ++k) {
                    result.matrix[i][j] += matrix[i][k] * other.matrix[k][j];
                }
            }
        }
        return result;
    }
    public Matrix transpose() {
        Matrix result = new Matrix(cols, rows);
        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < cols; j++) {
                result.matrix[j][i] = matrix[i][j];
            }
        }
        return result;
    }
    public double getValue(int row, int col) {
        return matrix[row][col];
    }
    public void setValue(int row, int col, double value) {
        matrix[row][col] = value;
    }
    public int getRows() {
        return rows;
    }
    public int getCols() {
        return cols;
    }
}

class IdentityMatrix extends Matrix{
    public IdentityMatrix(int n) {
        super(n, n);
        for (int i = 0; i < n; i++) {
            matrix[i][i] = 1;
        }
    }
}

class EliminationMatrix extends IdentityMatrix {
    public EliminationMatrix(int n, int i, int j, Matrix A) {
        super(n);
        if (A.getValue(i, j) != 0) {
            matrix[i][j] = - A.getValue(i, j) / A.getValue(j, j);
        }
    }
}

class PermutationMatrix extends IdentityMatrix {
    public PermutationMatrix(int n, int i1, int i2) {
        super(n);
        double[] temp = matrix[i1];
        matrix[i1] = matrix[i2];
        matrix[i2] = temp;
    }
}

class DiagonalNormalizationMatrix extends IdentityMatrix {
    public DiagonalNormalizationMatrix(int n, double[] pivots) {
        super(n);
        for (int i = 0; i < n; i++) {
            matrix[i][i] = 1 / pivots[i];
        }
    }
}

class AugmentedMatrix {
    public Matrix L;
    public Matrix R;
    int N;

    AugmentedMatrix(int n, Matrix A) throws IllegalArgumentException {
        if (A.getCols() != A.getRows()) {
            throw new IllegalArgumentException("Augmented matrix can only be built for square matrices");
        }
        L = A;
        R = new IdentityMatrix(n);
        N = n;
    }

    Matrix diagonalisation() {
        Matrix A = L;
        Matrix I = R;
        double[] pivots = new double[N];
        for (int i = 0; i < N; i++) {
            pivots[i] = A.getValue(i, i);
        }
        DiagonalNormalizationMatrix D = new DiagonalNormalizationMatrix(N, pivots);
        A = D.multiply(A);
        I = D.multiply(I);
        return I;
    }

    void display() {
        for (int i = 0; i < N; i++) {
            for (int j = 0; j < N; j++) {
                System.out.printf("%.2f ", L.getValue(i, j));
            }
            for (int j = 0; j < N; j++) {
                System.out.printf("%.2f ", R.getValue(i, j));
            }
            System.out.println();
        }
        System.out.println();
    }
}