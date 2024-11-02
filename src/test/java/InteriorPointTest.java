import org.junit.Test;

import static org.junit.Assert.assertEquals;

public class InteriorPointTest {

    @Test
    public void test1() {
        Matrix C = new Matrix(4, 1);
        C.matrix = new double[][]{
                {1},
                {1},
                {0},
                {0}
        };
        Matrix A = new Matrix(2, 4);
        A.matrix = new double[][]{
                {2, 4, 1, 0},
                {1, 3, 0, -1}
        };
        Matrix D = new Matrix(4, 4);
        D.matrix = new double[][]{
                {0.5, 0, 0, 0},
                {0, 3.5, 0, 0},
                {0, 0, 1, 0},
                {0, 0, 0, 2}
        };
        Matrix b = new Matrix(2, 1);
        b.matrix = new double[][]{
                {16},
                {9}
        };
        double eps = 0.0001;
        double alpha1 = 0.5;
        double alpha2 = 0.9;

        Matrix x = Main.computeIterations(C, A, D, alpha1, eps);
        assertEquals(x.getValue(0, 0), 5.99991, eps);
        assertEquals(x.getValue(1, 0), 1.00004, eps);
        assertEquals(x.getValue(2, 0), 0.00003, eps);
        assertEquals(x.getValue(3, 0), 0.00003, eps);
        x = Main.computeIterations(C, A, D, alpha2, eps);
        assertEquals(x.getValue(0, 0), 5.99998, eps);
        assertEquals(x.getValue(1, 0), 1.00001, eps);
        assertEquals(x.getValue(2, 0), 0.00000, eps);
        assertEquals(x.getValue(3, 0), 0.00001, eps);
    }

    @Test
    public void test2() {
        Matrix C = new Matrix(6, 1);
        C.matrix = new double[][]{
                {9},
                {10},
                {16},
                {0},
                {0},
                {0}
        };
        Matrix A = new Matrix(3, 6);
        A.matrix = new double[][]{
                {18, 15, 12, 1, 0, 0},
                {6, 4, 8, 0, 1, 0},
                {5, 3, 3, 0, 0, 1}
        };
        Matrix D = new Matrix(6, 6);
        D.matrix = new double[][]{
                {1, 0, 0, 0, 0, 0},
                {0, 1, 0, 0, 0, 0},
                {0, 0, 1, 0, 0, 0},
                {0, 0, 0, 315, 0, 0},
                {0, 0, 0, 0, 174, 0},
                {0, 0, 0, 0, 0, 169}
        };
        Matrix b = new Matrix(3, 1);
        b.matrix = new double[][]{
                {360},
                {192},
                {180}
        };
        double eps = 0.0001;
        double alpha1 = 0.5;
        double alpha2 = 0.9;

        Matrix x = Main.computeIterations(C, A, D, alpha1, eps);
        assertEquals(x.getValue(0, 0), 0.00002, eps);
        assertEquals(x.getValue(1, 0), 7.99993, eps);
        assertEquals(x.getValue(2, 0), 20.00001, eps);
        assertEquals(x.getValue(3, 0), 0.00053, eps);
        assertEquals(x.getValue(4, 0), 0.00007, eps);
        assertEquals(x.getValue(5, 0), 96.00007, eps);
        x = Main.computeIterations(C, A, D, alpha2, eps);
        assertEquals(x.getValue(0, 0), 0.00001, eps);
        assertEquals(x.getValue(1, 0), 7.99997, eps);
        assertEquals(x.getValue(2, 0), 20.00001, eps);
        assertEquals(x.getValue(3, 0), 0.00033, eps);
        assertEquals(x.getValue(4, 0), 0.00005, eps);
        assertEquals(x.getValue(5, 0), 96.00005, eps);
    }

    @Test
    public void test3() {
        Matrix C = new Matrix(6, 1);
        C.matrix = new double[][]{
                {-4},
                {1},
                {-2},
                {0},
                {0},
                {0}
        };
        Matrix A = new Matrix(3, 6);
        A.matrix = new double[][]{
                {1, 1, -1, 1, 0, 0},
                {2, -2, 1, 0, 1, 0},
                {6, 10, 0, 0, 0, 1}
        };
        Matrix D = new Matrix(6, 6);
        D.matrix = new double[][]{
                {1, 0, 0, 0, 0, 0},
                {0, 1, 0, 0, 0, 0},
                {0, 0, 1, 0, 0, 0},
                {0, 0, 0, 9, 0, 0},
                {0, 0, 0, 0, 29, 0},
                {0, 0, 0, 0, 0, 29}
        };
        Matrix b = new Matrix(3, 1);
        b.matrix = new double[][]{
                {10},
                {30},
                {45}
        };
        double eps = 0.0001;
        double alpha1 = 0.5;
        double alpha2 = 0.9;

        Matrix x = Main.computeIterations(C, A, D, alpha1, eps);
        assertEquals(x.getValue(0, 0), 0.00000, eps);
        assertEquals(x.getValue(1, 0), 4.49998, eps);
        assertEquals(x.getValue(2, 0), 0.00001, eps);
        assertEquals(x.getValue(3, 0), 5.50003, eps);
        assertEquals(x.getValue(4, 0), 38.99994, eps);
        assertEquals(x.getValue(5, 0), 0.00018, eps);
        x = Main.computeIterations(C, A, D, alpha2, eps);
        assertEquals(x.getValue(0, 0), 0.00000, eps);
        assertEquals(x.getValue(1, 0), 4.49999, eps);
        assertEquals(x.getValue(2, 0), 0.00000, eps);
        assertEquals(x.getValue(3, 0), 5.50003, eps);
        assertEquals(x.getValue(4, 0), 38.99999, eps);
        assertEquals(x.getValue(5, 0), 0.00005, eps);
    }

}
