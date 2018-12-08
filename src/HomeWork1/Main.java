package HomeWork1;

public class Main {


    private static class Cell {
        int state;
        Cell parent;
        double maxLikelihood;

        public Cell(int state, Cell parent, double maxLikelihood) {
            this.state = state;
            this.parent = parent;
            this.maxLikelihood = maxLikelihood;
        }
    }

    public static double[][] transitions;
    public static double[][] emissions;

    public static final int S1=0;
    public static final int S2=1;
    public static final int S3=2;
    public static final int S4=3;
    public static final int S5=4;
    public static final int S6=5;
    public static final int S7=6;
    public static final int S8=7;

    public static final int K_STATES = 8;


    public static final int A = 0;
    public static final int G = 1;
    public static final int C = 2;
    public static final int T = 3;

    public static final String SEQUENCE = "CCATCGCACTAGGGACGGTGGTCCGACGCACATGTTGCTCC";

    public static int baseToIndex(char base) {
        if(base == 'A') {
            return A;
        }
        else if(base == 'G') {
            return G;
        }
        else if (base == 'C') {
            return C;
        }

        return T;
    }

    public static void initTransitions() {
        transitions = new double[8][8];
        transitions[S1][S1] = 0.9;
        transitions[S1][S2] = 0.1;
        transitions[S2][S3] = 1.0;
        transitions[S3][S4] = 0.8;
        transitions[S3][S5] = 0.2;
        transitions[S4][S6] = 1.0;
        transitions[S5][S7] = 0.5;
        transitions[S5][S8] = 0.5;
        transitions[S6][S4] = 0.8;
        transitions[S6][S5] = 0.2;
        transitions[S7][S4] = 0.8;
        transitions[S7][S5] = 0.2;
        transitions[S8][S1] = 1.0;
    }

    public static void initEmissions() {
        emissions = new double[8][4];
        emissions[S1][A] = 0.3;
        emissions[S1][G] = 0.2;
        emissions[S1][C] = 0.2;
        emissions[S1][T] = 0.3;

        emissions[S2][A] = 1;
        emissions[S2][G] = 0;
        emissions[S2][C] = 0;
        emissions[S2][T] = 0;

        emissions[S3][A] = 0;
        emissions[S3][G] = 0;
        emissions[S3][C] = 1;
        emissions[S3][T] = 0;

        emissions[S4][A] = 0.25;
        emissions[S4][G] = 0.375;
        emissions[S4][C] = 0.375;
        emissions[S4][T] = 0;

        emissions[S5][A] = 0;
        emissions[S5][G] = 0;
        emissions[S5][C] = 0;
        emissions[S5][T] = 1;

        emissions[S6][A] = 0.2;
        emissions[S6][G] = 0.3;
        emissions[S6][C] = 0.3;
        emissions[S6][T] = 0.2;

        emissions[S7][A] = 0.5;
        emissions[S7][G] = 0;
        emissions[S7][C] = 0;
        emissions[S7][T] = 0.5;

        emissions[S8][A] = 0;
        emissions[S8][G] = 0.5;
        emissions[S8][C] = 0.5;
        emissions[S8][T] = 0;
    }

    public static void initModel() {
        initTransitions();
        initEmissions();
    }

    private static void initViterbiMatrix(Cell[][] viterbiMatrix) {
        viterbiMatrix[0][0] = new Cell(0, null, 1);

        for(int j=1;j<K_STATES;j++) {
            viterbiMatrix[0][j] = new Cell(0, null, 0);
        }
    }

    private static Cell calculateMaxParent(int currentIndex, int state, Cell[][] viterbiMatrix, String sequence) {

        // TODO: add log probabilities

        int x_i = baseToIndex(sequence.charAt(currentIndex));

        double maxValue = -Double.MAX_VALUE;
        Cell maxParent = null;

        for (int j = 0; j < K_STATES; j++) {

           double value = viterbiMatrix[currentIndex - 1][j].maxLikelihood * transitions[j][state] * emissions[state][x_i];
           if(value > maxValue) {
               maxValue = value;
               maxParent = viterbiMatrix[currentIndex - 1][j];
           }
        }

        return new Cell(state, maxParent, maxValue);
    }

    private static void updateViterbiMatrix(Cell[][] viterbiMatrix, String sequence) {
        int n = sequence.length();

        for(int i = 1; i < n;i++) {
            for(int j=0; j < K_STATES ; j++) {
                viterbiMatrix[i][j] = calculateMaxParent(i, j, viterbiMatrix, sequence);
            }
        }
        System.out.println();
    }

    // V[i,j] = the probability of an annotation of
    // the prefix X1...Xi  that has the highest probability among
    // those that end with state sj
    public static void viterbi(String sequence) {
        int n = sequence.length(); // Size of the sequence.
        Cell[][] viterbiMatrix = new Cell[n][K_STATES];

        initViterbiMatrix(viterbiMatrix);
        updateViterbiMatrix(viterbiMatrix, sequence);

    }



    public static void main(String[] args) {
	    initModel();

        viterbi(SEQUENCE);
    }




}
