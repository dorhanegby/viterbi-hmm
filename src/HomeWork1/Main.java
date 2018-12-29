package HomeWork1;

public class Main {


    private static class Cell {
        int state;
        Cell parent;
        double value;

        public Cell(int state, Cell parent, double value) {
            this.state = state;
            this.parent = parent;
            this.value = value;
        }
    }

    public static double[][] transitions;
    public static double[][] emissions;

    public static int OFFSET = 1;

    public static final int S0=0;
    public static final int S1=1;
    public static final int S2=2;
    public static final int S3=3;
    public static final int S4=4;
    public static final int S5=5;
    public static final int S6=6;
    public static final int S7=7;
    public static final int S8=8;
    public static final int S9=9;

    public static final int K_STATES = 9;


    public static final int A = 0;
    public static final int G = 1;
    public static final int C = 2;
    public static final int T = 3;

    public static double p_1;
    public static double p_2;
    public static double p_3;
    public static double p_4;

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
        transitions = new double[10][10];
        transitions[S0][S1] = 1; // We assume the sequence starts in intergenic
        for(int j=2;j<K_STATES + 1;j++) {
            transitions[S0][j] = 0;
        }
        transitions[S1][S1] = 1 - p_3;
        transitions[S1][S2] = p_3;
        transitions[S2][S3] = 1.0;
        transitions[S3][S4] = (1 - p_4) * (1 - p_1);
        transitions[S3][S5] = (1 - p_4) * p_1;
        transitions[S3][S6] = p_4;
        transitions[S4][S7] = 1.0;
        transitions[S5][S8] = 1.0;
        transitions[S6][S9] = 1.0;
        transitions[S7][S4] = (1 - p_4) * (1 - p_1);
        transitions[S7][S5] = (1 - p_4) * (1 - p_1) * p_1;
        transitions[S7][S6] = p_4;
        transitions[S8][S4] = (1 - p_4) * (1 - p_1);
        transitions[S8][S5] = (1 - p_4) * (1 - p_1) * p_1;
        transitions[S8][S6] = p_4;
        transitions[S9][S1] = 1.0;
    }

    public static void initEmissions() {
        emissions = new double[10][4];
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

        emissions[S4][A] = p_1 /(1 - p_1);
        emissions[S4][G] = (0.5 - p_1) / (1 - p_1);
        emissions[S4][C] = (0.5 - p_1) / (1 - p_1);
        emissions[S4][T] = 0;

        emissions[S5][A] = 0;
        emissions[S5][G] = 0;
        emissions[S5][C] = 0;
        emissions[S5][T] = 1;

        emissions[S6][A] = 0;
        emissions[S6][G] = 0;
        emissions[S6][C] = 0;
        emissions[S6][T] = 1;

        emissions[S7][A] = p_1;
        emissions[S7][G] = 0.5 - p_1;
        emissions[S7][C] = 0.5 - p_1;
        emissions[S7][T] = p_1;

        emissions[S8][A] = 0;
        emissions[S8][G] = 0.5;
        emissions[S8][C] = 0.5;
        emissions[S8][T] = 0;

        emissions[S9][A] = 0;
        emissions[S9][G] = 1 - p_2;
        emissions[S9][C] = p_2;
        emissions[S9][T] = 0;
    }

    public static void initModel() {
        initTransitions();
        initEmissions();
    }

    private static void initViterbiMatrix(Cell[][] viterbiMatrix) {
        viterbiMatrix[0][0] = new Cell(0, null, 0); // Dummy State.

        for(int j=1;j<K_STATES + 1;j++) {
            viterbiMatrix[0][j] = new Cell(0, null, Double.NEGATIVE_INFINITY);
        }
    }

    private static Cell calculateMaxParent(int currentIndex, int state, Cell[][] viterbiMatrix, String sequence) {

        int x_i = baseToIndex(sequence.charAt(currentIndex - 1));

        double maxValue = -Double.POSITIVE_INFINITY;
        Cell maxParent = null;

        for (int j = 0; j < K_STATES + 1; j++) {
            double value = viterbiMatrix[currentIndex - 1][j].value + Math.log(transitions[j][state]);

           if(value > maxValue) {
               maxValue = value;
               maxParent = viterbiMatrix[currentIndex - 1][j];
           }
        }

        return new Cell(state, maxParent, maxValue + Math.log(emissions[state][x_i]));
    }

    private static void updateViterbiMatrix(Cell[][] viterbiMatrix, String sequence) {
        int n = sequence.length();

        for(int i = 1; i < n + 1;i++) {
            for(int j=0; j < K_STATES + 1 ; j++) {
                viterbiMatrix[i][j] = calculateMaxParent(i, j, viterbiMatrix, sequence);
            }
        }
    }


    private static void outputViterbiMatrix(Cell[][] viterbiMatrix, String sequence) {

        int n = sequence.length();

        Cell maxCell = viterbiMatrix[n][1];

        for(int j=2; j < K_STATES + 1; j++) {
            Cell current = viterbiMatrix[n][j];

            if(current.value > maxCell.value) {
                maxCell = current;
            }
        }

        String hmm = "";
        Cell current = maxCell;
        double maxLikelihood = maxCell.value;

        while(current != null) {
            hmm = current.state + hmm;
            current = current.parent;
        }
        System.out.println(hmm);
        System.out.println(" "+ sequence);

        System.out.println("log(P(X, S | HMM)) = " + maxLikelihood);

    }

    // V[i,j] = the probability of an annotation of
    // the prefix X1...Xi  that has the highest probability among
    // those that end with state sj
    public static void viterbi(String sequence) {
        int n = sequence.length(); // Size of the sequence.
        Cell[][] viterbiMatrix = new Cell[n + 1][K_STATES + 1];

        initViterbiMatrix(viterbiMatrix);
        updateViterbiMatrix(viterbiMatrix, sequence);
        outputViterbiMatrix(viterbiMatrix, sequence);
    }

    private static void initForwardMatrix(Cell[][] forwardMatrix) {
        initViterbiMatrix(forwardMatrix);
    }


    private static double findMax(double[] a) {
        double maxValue = a[0];
        for(int i=1;i<a.length;i++) {
            if(a[i] > maxValue) {
                maxValue = a[i];
            }
        }

        return maxValue;
    }

    private static Cell calculateSumOfLastColumn(int currentIndex, int state, Cell[][] forwardMatrix, String sequence) {

        int x_i = baseToIndex(sequence.charAt(currentIndex - 1));

        // Calculate sum of Logs - numeric stable.
        double [] a = new double[K_STATES + 1];
        double [] b = new double[K_STATES + 1];

        for (int l=0;l<a.length;l++) {
            a[l] = forwardMatrix[currentIndex - 1][l].value + Math.log(transitions[l][state]);
        }

        double maxA = findMax(a);
        double sum = 0;
        for (int i =0;i<b.length;i++) {

            if(Double.isInfinite(maxA)) {
                return new Cell(state, null, maxA);
            }

            b[i] = a[i] - maxA;
            sum += Math.exp(b[i]);

        }


        double sumOfLogs = maxA + Math.log(sum);


        return new Cell(state, null, sumOfLogs + Math.log(emissions[state][x_i]));

    }


    private static void updateForwardMatrix(Cell[][] forwardMatrix, String sequence) {

        int n = sequence.length();

        for(int i = 1; i < n + 1;i++) {
            for(int j=0; j < K_STATES + 1; j++) {
                forwardMatrix[i][j] = calculateSumOfLastColumn(i, j, forwardMatrix, sequence);
            }
        }

    }



    public static Cell[][] forward(String sequence) {
        int n = sequence.length(); // Size of the sequence.
        Cell[][] forwardMatrix = new Cell[n + 1][K_STATES + 1];

        initForwardMatrix(forwardMatrix);
        updateForwardMatrix(forwardMatrix, sequence);

        return forwardMatrix;
    }

    private static void initBackwardMatrix(Cell[][] backwardMatrix,String sequence) {
        int n = sequence.length(); // Size of the sequence.

        for(int j=0;j<K_STATES + 1;j++) {
            backwardMatrix[n][j] = new Cell(j, null, 0);
        }
    }

    private static Cell calculateSumOfFrontColumn(int currentIndex, int state, Cell[][] backardMatrix, String sequence) {

        int x_i_1 = baseToIndex(sequence.charAt(currentIndex));

        // Calculate sum of Logs - numeric stable.
        double [] a = new double[K_STATES + 1];
        double [] b = new double[K_STATES + 1];

        for (int l=0;l<a.length;l++) {
            a[l] = backardMatrix[currentIndex + 1][l].value + Math.log(transitions[state][l]) + Math.log(emissions[l][x_i_1]);
        }

        double maxA = findMax(a);
        double sum = 0;
        for (int i =0;i<b.length;i++) {

            if(Double.isInfinite(maxA)) {
                return new Cell(state, null,maxA);
            }

            b[i] = a[i] - maxA;
            sum += Math.exp(b[i]);
        }

        return new Cell(state, null,maxA + Math.log(sum));
    }

    private static void updateBackwardMatrix(Cell[][] backwardMatrix, String sequence) {
        int n = sequence.length(); // Size of the sequence.

        for(int i = n - 1; i >= 0; i--) {
            for(int j=0;j<K_STATES + 1; j++) {
                backwardMatrix[i][j] = calculateSumOfFrontColumn(i, j, backwardMatrix, sequence);
            }
        }
    }

    public static Cell[][] backward(String sequence) {
        int n = sequence.length(); // Size of the sequence.
        Cell[][] backwardMatrix = new Cell[n + 1][K_STATES + 1];

        initBackwardMatrix(backwardMatrix, sequence);
        updateBackwardMatrix(backwardMatrix, sequence);

        return backwardMatrix;
    }


    private static int findMaxInColumn(Cell[] col) {

        int argMax = 1;
        double maxValue = col[1].value;

        for(int j=2;j<col.length;j++) {
            if (col[j].value > maxValue) {
                argMax = j;
            }
        }

        return argMax;
    }

    private static Cell[] [] calculateMAP(String sequence) {
        int n = sequence.length(); // Size of the sequence.

        Cell[] [] forwardMatrix = forward(sequence);
        Cell[] [] backwardMatrix = backward(sequence);

        Cell[] [] map = new Cell[n + 1] [K_STATES + 1];

        for(int i =0;i<n+1;i++) {
            for(int j=0;j<K_STATES + 1;j++) {
                double value =  Math.log(Math.exp(forwardMatrix[i][j].value) * Math.exp(backwardMatrix[i][j].value));
                if(Double.isInfinite(value) || Double.isNaN(value)) {
                    map[i][j] = new Cell(j, null, Double.NEGATIVE_INFINITY);
                }
                else {
                    map[i][j] = new Cell(j, null, value);
                }

            }
        }

        String hmm = "";
        for(int i=1;i<n + 1;i++) {
            int argMax = findMaxInColumn(map[i]);
            hmm = hmm + argMax;
        }

        System.out.println(hmm);
        System.out.println(sequence);

        return map;

    }


    public static void main(String[] args) {
	    initModel();
        String seq = args[0];
        String method = args[1];
        p_1 = Double.parseDouble(args[2]);
        p_2 = Double.parseDouble(args[3]);
        p_3 = Double.parseDouble(args[4]);
        p_4 = Double.parseDouble(args[5]);
//        viterbi(SEQUENCE);

        calculateMAP(SEQUENCE);
    }

 /**  DEBUG **/


//    public static void printMatrix(Cell[][] matrix) {
//        int n = SEQUENCE.length(); // Size of the sequence.
//
//        String seq = "";
//        String hmm = "";
//        String prob = "";
//        for (int j = 0; j < K_STATES + 1; j++) {
//            for(int i =0;i<n+1;i++) {
//                if(Double.isInfinite(matrix[i][j].value)) {
//                    System.out.printf("-inf");
//                }
//                else {
//                    System.out.printf("%.2f", matrix[i][j].value);
//
//                }
//
//                System.out.print("\t");
//            }
//            System.out.println();
//        }
//    }
//
//    public static void printStringWithTab(String seq) {
//        for(int i=0;i<seq.length();i++) {
//            System.out.print(seq.charAt(i));
//            System.out.print("\t\t");
//        }
//        System.out.println();
//    }
//
//    public static void printProbForGene(Cell[][] map) {
//        int n = SEQUENCE.length(); // Size of the sequence.
//        double[] probs = new double[42];
//        for (int i = 1; i < n + 1; i++) {
//            double geneProb = 0;
//            double sumOfCol = sumOfCol(map[i]);
//            geneProb = (1 - ((Math.exp(map[i][1].value)) / sumOfCol));
//            probs[i] = geneProb;
//        }
//
//        for(int i=1;i<probs.length;i++) {
//            System.out.println("The probability that base " + i + " is inside a gene: " + probs[i]);
//        }
//    }
//
//    public static double sumOfCol(Cell[] arr) {
//        double sum = 0;
//        for(int i =0;i<arr.length;i++) {
//            if(Double.isInfinite(arr[i].value)) {
//                continue;
//            }
//            sum += Math.exp(arr[i].value);
//        }
//
//        return sum;
//    }

}

