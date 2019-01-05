package HomeWork1;

public class Main {


    /** Dynamic Programming Matrix Cell **/

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

    /** HMM Model **/

    private static double[][] transitions;
    private static double[][] emissions;

    /** Inputs **/

    private static Method method;
    private static String SEQUENCE;

    /** Constants **/

    private static final int S0=0;
    private static final int S1=1;
    private static final int S2=2;
    private static final int S3=3;
    private static final int S4=4;
    private static final int S5=5;
    private static final int S6=6;
    private static final int S7=7;
    private static final int S8=8;
    private static final int S9=9;

    private static final int K_STATES = 9;

    private static final int A = 0;
    private static final int G = 1;
    private static final int C = 2;
    private static final int T = 3;

    private static final double EPSILON = 0.00001;

    /** Statistics **/

    private static double [][] N_trans = new double [10][10];
    private static double [][] N_emt = new double [10][4];

    /** Iterations **/
    private static double currentLikelihood = 0;
    private static double previousLikelihood = Double.NEGATIVE_INFINITY;
    private static boolean isFirstIteration = true;

    /** Helpers **/

    private static int baseToIndex(char base) {
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

    private static double findMax(double[] a) {
        double maxValue = a[0];
        for(int i=1;i<a.length;i++) {
            if(a[i] > maxValue) {
                maxValue = a[i];
            }
        }

        return maxValue;
    }

    private static double sumOfExponents(double [] a, double [] b, int state, double maxA) {
        double sum = 0;

        for (int i =0;i<b.length;i++) {

            b[i] = a[i] - maxA;
            sum += Math.exp(b[i]);

        }

        return sum;
    }

    /** Model Helpers **/

    private static void initTransitions(double p_1, double p_2, double p_3, double p_4) {
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
        transitions[S7][S5] = (1 - p_4) * p_1;
        transitions[S7][S6] = p_4;
        transitions[S8][S4] = (1 - p_4) * (1 - p_1);
        transitions[S8][S5] = (1 - p_4) * p_1;
        transitions[S8][S6] = p_4;
        transitions[S9][S1] = 1.0;
    }

    private static void initEmissions(double p_1, double p_2, double p_3, double p_4) {
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

        emissions[S4][A] = p_1 / (1 - p_1);
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

        emissions[S8][A] = 0.5;
        emissions[S8][G] = 0;
        emissions[S8][C] = 0;
        emissions[S8][T] = 0.5;

        emissions[S9][A] = 0;
        emissions[S9][G] = 1 - p_2;
        emissions[S9][C] = p_2;
        emissions[S9][T] = 0;
    }

    private static void initModel(double p_1, double p_2, double p_3, double p_4) {
        System.out.println("p_1: " + p_1);
        System.out.println("p_2: " + p_2);
        System.out.println("p_3: " + p_3);
        System.out.println("p_4: " + p_4);
        initTransitions(p_1, p_2, p_3, p_4);
        initEmissions(p_1, p_2, p_3, p_4);
    }

    /** Viterbi Helpers **/

    private static Cell calculateMaxParent(int currentIndex, int state, Cell[][] viterbiMatrix, String sequence) {

        int x_i = baseToIndex(sequence.charAt(currentIndex - 1));

        double maxValue = Double.NEGATIVE_INFINITY;
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

    /** Forward Helpers **/

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

        if(Double.isInfinite(maxA)) {
            return new Cell(state, null, maxA);
        }

        sum = sumOfExponents(a, b, state, maxA);


        double sumOfLogs = maxA + Math.log(sum);


        return new Cell(state, null, sumOfLogs + Math.log(emissions[state][x_i]));

    }

    /** Backward Helpers **/

    private static Cell calculateSumOfFrontColumn(int currentIndex, int state, Cell[][] backwardMatrix, String sequence) {

        int x_i_1 = baseToIndex(sequence.charAt(currentIndex));

        // Calculate sum of Logs - numeric stable.
        double [] a = new double[K_STATES + 1];
        double [] b = new double[K_STATES + 1];

        for (int l=0;l<a.length;l++) {
            a[l] = backwardMatrix[currentIndex + 1][l].value + Math.log(transitions[state][l]) + Math.log(emissions[l][x_i_1]);
        }

        double maxA = findMax(a);
        if(Double.isInfinite(maxA)) {
            return new Cell(state, null,maxA);
        }

        double sum = sumOfExponents(a, b, state, maxA);

        return new Cell(state, null,maxA + Math.log(sum));
    }

    /** MAP Helpers **/

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

    /** Viterbi **/

    private static void initViterbiMatrix(Cell[][] viterbiMatrix) {
        viterbiMatrix[0][0] = new Cell(0, null, 0); // Dummy State.

        for(int j=1;j<K_STATES + 1;j++) {
            viterbiMatrix[0][j] = new Cell(0, null, Double.NEGATIVE_INFINITY);
        }
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

        N_trans = new double [10][10];
        N_emt = new double [10][4];


        int wordPointer = sequence.length() - 1;
        while(current != null) {
            if(current.parent != null) {
                N_trans[current.parent.state][current.state]++;
            }
            if(wordPointer >= 0) {
                int base = baseToIndex(sequence.charAt(wordPointer));
                N_emt[current.state][base]++;
            }
            hmm = current.state + hmm;
            current = current.parent;
            wordPointer--;
        }
        System.out.println(hmm);
        System.out.println(" "+ sequence);

        System.out.println("log(P(X, S | HMM)) = " + maxLikelihood);
        if(!isFirstIteration) {
            previousLikelihood = currentLikelihood;
        }
        isFirstIteration = false;
        currentLikelihood = maxLikelihood;

    }

    /** Baum-Welch Helpers **/


    private static void updatedExpectedEmission(Cell[][] forward, Cell[][] backward) {
        double likelihood = calculateLikelihoodFromForward(forward);

        for(int j=1;j<K_STATES + 1; j++) {
            for(int sigma=0;sigma<4;sigma++) {
                double sumOfPostEmissions = calculatePotEmissions(forward, backward, sigma, j);
                N_emt[j][sigma] = sumOfPostEmissions;
            }
        }

        System.out.println("log(P(X, S | HMM)) = " + likelihood);
        if(!isFirstIteration) {
            previousLikelihood = currentLikelihood;
        }
        isFirstIteration = false;
        currentLikelihood = likelihood;
    }

    private static double calculatePotEmissions(Cell[][] forward, Cell[][] backward, int sigma, int state) {
        int n = SEQUENCE.length();
        double sum = 0;
        for(int i=1;i<n + 1;i++) {
            if(baseToIndex(SEQUENCE.charAt(i - 1)) != sigma) {
                continue;
            }
            sum += Math.exp(forward[i][state].value) * Math.exp(backward[i][state].value);
        }

        return sum;
    }


    private static void updateExpectedTransition(Cell[][] forward, Cell[][] backward) {
        for(int j=1;j<K_STATES + 1; j++) {
            for(int l=1;l<K_STATES + 1;l++) {
                if(j == 1 && l == 2) {
                    System.out.println();
                }
                double sumOfPaths = calculateSumOfPaths(forward, backward, j, l);

                N_trans[j][l] = transitions[j][l] * sumOfPaths;
            }
        }

    }

    private static double calculateSumOfPaths(Cell[][] forward, Cell[][] backward, int startState, int endState) {
        int n = SEQUENCE.length();
        double sum = 0;
        for(int i=1;i<n;i++) {
            int x_i = baseToIndex(SEQUENCE.charAt(i));
            sum += Math.exp(forward[i][startState].value) * emissions[endState][x_i] * Math.exp(backward[i + 1][endState].value);
        }

        return sum;
    }

    private static double calculateLikelihoodFromForward(Cell[][] forward) {
        double sum = 0;
        Cell[] lastColumn = forward[SEQUENCE.length()];

        for(int i=0;i<lastColumn.length;i++) {
            if(Double.isInfinite(lastColumn[i].value)) {
                continue;
            }
            sum += lastColumn[i].value;
        }

        return sum;
    }

    // V[i,j] = the probability of an annotation of
    // the prefix X1...Xi  that has the highest probability among
    // those that end with state sj
    private static void viterbi(String sequence) {
        int n = sequence.length(); // Size of the sequence.
        Cell[][] viterbiMatrix = new Cell[n + 1][K_STATES + 1];

        initViterbiMatrix(viterbiMatrix);
        updateViterbiMatrix(viterbiMatrix, sequence);
        outputViterbiMatrix(viterbiMatrix, sequence);
    }

    /** Forward **/
    private static void initForwardMatrix(Cell[][] forwardMatrix) {
        initViterbiMatrix(forwardMatrix);
    }

    private static void updateForwardMatrix(Cell[][] forwardMatrix, String sequence) {

        int n = sequence.length();

        for(int i = 1; i < n + 1;i++) {
            for(int j=0; j < K_STATES + 1; j++) {
                forwardMatrix[i][j] = calculateSumOfLastColumn(i, j, forwardMatrix, sequence);
            }
        }

    }

    private static Cell[][] forward(String sequence) {
        int n = sequence.length(); // Size of the sequence.
        Cell[][] forwardMatrix = new Cell[n + 1][K_STATES + 1];

        initForwardMatrix(forwardMatrix);
        updateForwardMatrix(forwardMatrix, sequence);

        return forwardMatrix;
    }

    /** Backward **/

    private static void initBackwardMatrix(Cell[][] backwardMatrix,String sequence) {
        int n = sequence.length(); // Size of the sequence.

        for(int j=0;j<K_STATES + 1;j++) {
            backwardMatrix[n][j] = new Cell(j, null, 0);
        }
    }

    private static void updateBackwardMatrix(Cell[][] backwardMatrix, String sequence) {
        int n = sequence.length(); // Size of the sequence.

        for(int i = n - 1; i >= 0; i--) {
            for(int j=0;j<K_STATES + 1; j++) {
                backwardMatrix[i][j] = calculateSumOfFrontColumn(i, j, backwardMatrix, sequence);
            }
        }
    }

    private static Cell[][] backward(String sequence) {
        int n = sequence.length(); // Size of the sequence.
        Cell[][] backwardMatrix = new Cell[n + 1][K_STATES + 1];

        initBackwardMatrix(backwardMatrix, sequence);
        updateBackwardMatrix(backwardMatrix, sequence);

        return backwardMatrix;
    }

    /** MAP **/

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

    /** Parameter inference **/

    private static void runParamInfer() {
        if(method == Method.Viterbi) {
            runViterbiTraining();
        }

        else {
            runBaumWelch();
        }

    }

    private static void updateParameters() {

        double p_1;
        double numerator = (N_emt[7][A] + N_emt[7][T] + N_trans[3][5] + N_trans[7][5] + N_trans[8][5] + N_emt[4][A]);

        if(numerator == 0) {
            p_1 = 0;
        }
        else {
            p_1 =  (1.0 / 2.0) * (numerator) / (numerator + N_emt[7][C] + N_emt[7][G] + N_emt[4][C] + N_emt[4][G]);
        }

        double p_2;

        if(N_emt[9][C] == 0) {
            p_2 = 0;
        }
        else {
            p_2 = (N_emt[9][C]) / (N_emt[9][C] + N_emt[9][G]);
        }

        double p_3;

        if(N_trans[1][2] == 0) {
            p_3 = 0;
        }
        else {
            p_3 = (N_trans[1][2]) / (N_trans[1][1] + N_trans[1][2]);
        }

        double p_4;
        numerator = N_trans[3][6] + N_trans[7][6] + N_trans[8][6];

        if(numerator == 0) {
            p_4 = 0;
        }
        else {
            p_4 = numerator / (numerator + N_trans[3][4] + N_trans[3][5] + N_trans[7][4] + N_trans[7][5] + N_trans[8][4] + N_trans[8][5]);
        }


        initModel(p_1, p_2, p_3, p_4);
    }

    /** Viterbi training **/

    private static void runViterbiTraining() {
        while(currentLikelihood - previousLikelihood > EPSILON) {
            viterbi(SEQUENCE);
            updateParameters();
        }
    }

    /** Baum-Welch **/

    private static void runBaumWelch() {
        while(currentLikelihood - previousLikelihood > EPSILON) {
            Cell[][] forward = forward(SEQUENCE);
            Cell[][] backward = backward(SEQUENCE);
            updateExpectedTransition(forward, backward);
            updatedExpectedEmission(forward, backward);
            updateParameters();
        }
    }

    public static void main(String[] args) {

        SEQUENCE = args[0];
        method = args[1].equals("V") ? Method.Viterbi : Method.BW;
        double p_1 = Double.parseDouble(args[2]);
        double p_2 = Double.parseDouble(args[3]);
        double p_3 = Double.parseDouble(args[4]);
        double p_4 = Double.parseDouble(args[5]);

        initModel(p_1, p_2, p_3, p_4);

        runParamInfer();
    }


    enum Method {
        Viterbi,
        BW
    }

}

