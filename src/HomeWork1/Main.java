package HomeWork1;

public class Main {


    private class Cell {
        int state;
        Cell parent;
        double maxLikelihood;
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


    public static final int A = 0;
    public static final int G = 1;
    public static final int C = 2;
    public static final int T = 3;


    public static void main(String[] args) {
	    initModel();
    }

    public static void initModel() {
        initTransitions();
        initEmissions();
    }

    public static void initTransitions() {
        transitions = new double[8][8];
        transitions[0][0] = 0.9;
        transitions[0][1] = 0.1;
        transitions[1][2] = 1;
        transitions[2][3] = 0.8;
        transitions[2][4] = 0.2;
        transitions[3][5] = 1;
        transitions[4][6] = 0.5;
        transitions[4][7] = 0.5;
        transitions[5][3] = 0.8;
        transitions[2][4] = 0.2;
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
}
