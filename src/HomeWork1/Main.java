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
        System.out.println("hey");
    }
}
