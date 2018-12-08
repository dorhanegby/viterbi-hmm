package HomeWork1;

public class Main {


    private class Cell {
        int state;
        Cell parent;
        double maxLikelihood;
    }

    public static double[][] transitions;
    public double[][] emissions;


    public static void main(String[] args) {
	    initModel();
    }

    public static void initModel() {
        initTransitions();
        initEmissions();
    }

    public static void initTransitions() {
        transitions = new double[8][8];
        transitions[1][2] = 0.1;
    }

    public static void initEmissions() {
        emissions = new double[8][8];
    }
}
