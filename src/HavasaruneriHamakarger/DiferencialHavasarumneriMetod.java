package HavasaruneriHamakarger;

import org.la4j.Vector;
import org.la4j.vector.dense.BasicVector;
import org.mariuszgromada.math.mxparser.Argument;
import org.mariuszgromada.math.mxparser.Expression;

public class DiferencialHavasarumneriMetod {
    private Argument[] xArray;
    private Expression[] function;
    private Vector lyamna;
    private Vector kArray;
    private double h;
    private double delta;
    private int countOfIteration = 1;

    public DiferencialHavasarumneriMetod(String[] function, double[] xArray, double[] kArray, double h, double delta) {
        this.xArray = new Argument[xArray.length];
        this.function = new Expression[function.length];
        this.lyamna = new BasicVector(xArray.length);
        this.kArray = new BasicVector(kArray);
        this.h = h;
        this.delta = delta;

        for (int i = 0; i < xArray.length; i++) {
            this.xArray[i] = new Argument("x" + i, xArray[i]);
            this.function[i] = new Expression(function[i]);
        }

        for (int i = 0; i < xArray.length; i++) {
            for (int j = 0; j < xArray.length; j++) {
                this.function[i].addArguments(this.xArray[j]);
            }
        }
    }

    public void iterations(int precision) {
        System.out.println("Դիֆերենցյալ հավասարումների մեթոդ համակարգերի համար\n");

        for (int j = 0; j < xArray.length; j++) {
            System.out.print(function[j].getExpressionString() +
                    " = " + function[j].calculate() + ", ");
        }

        System.out.println();

        for (int j = 0; j < xArray.length; j++) {
            System.out.print("x" + j + " = " + xArray[j].getArgumentValue() + ", ");
        }

        System.out.println("");

        Vector functionValue = new BasicVector(xArray.length);
        Vector previousX = new BasicVector(xArray.length);
        Vector tmpX = new BasicVector(xArray.length);

        lyamna = getLyamna();

        for (int i = 0; i < xArray.length; i++) {
            previousX.set(i, (4 * xArray[i].getArgumentValue() - 2 * Math.pow(h, 2) * lyamna.get(i)) /
                    (2 + delta * h));
//            System.out.print(previousX.get(i) );
        }

        do {
            System.out.println("--------------------------------------" + countOfIteration++ + "------------------------------------------------");
            for (int i = 0; i < xArray.length; i++) {
                functionValue.set(i, function[i].calculate());
            }

            for (int i = 0; i < xArray.length; i++) {
                tmpX.set(i, (4 * xArray[i].getArgumentValue() - (2 - delta * h) * previousX.get(i) -
                        2 * Math.pow(h, 2) * lyamna.get(i)) / (2 + delta * h));
                previousX.set(i, xArray[i].getArgumentValue());
                xArray[i].setArgumentValue(tmpX.get(i));
            }

            lyamna = getLyamna();

            for (int j = 0; j < xArray.length; j++)
                System.out.print(function[j].getExpressionString() + " = " + function[j].calculate() + ", ");

            System.out.println();

            for (int j = 0; j < xArray.length; j++)
                System.out.print("x" + j + " = " + xArray[j].getArgumentValue() + ", ");
            System.out.println();
        }
        while (!end(function, precision));
    }

    private boolean end(Expression[] function, int precision) {
        for (int i = 0; i < function.length; i++) {
            if (function[i].calculate() <= -Math.pow(10, -precision) || function[i].calculate() >= Math.pow(10, -precision) || Math.round(function[i].calculate()) != 0) {
                System.out.println("                    function value  = " + "   " + function[i].calculate());
                return false;
            }
        }
        return true;
    }

    private Vector getLyamna() {
        double partOfLyamna = 0;
        Expression exp;
        Vector lyamna = new BasicVector(xArray.length);
        for (int i = 0; i < xArray.length; i++) {
            for (int j = 0; j < xArray.length; j++) {
                exp = new Expression("der(" + kArray.get(j) + "*(" + function[j].getExpressionString() + ")^2,x" + i + ")");
                for (int k = 0; k < xArray.length; k++) {
                    exp.addArguments(xArray[k]);
                }
                partOfLyamna += exp.calculate();
            }
            lyamna.set(i, partOfLyamna);
        }
        return lyamna;
    }
}