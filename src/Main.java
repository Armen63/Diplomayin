import HavasaruneriHamakarger.DiferencialHavasarumneriMetod;

public class Main {
    public static void main(String[] args) {
//        HasarakIteracianeriExanak hasarakIteracianeriExanak =
//                new HasarakIteracianeriExanak(
//                        "2*x+log10(2*x+3)-1",
//                        "-0.5*log10(2*x+3)+0.5",
//                        0);
//        hasarakIteracianeriExanak.iterations(3);
//
//        NewtoniExanak newtoniExanak = new NewtoniExanak("x^3-2*x^2+1", 0.5);
//        newtoniExanak.iterations(3, 1);
//        newtoniExanak.iterations(3, 2);
//        newtoniExanak.iterations(3, 3);
//        newtoniExanak.iterations(3, 4);
//        newtoniExanak.iterations(3, 5);
//
//        HavasaruneriHamakarger.HasarakIteracianeriExanak hasarakIteracianeriExanak1 =
//                new HavasaruneriHamakarger.HasarakIteracianeriExanak(
//                        new String[]{"3*x0-cos(x1)-0.9", "sin(x0-0.6)-x1-1.6"},
//                        new String[]{"cos(x1)/3+0.3", "sin(x0-0.6)-1.6"},
//                        new double[]{-1, -1});
//        hasarakIteracianeriExanak1.iterations(3);

//        HavasaruneriHamakarger.NewtoniExanak NewtoniExanak2 =
//                new HavasaruneriHamakarger.NewtoniExanak(
//                        new String[]{"3*x0-cos(x1)-0.9", "sin(x0-0.6)-x1-1.6"},
//                        new double[]{2, 1});
//        NewtoniExanak2.iterations(3);

        DiferencialHavasarumneriMetod diferencialHavasarumneriMetod =
                new DiferencialHavasarumneriMetod(
                        new String[]{"3*x0-cos(x1)-0.9", "sin(x0-0.6)-x1-1.6"},
                        new double[]{2, 1},
                        new double[]{0.015, 1},
                        0.7,
                        1
                );
        diferencialHavasarumneriMetod.iterations(3);
    }
}
