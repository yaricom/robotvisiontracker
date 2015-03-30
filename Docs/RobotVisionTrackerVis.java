import javax.imageio.ImageIO;
import javax.swing.*;
import java.awt.*;
import java.awt.image.BufferedImage;
import java.awt.image.DataBufferByte;
import java.awt.event.KeyAdapter;
import java.awt.event.KeyEvent;
import java.awt.event.WindowAdapter;
import java.awt.event.WindowEvent;
import java.io.*;
import java.util.Arrays;
import java.util.*;
import java.io.BufferedReader;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.PrintWriter;
import java.security.SecureRandom;

public class RobotVisionTrackerVis {

    public static int NUM_OF_TESTING_FRAMES = 50;
    public static long seed = 1;
    public static int delay = 100;
    public static boolean debug = true;
    public static boolean visualize = true;
    public static String execCommand = null;
    public static String testingFile = null;
    public static String trainingFile = null;
    public static String folder = "";
    public SecureRandom rnd;
    public static int W,H;
    public static Process solution;
    public final Object worldLock = new Object();

    public void printMessage(String s) {
        if (debug) {
            System.out.println(s);
        }
    }

    class Drawer extends JFrame {

        public DrawerPanel panel;
        public int width, height;
        public boolean pauseMode = false;

        class DrawerKeyListener extends KeyAdapter {
            public void keyPressed(KeyEvent e) {
                synchronized (keyMutex) {
                    if (e.getKeyChar() == ' ') {
                        pauseMode = !pauseMode;
                    }
                    keyPressed = true;
                    keyMutex.notifyAll();
                }
            }
        }

        public BufferedImage imgLeft = null;
        public BufferedImage imgRight = null;

        public int gtfLeftX, gtfLeftY, gtfRightX, gtfRightY;
        public int userLeftX, userLeftY, userRightX, userRightY;

        public void update(String leftFile, String rightFile, int gtfLeftX, int gtfLeftY, int gtfRightX, int gtfRightY,
                           int userLeftX, int userLeftY, int userRightX, int userRightY) throws Exception {
           synchronized (worldLock) {
             imgLeft = ImageIO.read(new File(leftFile));
             imgRight = ImageIO.read(new File(rightFile));
             this.gtfLeftX = gtfLeftX;
             this.gtfLeftY = gtfLeftY;
             this.gtfRightX = gtfRightX;
             this.gtfRightY = gtfRightY;
             this.userLeftX = userLeftX;
             this.userLeftY = userLeftY;
             this.userRightX = userRightX;
             this.userRightY = userRightY;
           }

        }

        class DrawerPanel extends JPanel {

            public void paint(Graphics g) {
                synchronized (worldLock) {

                    if (imgLeft!=null)
                    {
                        g.drawImage(imgLeft,0,0,null);
                        if (gtfLeftX>=0) {
                            g.setColor(Color.RED);
                            g.fillOval(gtfLeftX-5, gtfLeftY-5, 10, 10);
                        }
                        if (userLeftX>=0) {
                            g.setColor(Color.GREEN);
                            g.fillOval(userLeftX-5, userLeftY-5, 10, 10);
                        }
                    }

                    if (imgRight!=null)
                    {
                        g.drawImage(imgRight,640,0,null);
                        if (gtfRightX>=0) {
                            g.setColor(Color.RED);
                            g.fillOval(640+gtfRightX-5, gtfRightY-5, 10, 10);
                        }
                        if (userRightX>=0) {
                            g.setColor(Color.GREEN);
                            g.fillOval(640+userRightX-5, userRightY-5, 10, 10);
                        }
                    }
                }
            }
        }

        class DrawerWindowListener extends WindowAdapter {
            public void windowClosing(WindowEvent event) {
                RobotVisionTrackerVis.stopSolution();
                System.exit(0);
            }
        }

        final Object keyMutex = new Object();
        boolean keyPressed;

        public void processPause() {
            synchronized (keyMutex) {
                if (!pauseMode) {
                    return;
                }
                keyPressed = false;
                while (!keyPressed) {
                    try {
                        keyMutex.wait();
                    } catch (InterruptedException e) {
                        // do nothing
                    }
                }
            }
        }

        public Drawer() {
            super();

            panel = new DrawerPanel();
            getContentPane().add(panel);

            addWindowListener(new DrawerWindowListener());

            width = 640*2;
            height = 480;

            addKeyListener(new DrawerKeyListener());

            setSize(width, height);
            setTitle("Handle Tracking Marathon Match");

            setResizable(false);
            setVisible(true);
        }
    }


    public int[] imageToArray(String imageFile) throws Exception {
        printMessage("Reading image from " + imageFile);
        BufferedImage img = ImageIO.read(new File(imageFile));

        H = img.getHeight();
        W = img.getWidth();
        int[] res = new int[H * W];

        int pos = 0;
        byte[] pixels = ((DataBufferByte) img.getRaster().getDataBuffer()).getData();
        for (int i = 0; i < pixels.length; i+=3) {
            int v0 = (int)(pixels[i]);
            if (v0<0) v0 += 256;
            int v1 = (int)(pixels[i+1]);
            if (v1<0) v1 += 256;
            int v2 = (int)(pixels[i+2]);
            if (v2<0) v2 += 256;
            res[pos++] = v0 | (v1<<8) | (v2<<16);
        }
        return res;
    }

    class VideoFrame {
        int[] leftImg;
        int[] rightImg;

        public void Load(String fileNameLeft, String fileNameRight) throws Exception {
            this.leftImg = imageToArray(fileNameLeft);
            this.rightImg = imageToArray(fileNameRight);
        }
    }

    class VideoAnnotations {
        ArrayList<Integer> leftX = new ArrayList<Integer>();
        ArrayList<Integer> rightX = new ArrayList<Integer>();
        ArrayList<Integer> leftY = new ArrayList<Integer>();
        ArrayList<Integer> rightY = new ArrayList<Integer>();
        int numFrames;

        public void Load(String fileName, boolean bLeft) throws Exception {

            if (bLeft) {
                leftX.clear();
                leftY.clear();
                rightX.clear();
                rightY.clear();
            }

            BufferedReader br = new BufferedReader(new FileReader(fileName));
            printMessage("Loading annotation file " + fileName);

            // read annotations
            for (;;)
            {
                String sname = br.readLine();
                if (sname==null) break;
                int num = Integer.parseInt(sname.substring(0, 5));
                while (leftX.size()<=num) {
                    leftX.add(-1);
                    leftY.add(-1);
                    rightX.add(-1);
                    rightY.add(-1);
                }
                String scnt = br.readLine();
                String spnt = br.readLine();
                String[] items = spnt.split(" ");
                if (bLeft)
                {
                    leftX.set(num, (int)Double.parseDouble(items[0]));
                    leftY.set(num, (int)Double.parseDouble(items[1]));
                } else
                {
                    rightX.set(num, (int)Double.parseDouble(items[0]));
                    rightY.set(num, (int)Double.parseDouble(items[1]));
                }
            }
            numFrames = leftX.size();
        }
    }

    public void correctOutOfBounds(int[] userX, int[] userY) {
        for (int i=0;i<userX.length;i++) {
            if (userX[i]<0 || userX[i]>639 || userY[i]<0 || userY[i]>479) {
                userX[i] = -1;
                userY[i] = -1;
            }
        }
    }

    public double countFrameCorrect(int[] gtfX, int[] gtfY, int[] userX, int[] userY, int radi) {
        double p = 0;
        for (int i=0;i<gtfX.length;i++) {
            if (gtfX[i]<0 && userX[i]<0) {
                // correctly detected that object is not in view
                p += 1.0 / NUM_OF_TESTING_FRAMES;
                continue;
            }
            if (gtfX[i]<0 && userX[i]>=0) {
                // incorrectly detected an object
                continue;
            }
            if (gtfX[i]>=0 && userX[i]<0) {
                // incorrectly detected no object
                continue;
            }
            double dx = (gtfX[i] - userX[i]);
            double dy = (gtfY[i] - userY[i]);
            if (dx*dx+dy*dy <= radi) p += 1.0 / NUM_OF_TESTING_FRAMES;
        }
        return p;
    }

    public double doExec() throws Exception {

        try {
            rnd = SecureRandom.getInstance("SHA1PRNG");
        } catch (Exception e) {
            System.err.println("ERROR: unable to generate test case.");
            System.exit(1);
        }
        rnd.setSeed(seed);        
        
        // launch solution
        printMessage("Executing your solution: " + execCommand);
        solution = Runtime.getRuntime().exec(execCommand);

        BufferedReader reader = new BufferedReader(new InputStreamReader(solution.getInputStream()));
        PrintWriter writer = new PrintWriter(solution.getOutputStream());
        new ErrorStreamRedirector(solution.getErrorStream()).start();

        VideoFrame aFrame = new VideoFrame();
        VideoAnnotations aAnnotations = new VideoAnnotations();
        if (trainingFile != null) {
            BufferedReader br = new BufferedReader(new FileReader(trainingFile));
            int NumTrainingVideos = Integer.parseInt(br.readLine());
            printMessage("Training with " + NumTrainingVideos + " videos");

            writer.println(NumTrainingVideos);
            writer.flush();
            for (int i=0;i<NumTrainingVideos;i++)
            {
                String sVideoName = br.readLine();
                // load annotations
                aAnnotations.Load(folder+sVideoName+"_Left_annotation_pt.txt", true);
                aAnnotations.Load(folder+sVideoName+"_Right_annotation_pt.txt", false);

                writer.println(aAnnotations.numFrames-1);
                writer.flush();
                boolean bExitTraining = false;
                for (int f=1;f<aAnnotations.numFrames;f++)
                {
                    String leftFilename = folder+sVideoName+"_Left_"+String.format("%05d.png",f);
                    String rightFilename = folder+sVideoName+"_Right_"+String.format("%05d.png",f);

                    aFrame.Load(leftFilename, rightFilename);
                    for (int v : aFrame.leftImg) {
                        writer.println(v);
                    }
                    writer.flush();
                    for (int v : aFrame.rightImg) {
                        writer.println(v);
                    }
                    writer.flush();
                    int leftX = aAnnotations.leftX.get(f);
                    int leftY = aAnnotations.leftY.get(f);
                    int rightX = aAnnotations.rightX.get(f);
                    int rightY = aAnnotations.rightY.get(f);
                    writer.println(leftX);
                    writer.println(leftY);
                    writer.println(rightX);
                    writer.println(rightY);
                    writer.flush();
                    // call training function
                    // training(i+1, f, aFrame.leftImg, aFrame.rightImg, leftX, leftY, rightX, rightY);
                    int ret = Integer.parseInt(reader.readLine());
                    if (ret==1)
                    {
                        bExitTraining = true; // stop receiving training images
                        break;
                    }
                }
                if (bExitTraining) break;
            }
            br.close();
        } else
        {
            System.out.println("ERROR: Training file not provided");
            System.exit(0);
        }

        // call doneTesting function

        double score = 0;

        int[] gtfLeftX = new int[NUM_OF_TESTING_FRAMES];
        int[] gtfLeftY = new int[NUM_OF_TESTING_FRAMES];
        int[] gtfRightX = new int[NUM_OF_TESTING_FRAMES];
        int[] gtfRightY = new int[NUM_OF_TESTING_FRAMES];
        int[] userLeftX = new int[NUM_OF_TESTING_FRAMES];
        int[] userLeftY = new int[NUM_OF_TESTING_FRAMES];
        int[] userRightX = new int[NUM_OF_TESTING_FRAMES];
        int[] userRightY = new int[NUM_OF_TESTING_FRAMES];

        Drawer drawer = null;
        if (visualize) {
            drawer = new Drawer();
            drawer.pauseMode = false;
        }

        VideoAnnotations[] aTestAnnotations = new VideoAnnotations[5];
        for (int i=0;i<5;i++) aTestAnnotations[i] = new VideoAnnotations();

        if (testingFile != null) {
            BufferedReader br = new BufferedReader(new FileReader(testingFile));
            int N = Integer.parseInt(br.readLine());
            // pick random video
            ArrayList<String> testVideos = new ArrayList<String>();
            for (int i=0;i<N;i++) {
                testVideos.add(br.readLine());
            }
            br.close();

            String[] sTestVideoName = new String[5];
            int[] frameStart = new int[5];
            for (int i=0;i<5;i++) {
                int iN = rnd.nextInt(N);
                sTestVideoName[i] = testVideos.get(iN);
                // load annotations for video
                aTestAnnotations[i].Load(folder+sTestVideoName[i]+"_Left_annotation_pt.txt", true);
                aTestAnnotations[i].Load(folder+sTestVideoName[i]+"_Right_annotation_pt.txt", false);
                // pick random starting frame
                frameStart[i] = Math.max(1, 1+rnd.nextInt(aTestAnnotations[i].numFrames-NUM_OF_TESTING_FRAMES/5));
                printMessage("Testing with video " + sTestVideoName[i] + " starting at frame " + frameStart[i]);
            }


            long processingTime = 0;

            writer.println(NUM_OF_TESTING_FRAMES);
            writer.flush();
            for (int f=0;f<NUM_OF_TESTING_FRAMES;f++)
            {
                userLeftX[f] = userLeftY[f] = userRightX[f] = userRightY[f] = -1;
            }
            for (int f=0;f<NUM_OF_TESTING_FRAMES;f++)
            {
                int videoIndex = f/10;
                int frameIndex = f%10;
                writer.println(videoIndex); 
                writer.println(frameIndex);
                writer.flush();
                String leftFilename = folder+sTestVideoName[videoIndex]+"_Left_"+String.format("%05d.png",frameIndex+frameStart[videoIndex]);
                String rightFilename = folder+sTestVideoName[videoIndex]+"_Right_"+String.format("%05d.png",frameIndex+frameStart[videoIndex]);
                aFrame.Load(leftFilename, rightFilename);
                gtfLeftX[f] = aTestAnnotations[videoIndex].leftX.get(frameIndex+frameStart[videoIndex]);
                gtfLeftY[f] = aTestAnnotations[videoIndex].leftY.get(frameIndex+frameStart[videoIndex]);
                gtfRightX[f] = aTestAnnotations[videoIndex].rightX.get(frameIndex+frameStart[videoIndex]);
                gtfRightY[f] = aTestAnnotations[videoIndex].rightY.get(frameIndex+frameStart[videoIndex]);

                for (int v : aFrame.leftImg) {
                    writer.println(v);
                }
                writer.flush();
                for (int v : aFrame.rightImg) {
                    writer.println(v);
                }
                writer.flush();

                long startTime = System.currentTimeMillis();
                // call testing function
                // ret = testing(f, aFrame.leftImg, aFrame.rightImg);
                userLeftX[f] = Integer.parseInt(reader.readLine());
                userLeftY[f] = Integer.parseInt(reader.readLine());
                userRightX[f] = Integer.parseInt(reader.readLine());
                userRightY[f] = Integer.parseInt(reader.readLine());
                long endTime = System.currentTimeMillis();

                processingTime += endTime - startTime;

                if (visualize) {
                    drawer.update(leftFilename, rightFilename, gtfLeftX[f],gtfLeftY[f],gtfRightX[f],gtfRightY[f],userLeftX[f],userLeftY[f],userRightX[f],userRightY[f]);
                    if (f==0) drawer.pauseMode = true;
                    drawer.repaint();
                    drawer.processPause();
                    try {
                        Thread.sleep(delay);
                    } catch (Exception e) {
                        // do nothing
                    }
                }

            }

            // calculate score
            correctOutOfBounds(gtfLeftX, gtfLeftY);
            correctOutOfBounds(gtfRightX, gtfRightY);
            
            double l10 = countFrameCorrect(gtfLeftX, gtfLeftY, userLeftX, userLeftY, 10*10);
            double l20 = countFrameCorrect(gtfLeftX, gtfLeftY, userLeftX, userLeftY, 20*20);
            double l50 = countFrameCorrect(gtfLeftX, gtfLeftY, userLeftX, userLeftY, 50*50);
            double r10 = countFrameCorrect(gtfRightX, gtfRightY, userRightX, userRightY, 10*10);
            double r20 = countFrameCorrect(gtfRightX, gtfRightY, userRightX, userRightY, 20*20);
            double r50 = countFrameCorrect(gtfRightX, gtfRightY, userRightX, userRightY, 50*50);
            printMessage("R[10] = " + l10 + " " + r10);
            printMessage("R[20] = " + l20 + " " + r20);
            printMessage("R[50] = " + l50 + " " + r50);
            double accuracyScore = 10000.0 * (50.0*(l10+r10) + 35.0*(l20+r20) + 15.0*(l50+r50));
            System.out.println("Accuracy Score = " + accuracyScore);
            System.out.println("Processing Time = " + processingTime + "ms");

            double timeMultiplier = 1.0;
            double T = 0.001*processingTime;
            if (T>3.33) timeMultiplier = 1.3536 - 0.2939 * Math.log(T);
            if (T>100.0) timeMultiplier = 0;
            timeMultiplier += 1.0;
            System.out.println("Time Multiplier = " + timeMultiplier);

            score = timeMultiplier * accuracyScore;
        } else
        {
            System.out.println("ERROR: Testing file not provided");
            System.exit(0);
        }
        stopSolution();
        return score;
    }

    public static void stopSolution() {
        if (solution != null) {
            try {
                solution.destroy();
            } catch (Exception e) {
                // do nothing
            }
        }
    }

    public static void main(String[] args) throws Exception {

       for (int i = 0; i < args.length; i++) {
            if (args[i].equals("-folder")) {
                folder = args[++i];
            } else if (args[i].equals("-exec")) {
                execCommand = args[++i];
            } else if (args[i].equals("-seed")) {
                seed = Long.parseLong(args[++i]);
            } else if (args[i].equals("-train")) {
                trainingFile = args[++i];
            } else if (args[i].equals("-test")) {
                testingFile = args[++i];
            } else if (args[i].equals("-silent")) {
                debug = false;
            } else if (args[i].equals("-delay")) {
                delay = Integer.parseInt(args[++i]);
            } else if (args[i].equals("-novis")) {
                visualize = false;
            } else {
                System.out.println("WARNING: unknown argument " + args[i] + ".");
            }
        }

        RobotVisionTrackerVis vis = new RobotVisionTrackerVis();
        try {
            double score = vis.doExec();
            System.out.println("Score  = " + score);
        } catch (Exception e) {
            System.out.println("FAILURE: " + e.getMessage());
            e.printStackTrace();
            RobotVisionTrackerVis.stopSolution();
        }
    }

    class ErrorStreamRedirector extends Thread {
        public BufferedReader reader;

        public ErrorStreamRedirector(InputStream is) {
            reader = new BufferedReader(new InputStreamReader(is));
        }

        public void run() {
            while (true) {
                String s;
                try {
                    s = reader.readLine();
                } catch (Exception e) {
                    // e.printStackTrace();
                    return;
                }
                if (s == null) {
                    break;
                }
                System.out.println(s);
            }
        }
    }
}
