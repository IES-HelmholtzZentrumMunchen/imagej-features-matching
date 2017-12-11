//import Jama.Matrix;
import ij.IJ;
import ij.ImageJ;
import ij.ImagePlus;
import ij.plugin.filter.PlugInFilter;
import ij.plugin.frame.RoiManager;
import ij.process.ImageProcessor;

import ij.gui.*;
import ij.*;

import java.awt.*;
import java.awt.geom.Point2D;
import java.util.ArrayList;
import java.util.Vector;
import java.util.concurrent.*;

/**
 * Points of interest tracking
 *
 * A plugin for ImageJ 1.x for points of interest tracking
 *
 * @author Julien Pontabry
 */
public class Features_Matching implements PlugInFilter, DialogListener {
    /**
     * Size of the patch for correlation search (default: 25 pixels).
     */
    private int m_patchSize = 25;

    /**
     * Size of the window for correlation search (default: 100 pixels).
     */
    private int m_windowSize = 70;

    /**
     * Image object.
     */
    private ImagePlus m_image = null;

    /**
     * Time-lapse used for tracking.
     */
    private ImageStack m_timeLapse = null;

    /**
     * Initial position of the points of interest on the first frame.
     */
    private Polygon m_startingPoints = null;

    /**
     * Number of points of interest to track.
     */
    private int m_numberOfPointsOfInterest = 0;

    /**
     * Number of frames in time-lapse.
     */
    private int m_numberOfTimePoints = 0;

    /**
     * Reduction factor of patch sampling (sub-sampling).
     */
    private int m_reductionFactor = 1;

//    /**
//     * True if the trajectory has to be filtered (Kalman filtering), false otherwise.
//     */
//    private boolean m_filterTrajectory = false;
//
//    /**
//     * Std deviation of the input acceleration.
//     */
//    private double m_input_std = 1.0;
//
//    /**
//     * Std deviation of the measurement.
//     */
//    private double m_measurement_std = 10.0;

    /**
     * Normalisation constant of the prior function.
     */
    private double m_normalisation_constant;

    /**
     * Exponent constant of the prior function.
     */
    private double m_exponent_constant;


	/**
	 * @see ij.plugin.filter.PlugInFilter#setup(java.lang.String, ij.ImagePlus)
	 */
	@Override
	public int setup(String arg, ImagePlus imp) {
        // Check if there is an image
        if(imp == null)
        {
            IJ.error("Points of interest tracking", "An image is required for this plugin !");
            return DONE;
        }
        // Check if the image is a time-lapse
        if(imp.getNFrames() <= 1)
        {
            IJ.error("Points of interest tracking", "A time-lapse stack is required for time tracking !");
            return DONE;
        }

        // Check for the number of channels
        if(imp.getNChannels() > 1)
        {
            IJ.error("Points of interest tracking", "Multi-channel is not supported yet !");
            return DONE;
        }

        // Check if there is a ROI
        if(imp.getRoi() == null)
        {
            IJ.error("Points of interest tracking", "A ROI defining points to track is required !");
            return DONE;
        }

        // Get image
        m_image = imp;

        // Get the time-lapse stack
        m_timeLapse = m_image.getStack();

        // Get the points (as Multiple points ROI)
        m_startingPoints = m_image.getRoi().getPolygon();

        // Get the number of time points
        m_numberOfTimePoints = m_image.getNFrames();

        // Get the number of points of interest to track
        m_numberOfPointsOfInterest = m_startingPoints.npoints;

        // Display patch and window on first point
        Roi roiPatch = new Roi(m_startingPoints.xpoints[0] - m_patchSize / 2, m_startingPoints.ypoints[0] - m_patchSize / 2, m_patchSize, m_patchSize);
        roiPatch.setStrokeColor(Color.RED);

        Roi roiWindow = new Roi(m_startingPoints.xpoints[0] - m_windowSize / 2, m_startingPoints.ypoints[0] - m_windowSize / 2, m_windowSize, m_windowSize);
        roiWindow.setStrokeColor(Color.BLUE);

        Overlay previewOverlay = new Overlay();
        previewOverlay.add(roiPatch);
        previewOverlay.add(roiWindow);

        m_image.setOverlay(previewOverlay);

        // Display GUI
        GenericDialog gui = new GenericDialog("Points of interest tracking");
        gui.addNumericField("Patch size", m_patchSize, 0);
        gui.addNumericField("Search window size", m_windowSize, 0);
        gui.addNumericField("Reduction factor", m_reductionFactor, 0);
//        gui.addCheckbox("Filter trajectory", m_filterTrajectory);
//        gui.addNumericField("Acceleration input std", m_input_std, 1);
//        gui.addNumericField("Measurement std", m_measurement_std, 1);
        gui.addNumericField("Prior strength", 0.5, 1);

        gui.addDialogListener(this);
        gui.showDialog();

        // Checking GUI events
        if(gui.wasCanceled())
        {
            m_image.setOverlay(null);
            m_image.setRoi(new PointRoi(m_startingPoints.xpoints, m_startingPoints.ypoints, m_startingPoints.npoints));
            return DONE;
        }

        if(gui.invalidNumber())
        {
            IJ.error("Points of interest tracking", "A numerical field is not filled correctly !");
            return DONE;
        }

        // Get back the parameters from the GUI
        m_patchSize        = (int)gui.getNextNumber();  // Patch size parameter
        m_windowSize       = (int)gui.getNextNumber();  // Search window size parameter
        m_reductionFactor  = (int)gui.getNextNumber();  // Reduction factor of exhaustive search
//        m_filterTrajectory = gui.getNextBoolean();      // Use Kalman filter or not
//        m_input_std        = gui.getNextNumber();       // Input acceleration std deviation
//        m_measurement_std  = gui.getNextNumber();       // Measurement std deviation
        double priorSigma = gui.getNextNumber()*m_windowSize;
        double priorSigma2 =  priorSigma * priorSigma;
        m_normalisation_constant = 1.0 / (2.0 * Math.PI * priorSigma2);
        m_exponent_constant = -0.5 / priorSigma2;

        // Check the parameters
        if(m_patchSize < 3)
        {
            IJ.error("Points of interest tracking", "The patch size should contains at least 3x3 pixels !");
            return DONE;
        }

        if(m_windowSize <= m_patchSize)
        {
            IJ.error("Points of interest tracking", "The search window size should be larger than patch size !");
            return DONE;
        }

        // Return process
		return DOES_8G+DOES_16+DOES_32+ROI_REQUIRED+NO_CHANGES;
	}

	/**
	 * @see ij.plugin.filter.PlugInFilter#run(ij.process.ImageProcessor)
	 */
	@Override
	public void run(ImageProcessor ip) {
//        if(m_filterTrajectory)
//            this.runWithFilter();
//        else // !m_filterTrajectory
            this.runWithoutFilter();
    }

    /**
     * Run the features matching tracking algorithm.
     */
    private void runWithoutFilter() {
        // Create a traces structures to store coordinates of points of interest
        Vector< Point2D[] > pointsTraces = new Vector< Point2D[] >();

        // Display a progress bar
        int currentProgress = 0, maxProgress = m_numberOfPointsOfInterest * m_numberOfTimePoints;
        IJ.showProgress(currentProgress, maxProgress);

        // For each points of interest to track
        for(int pIndex = 0; pIndex < m_numberOfPointsOfInterest; pIndex++)
        {
            // Create a trace structure to store coordinates of current point and initialize it
            Point2D[] pointTrace = new Point2D[m_numberOfTimePoints];
            pointTrace[0]        = new Point2D.Double(m_startingPoints.xpoints[pIndex], m_startingPoints.ypoints[pIndex]);

            // Get initial patch from first frame
            int[] initialExtent = getExtent(pointTrace[0], m_patchSize, m_timeLapse.getWidth(), m_timeLapse.getHeight());

            // Compute the displacement for each time point
            for(int currentFrameIndex = 0; currentFrameIndex < m_numberOfTimePoints-1; currentFrameIndex++)
            {
/*
                // Get patch from current frame around point
                int[] patchExtent = getExtent(pointTrace[currentFrameIndex], m_patchSize, m_timeLapse.getWidth(), m_timeLapse.getHeight());

                // Get search window from next frame around point
                int[] windowExtent = this.getExtent(pointTrace[currentFrameIndex], m_windowSize, m_timeLapse.getWidth(), m_timeLapse.getHeight());

                // Compute normalized cross-correlation of the patch/window couple and get the displacement vector to the maximal point
                Point2D nextPoint = this.getMovedPoint(m_timeLapse.getProcessor(currentFrameIndex + 1), patchExtent, m_timeLapse.getProcessor(currentFrameIndex + 2), windowExtent);

                // Add the location of the point in next frame
                pointTrace[currentFrameIndex+1] = nextPoint;
//*/

//*
                // Search point with previous frame
                // Get patch from current frame around point
                int[] patchExtent = getExtent(pointTrace[currentFrameIndex], m_patchSize, m_timeLapse.getWidth(), m_timeLapse.getHeight());

                // Get search window from next frame around point
                int[] windowExtent = this.getExtent(pointTrace[currentFrameIndex], m_windowSize, m_timeLapse.getWidth(), m_timeLapse.getHeight());

                // Compute normalized cross-correlation of the patch/window couple and get the displacement vector to the maximal point
                Point2D nextPoint = this.getMovedPoint(m_timeLapse.getProcessor(currentFrameIndex + 1), patchExtent, m_timeLapse.getProcessor(currentFrameIndex + 2), windowExtent, pointTrace[currentFrameIndex]);

                // Correcting small drift of point with initial frame
                // Get search window from next frame around point
                int[] windowExtentCorrected = this.getExtent(nextPoint, m_windowSize, m_timeLapse.getWidth(), m_timeLapse.getHeight());

                // Compute normalized cross-correlation of the patch/window couple and get the displacement vector to the maximal point
                Point2D nextPointCorrected = this.getMovedPoint(m_timeLapse.getProcessor(1), initialExtent, m_timeLapse.getProcessor(currentFrameIndex + 2), windowExtentCorrected, pointTrace[currentFrameIndex]);

                Point2D diff = new Point2D.Double(nextPointCorrected.getX()-nextPoint.getX(), nextPointCorrected.getY()-nextPoint.getY());

                if(Math.sqrt(diff.getX()*diff.getX() + diff.getY()*diff.getY()) < 0.05*m_windowSize) {
                    pointTrace[currentFrameIndex+1] = nextPointCorrected;
                }
                else {
                    pointTrace[currentFrameIndex+1] = nextPoint;
                }
//*/

                // Update progress bar
                IJ.showProgress(++currentProgress, maxProgress);
            }

            // Add current trace to traces
            pointsTraces.add(pointTrace);
        }

        // Write points trace to ROI manager
        this.addPointsToRoiManager(pointsTraces);
	}

    private void runWithProbabilisticFilter() {
        // Create a traces structures to store coordinates of points of interest
        Vector< Point2D[] > pointsTraces = new Vector< Point2D[] >();

        // Display a progress bar
        int currentProgress = 0, maxProgress = m_numberOfPointsOfInterest * m_numberOfTimePoints;
        IJ.showProgress(currentProgress, maxProgress);

        // For each points of interest to track
        for(int pIndex = 0; pIndex < m_numberOfPointsOfInterest; pIndex++)
        {
            // Create a trace structure to store coordinates of current point and initialize it
            Point2D[] pointTrace = new Point2D[m_numberOfTimePoints];
            pointTrace[0]        = new Point2D.Double(m_startingPoints.xpoints[pIndex], m_startingPoints.ypoints[pIndex]);

            // Get initial patch from first frame
            int[] initialExtent = getExtent(pointTrace[0], m_patchSize, m_timeLapse.getWidth(), m_timeLapse.getHeight());

            // Compute the displacement for each time point
            for(int currentFrameIndex = 0; currentFrameIndex < m_numberOfTimePoints-1; currentFrameIndex++)
            {
                // Search point with previous frame
                // Get patch from current frame around point
                int[] patchExtent = getExtent(pointTrace[currentFrameIndex], m_patchSize, m_timeLapse.getWidth(), m_timeLapse.getHeight());

                // Get search window from next frame around point
                int[] windowExtent = this.getExtent(pointTrace[currentFrameIndex], m_windowSize, m_timeLapse.getWidth(), m_timeLapse.getHeight());

                // Compute normalized cross-correlation of the patch/window couple and get the displacement vector to the maximal point
                Point2D nextPoint = this.getMovedPoint(m_timeLapse.getProcessor(currentFrameIndex + 1), patchExtent, m_timeLapse.getProcessor(currentFrameIndex + 2), windowExtent, pointTrace[currentFrameIndex]);

                // Correcting small drift of point with initial frame
                // Get search window from next frame around point
                int[] windowExtentCorrected = this.getExtent(nextPoint, m_windowSize, m_timeLapse.getWidth(), m_timeLapse.getHeight());

                // Compute normalized cross-correlation of the patch/window couple and get the displacement vector to the maximal point
                Point2D nextPointCorrected = this.getMovedPoint(m_timeLapse.getProcessor(1), initialExtent, m_timeLapse.getProcessor(currentFrameIndex + 2), windowExtentCorrected, pointTrace[currentFrameIndex]);
                Point2D               diff = new Point2D.Double(nextPointCorrected.getX()-nextPoint.getX(), nextPointCorrected.getY()-nextPoint.getY());

                if(Math.sqrt(diff.getX()*diff.getX() + diff.getY()*diff.getY()) < 0.05*m_windowSize) {
                    pointTrace[currentFrameIndex+1] = nextPointCorrected;
                }
                else {
                    pointTrace[currentFrameIndex+1] = nextPoint;
                }

                // Pour chaque particule
                // |    1. Calcul de la CC dans toute la fenêtre
                // |    2. Récupérer le max (c'est notre observation)
                // |    3. Calculer la constante de normalisation de CC (somme de toutes les valeurs)
                // |    4. Échantillonner la fonction d'importance (CC normalisée - 1 et 3)
                // |    5. Calculer le poids en parcourant la fenêtre et utilisant la formule
                // Fin
                // 6. Normaliser les poids
                // 7. Re-échantillonner si besoin

                // Update progress bar
                IJ.showProgress(++currentProgress, maxProgress);
            }

            // Add current trace to traces
            pointsTraces.add(pointTrace);
        }

        // Write points trace to ROI manager
        this.addPointsToRoiManager(pointsTraces);
    }

//    /**
//     * Run the tracking algorithm with a Kalman filter.
//     */
//    private void runWithFilter() {
//        // Create a traces structures to store coordinates of points of interest
//        Vector< Point2D[] > pointsTraces = new Vector< Point2D[] >();
//
//        // Display a progress bar
//        int currentProgress = 0, maxProgress = m_numberOfPointsOfInterest * m_numberOfTimePoints;
//        IJ.showProgress(currentProgress, maxProgress);
//
//        // For each points of interest to track
//        for(int pIndex = 0; pIndex < m_numberOfPointsOfInterest; pIndex++)
//        {
//            // Create a trace structure to store coordinates of current point and initialize it
//            Point2D[] pointTrace = new Point2D[m_numberOfTimePoints];
//            pointTrace[0]        = new Point2D.Double(m_startingPoints.xpoints[pIndex], m_startingPoints.ypoints[pIndex]);
//
//            // Get initial patch from first frame
//            int[] initialExtent = getExtent(pointTrace[0], m_patchSize, m_timeLapse.getWidth(), m_timeLapse.getHeight());
//
//            // Initialise input acceleration
//            Matrix input_vector = new Matrix(2,1);
//            input_vector.set(0,0, 0.005);
//            input_vector.set(1,0, 0.005);
//
//            // Initialise input state matrix
//            Matrix input_matrix = new Matrix(4,2);
//            input_matrix.set(0,0, 0.5);
//            input_matrix.set(1,1, 0.5);
//            input_matrix.set(2,0, 1);
//            input_matrix.set(3,1, 1);
//
//            // Initialise state vector
//            Matrix corrected_state_vector = new Matrix(4,1);
//            corrected_state_vector.set(0, 0, pointTrace[0].getX());
//            corrected_state_vector.set(1, 0, pointTrace[0].getY());
//            Matrix predicted_state_vector;
//
//            // Initialise state matrices
//            Matrix state_matrix = new Matrix(4,4);
//            state_matrix.set(0,0, 1.0);
//            state_matrix.set(0,2, 1.0);
//            state_matrix.set(1,1, 1.0);
//            state_matrix.set(1,3, 1.0);
//            state_matrix.set(2,2, 1.0);
//            state_matrix.set(3,3, 1.0);
//
//            // Initialise state covariance matrix
//            Matrix state_covariance = new Matrix(4,4);
//            state_covariance.set(0, 0, 0.25);
//            state_covariance.set(0, 2, 0.5);
//            state_covariance.set(1, 1, 0.25);
//            state_covariance.set(1, 3, 0.5);
//            state_covariance.set(2, 0, 0.5);
//            state_covariance.set(2, 2, 1.0);
//            state_covariance.set(3, 1, 0.5);
//            state_covariance.set(3, 3, 1.0);
//            state_covariance.timesEquals(m_input_std * m_input_std);
//            Matrix corrected_state_covariance = state_covariance.copy();
//            Matrix predicted_state_covariance;
//
//            // Initialise measurement vector
//            Matrix measurement_vector = new Matrix(2,1);
//
//            // Initialise measurement matrix
//            Matrix measurement_matrix = new Matrix(2,4);
//            measurement_matrix.set(0,0, 1.0);
//            measurement_matrix.set(1,1, 1.0);
//
//            // Initialise measurement covariance matrix
//            Matrix measurement_covariance = new Matrix(2,2);
//            measurement_covariance.set(0,0, 1.0);
//            measurement_covariance.set(1,1, 1.0);
//            measurement_covariance.timesEquals(m_measurement_std * m_measurement_std);
//
//            // Compute the displacement for each time point
//            for(int currentFrameIndex = 0; currentFrameIndex < m_numberOfTimePoints-1; currentFrameIndex++)
//            {
////                // Get patch from current frame around point
////                int[] patchExtent = getExtent(pointTrace[currentFrameIndex], m_patchSize, m_timeLapse.getWidth(), m_timeLapse.getHeight());
//
//                // Get search window from next frame around point
//                int[] windowExtent = this.getExtent(pointTrace[currentFrameIndex], m_windowSize, m_timeLapse.getWidth(), m_timeLapse.getHeight());
//
//                // Compute normalized cross-correlation of the patch/window couple and get the displacement vector to the maximal point
////                Point2D nextPoint = this.getMovedPoint(m_timeLapse.getProcessor(currentFrameIndex + 1), patchExtent, m_timeLapse.getProcessor(currentFrameIndex + 2), windowExtent);
//                Point2D nextPoint = this.getMovedPoint(m_timeLapse.getProcessor(currentFrameIndex + 1), initialExtent, m_timeLapse.getProcessor(currentFrameIndex + 2), windowExtent);
//
//                // Get measurement for filtering
//                measurement_vector.set(0,0, nextPoint.getX());
//                measurement_vector.set(1,0, nextPoint.getY());
//
//                // Get previous velocity
//                double previous_x_velocity = corrected_state_vector.get(2,0);
//                double previous_y_velocity = corrected_state_vector.get(3,0);
//
//                // Prediction
//                predicted_state_vector   = state_matrix.times(corrected_state_vector).plus(input_matrix.times(input_vector));                       // predict the current state vector
//                predicted_state_covariance = state_matrix.times(corrected_state_covariance).times(state_matrix.transpose()).plus(state_covariance); // predict the current state covariance
//
//                // Correction
//                Matrix       gain_matrix = predicted_state_covariance.times(measurement_matrix.transpose()).times( (measurement_matrix.times(predicted_state_covariance).times(measurement_matrix.transpose()).plus(measurement_covariance)).inverse() );   // Compute the kalman gain
//                corrected_state_vector     = predicted_state_vector.plus(gain_matrix.times( measurement_vector.minus(measurement_matrix.times(predicted_state_vector)) ));                                                                                  // Compute the corrected state
//                corrected_state_covariance = predicted_state_covariance.minus(gain_matrix.times(measurement_matrix).times(predicted_state_covariance));                                                                                                     // Compute the corrected state covariance
//
//                // Correct measurement
//                nextPoint = new Point2D.Double(corrected_state_vector.get(0,0),corrected_state_vector.get(1,0));
//
////                // Estimate current acceleration
////                input_vector.set(0, 0, corrected_state_vector.get(2, 0) - previous_x_velocity);
////                input_vector.set(1,0, corrected_state_vector.get(3,0) - previous_y_velocity);
//
//                // Add the location of the point in next frame
//                pointTrace[currentFrameIndex+1] = nextPoint;
//
//                // Update progress bar
//                IJ.showProgress(++currentProgress, maxProgress);
//            }
//
//            // Add current trace to traces
//            pointsTraces.add(pointTrace);
//        }
//
//        // Write points trace to ROI manager
//        this.addPointsToRoiManager(pointsTraces);
//    }

    /**
     * Get extent of a patch given a point, size of patch and size of image.
     * The extent is defined by four integers : min and max values of x coordinate
     * and min and max values of y coordinate.
     * @param point Center point of the patch
     * @param patchSize Size of the patch
     * @param width Width of the image
     * @param height Height of the image
     * @return An array containing four integers defining the extent
     */
    private int[] getExtent(Point2D point, int patchSize, int width, int height)
    {
        // Initialize extent
        int[] extent = new int[4];

        // Fill it
        extent[0] = (int)point.getX() - patchSize/2;
        extent[1] = (int)point.getX() + patchSize/2;
        extent[2] = (int)point.getY() - patchSize/2;
        extent[3] = (int)point.getY() + patchSize/2;

        // Check boundary and crop if necessary
        if(extent[0] < 0)
            extent[0] = 0;

        if(extent[1] >= width)
            extent[1] = width;

        if(extent[2] < 0)
            extent[2] = 0;

        if(extent[3] >= height)
            extent[3] = height;

        // Return the extent
        return extent;
    }

//    /**
//     * Get the new location of a point moved by the tracking process.
//     * @param currentFrameIp Processor of the current frame
//     * @param patchExtent Extent of the patch around point of interest
//     * @param nextFrameIp Processor of the next frame
//     * @param windowExtent Extent of the search window
//     * @return New point moved by tracking process
//     */
//    private Point2D getMovedPoint(ImageProcessor currentFrameIp, int[] patchExtent, ImageProcessor nextFrameIp, int[] windowExtent)
//    {
//        // Initialize the displacement vector to the null vector
//        Point2D movedPoint = new Point2D.Double();
//
//        // Offsets are defined for the moving patch to the window (in window space for convenience)
//        int  patchSizeX = patchExtent[1]  - patchExtent[0],   patchSizeY = patchExtent[3]  - patchExtent[2];
//        int windowSizeX = windowExtent[1] - windowExtent[0], windowSizeY = windowExtent[3] - windowExtent[2];
//        int  maxOffsetX = windowSizeX - patchSizeX;
//        int  maxOffsetY = windowSizeY - patchSizeY;
//
//        // We search for the maximal cross-correlation
//        // So initialize search variable to 0
//        double maxCrossCorrelation = 0;
//
//        // For each offsets between image and template
//        for(int offsetX = 0; offsetX < maxOffsetX; offsetX++)
//        {
//            for(int offsetY = 0; offsetY < maxOffsetY; offsetY++)
//            {
//                // Compute coordinates in window space
//                int minX = windowExtent[0] + offsetX, maxX = minX + patchSizeX;
//                int minY = windowExtent[2] + offsetY, maxY = minY + patchSizeY;
//
//                // Compute normalized cross-correlation and prior between image and template + offset
//                int[] movedPatchExtent = { minX, maxX, minY, maxY };
//                double currentCrossCorrelation = this.normalizedCrossCorrelation(currentFrameIp, patchExtent, nextFrameIp, movedPatchExtent);
//
//                // Seek for maximal correlation
//                if(currentCrossCorrelation > maxCrossCorrelation)
//                {
//                    maxCrossCorrelation = currentCrossCorrelation;
//                    movedPoint.setLocation(minX, minY);
//                }
//            }
//        }
//
//        // Displacement correction because we want the displacement from the center in image space
//        movedPoint.setLocation(movedPoint.getX() + patchSizeX/2, movedPoint.getY() + patchSizeX/2);
//
//
//        return movedPoint;
//    }

    private double priorFunction(Point2D previousPoint, Point2D nextPoint) {
//        double sigma = m_windowSize/2.0;
//        double sigma2 = sigma * sigma;
//        double normalisation_constant = 1.0 / (2.0 * Math.PI * sigma2);
//        double exponent_constant = -0.5 / sigma2;

        double diffX = previousPoint.getX() - nextPoint.getX();
        double diffY = previousPoint.getY() - nextPoint.getY();

        return m_normalisation_constant * Math.exp( m_exponent_constant * (diffX*diffX + diffY*diffY) );
    }

    /**
     * Get the new location of a point moved by the tracking process.
     * This method is multi-threaded.
     * @param currentFrameIp Processor of the current frame
     * @param patchExtent Extent of the patch around point of interest
     * @param nextFrameIp Processor of the next frame
     * @param windowExtent Extent of the search window
     * @return New point moved by tracking process
     */
    private Point2D getMovedPoint(ImageProcessor currentFrameIp, int[] patchExtent, ImageProcessor nextFrameIp, int[] windowExtent, Point2D previousPoint)
    {
        // Initialize the displacement vector to the null vector
        Point2D movedPoint = new Point2D.Double();

        // Offsets are defined for the moving patch to the window (in window space for convenience)
        int  patchSizeX = patchExtent[1]  - patchExtent[0],   patchSizeY = patchExtent[3]  - patchExtent[2];
        int windowSizeX = windowExtent[1] - windowExtent[0], windowSizeY = windowExtent[3] - windowExtent[2];
        int  maxOffsetX = windowSizeX - patchSizeX;
        int  maxOffsetY = windowSizeY - patchSizeY;

        // Decenter the previous point
        Point2D nonCenteredPreviousPoint = new Point2D.Double(previousPoint.getX() - patchSizeX/2, previousPoint.getY() - patchSizeY/2);

        // Create parallel inputs
        ArrayList< ParallelInput > inputs = new ArrayList< ParallelInput >();

        for(int offsetX = 0; offsetX < maxOffsetX; offsetX++) {
            for(int offsetY = 0; offsetY < maxOffsetY; offsetY++) {
                ParallelInput input = new ParallelInput();

                // Private data
                input.offsetX = offsetX;
                input.offsetY = offsetY;

                // Shared data
                input.patchExtent    = patchExtent;
                input.windowExtent   = windowExtent;
                input.currentFrameIp = currentFrameIp;
                input.nextFrameIp    = nextFrameIp;
                input.parent         = this;
                input.patchSizeX     = patchSizeX;
                input.patchSizeY     = patchSizeY;

                inputs.add(input);
            }
        }

        try {
            // Threaded computations
            ArrayList< ParallelOutput > outputs = processInputs(inputs);

/*
            // Search for point with maximal correlation
            double maxCrossCorrelation = 0.0;

            for(final ParallelOutput output : outputs) {
                if(maxCrossCorrelation < output.correlation) {
                    maxCrossCorrelation = output.correlation;
                    movedPoint = output.movedPoint;
                }
            }
//*/

//*
            // Search for point with maximal correlation
            double maxScore = 0.0;

            for(final ParallelOutput output : outputs) {
                double score = output.correlation * this.priorFunction(nonCenteredPreviousPoint, output.movedPoint);

                if(maxScore < score) {
                    maxScore = score;
                    movedPoint = output.movedPoint;
                }
            }
//*/
        } catch (Exception e){
            IJ.showMessage(e.getMessage());
        }

        // Displacement correction because we want the displacement from the center in image space
        movedPoint.setLocation(movedPoint.getX() + patchSizeX/2, movedPoint.getY() + patchSizeY/2);


        return movedPoint;
    }

    /**
     * Simple for loop parallel processing.
     * @param inputs Input parameters (local and not shared)
     * @return Output parameters
     * @throws InterruptedException
     * @throws java.util.concurrent.ExecutionException
     */
    private ArrayList< ParallelOutput > processInputs(ArrayList< ParallelInput > inputs) throws InterruptedException, ExecutionException {
        // Create a service for threads (same number as available processors)
        int threads = Runtime.getRuntime().availableProcessors() + 1;
        ExecutorService service = Executors.newFixedThreadPool(threads);

        // Create a list of future outputs
        ArrayList< Future< ParallelOutput > > futures = new ArrayList< Future< ParallelOutput > >();

        // For each input, do a parallel process using service
        for(final ParallelInput input : inputs) {
            // Create a callable defining the process
            final Callable< ParallelOutput > callable = new Callable< ParallelOutput >() {
                public ParallelOutput call() throws Exception {
                    // Allocate memory for output
                    ParallelOutput output = new ParallelOutput();

                    // Compute coordinates in window space
                    int minX = input.windowExtent[0] + input.offsetX, maxX = minX + input.patchSizeX;
                    int minY = input.windowExtent[2] + input.offsetY, maxY = minY + input.patchSizeY;

                    // Compute normalized cross-correlation and prior between image and template + offset
                    int[] movedPatchExtent = { minX, maxX, minY, maxY };
                    output.correlation = input.parent.normalizedCrossCorrelation(input.currentFrameIp, input.patchExtent, input.nextFrameIp, movedPatchExtent);

                    // Set output point
                    output.movedPoint = new Point2D.Double(minX, minY);

                    return output;
                }
            };

            // Add future results as outputs
            futures.add(service.submit(callable));
        }

        // Stop the threads service
        service.shutdown();

        // Create an output list
        ArrayList< ParallelOutput > outputs = new ArrayList< ParallelOutput >();

        // Fill the list with the results
        for(Future<ParallelOutput> future : futures) {
            outputs.add(future.get());
        }

        return outputs;
    }

    /**
     * Input data for parallel processing.
     * The offsets are private and the other are shared.
     */
    private class ParallelInput {
        // Private data
        public int offsetX;
        public int offsetY;

        // Shared data
        public int[] patchExtent;
        public int[] windowExtent;
        public ImageProcessor currentFrameIp;
        public ImageProcessor nextFrameIp;
        public Features_Matching parent;
        public int patchSizeX;
        public int patchSizeY;
    }

    /**
     * Output data for parallel processing.
     */
    private class ParallelOutput {
        public double correlation;
        public Point2D movedPoint;
    }

    /**
     * Compute the normalized cross-correlation of a template and an image.
     * Both extent of template and image should define the same number of pixels.
     * Note that due to processing saving, not check of the number of pixels is done.
     * @param templateIp Processor of the template
     * @param templateExtent Extent of the template
     * @param imageIp Processor of the image
     * @param imageExtent Extent of the image
     * @return The normalized cross-correlation value
     */
    private double normalizedCrossCorrelation(ImageProcessor templateIp, int[] templateExtent, ImageProcessor imageIp, int[] imageExtent)
    {
        // Initialize correlation
        double correlation = 0;

        // NOTE Check for the size of the template and images (should be the same, not checked because of processing saving)
        // Get the number of samples
        double numberOfSamples = (templateExtent[1]-templateExtent[0]) * (templateExtent[3]-templateExtent[2]) / (m_reductionFactor * m_reductionFactor);

        // Initialize means and std deviations
        double templateMean = 0.0, templateStd = 1.0;
        double    imageMean = 0.0,    imageStd = 1.0;

        // Compute the means
        for(int templateX = templateExtent[0], imageX = imageExtent[0]; templateX <= templateExtent[1] && imageX <= imageExtent[1]; templateX+=m_reductionFactor, imageX+=m_reductionFactor)
        {
            for(int templateY = templateExtent[2], imageY = imageExtent[2]; templateY <= templateExtent[3] && imageY <= imageExtent[3]; templateY+=m_reductionFactor, imageY+=m_reductionFactor)
            {
                templateMean += templateIp.getf(templateX,templateY);
                imageMean    += imageIp.getf(imageX,imageY);
            }
        }

        templateMean /= numberOfSamples;
        imageMean    /= numberOfSamples;

        // Compute the sum of products and the std deviations
        for(int templateX = templateExtent[0], imageX = imageExtent[0]; templateX <= templateExtent[1] && imageX <= imageExtent[1]; templateX+=m_reductionFactor, imageX+=m_reductionFactor)
        {
            for(int templateY = templateExtent[2], imageY = imageExtent[2]; templateY <= templateExtent[3] && imageY <= imageExtent[3]; templateY+=m_reductionFactor, imageY+=m_reductionFactor)
            {
                // Compute the differences to the means
                double templateDiff = templateIp.getf(templateX,templateY) - templateMean;
                double    imageDiff = imageIp.getf(imageX,imageY)          - imageMean;

                // Sum up the products
                correlation += templateDiff * imageDiff;

                // Compute std deviations
                templateStd += templateDiff * templateDiff;
                imageStd    += imageDiff    * imageDiff;
            }
        }

        templateStd = Math.sqrt(templateStd / (numberOfSamples-1));
        imageStd    = Math.sqrt(imageStd    / (numberOfSamples-1));


        // Normalize cross-correlation
        correlation /= numberOfSamples * templateStd * imageStd;


        return  correlation;
    }

    /**
     * Add tracking points to the ROI Manager
     * @param pointsTraces Traces of all points.
     */
    private void addPointsToRoiManager(Vector< Point2D[] > pointsTraces)
    {
        // Create the roi manager
        RoiManager manager = new RoiManager();

        // For each time point
        for(int t = 0; t < m_numberOfTimePoints; t++)
        {
            // Create arrays of points
            float[] coordinatesX = new float[m_numberOfPointsOfInterest];
            float[] coordinatesY = new float[m_numberOfPointsOfInterest];

            // For each point of interest
            for(int pointsIndex = 0; pointsIndex < m_numberOfPointsOfInterest; pointsIndex++)
            {
                // Add current point coordinates to arrays
                coordinatesX[pointsIndex] = (float)pointsTraces.get(pointsIndex)[t].getX();
                coordinatesY[pointsIndex] = (float)pointsTraces.get(pointsIndex)[t].getY();
            }

            // Create a point roi for the current time
            PointRoi roi = new PointRoi(coordinatesX, coordinatesY, m_numberOfPointsOfInterest);

            // Add the roi the manager
            m_image.setSlice(t + 1);
            manager.add(m_image, roi, t);
        }
    }

    /**
     * Event listener of the dialog gui.
     * @see ij.gui.DialogListener#dialogItemChanged(ij.gui.GenericDialog, java.awt.AWTEvent)
     * @param gui Graphical user interface object who raise an event
     * @param e Event raised
     * @return True if the parameters are allowed, false otherwise
     */
    @Override
    public boolean dialogItemChanged(GenericDialog gui, AWTEvent e)
    {
        boolean dialogOk = true;

        if(e != null && e.paramString().equals("TEXT_VALUE_CHANGED"))
        {
            // Get fields
            int  patchSize = (int)gui.getNextNumber();
            int windowSize = (int)gui.getNextNumber();

            // Activate or deactivate Ok button
            if(patchSize < 3 || windowSize <= patchSize)
            {
                dialogOk = false;
            }

            // Display on image an example of patch and window on the first point
            Roi roiPatch = new Roi(m_startingPoints.xpoints[0] - patchSize / 2, m_startingPoints.ypoints[0] - patchSize / 2, patchSize, patchSize);
            roiPatch.setStrokeColor(Color.RED);

            Roi roiWindow = new Roi(m_startingPoints.xpoints[0] - windowSize / 2, m_startingPoints.ypoints[0] - windowSize / 2, windowSize, windowSize);
            roiWindow.setStrokeColor(Color.BLUE);

            Overlay previewOverlay = new Overlay();
            previewOverlay.add(roiPatch);
            previewOverlay.add(roiWindow);

            m_image.setOverlay(previewOverlay);
        }
        else if(e == null)
        {
            m_image.setRoi(new PointRoi(m_startingPoints.xpoints, m_startingPoints.ypoints, m_startingPoints.npoints));
            m_image.setOverlay(null);
        }
        // else any other event

        return dialogOk;
    }

	/**
	 * Main method for debugging.
	 *
	 * For debugging, it is convenient to have a method that starts ImageJ, loads an
	 * image and calls the plugin, e.g. after setting breakpoints.
	 *
	 * @param args unused
	 */
	public static void main(String[] args) {
		// set the plugins.dir property to make the plugin appear in the Plugins menu
		Class<?> clazz = Features_Matching.class;
		String url = clazz.getResource("/" + clazz.getName().replace('.', '/') + ".class").toString();
		String pluginsDir = url.substring(5, url.length() - clazz.getName().length() - 6);
		System.setProperty("plugins.dir", pluginsDir);

		// start ImageJ
		new ImageJ();
	}
}
