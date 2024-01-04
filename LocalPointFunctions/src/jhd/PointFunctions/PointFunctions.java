package jhd.PointFunctions;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Random;

import jhd.ProgressBars.*;

/**Experimental Class for Point probability functions
 * @author John
 *
 */
public class PointFunctions
{
	/**Default Constructor*/
	public PointFunctions() {}

	/**A structure describing a 3D pixel
	 * @author John
	 *
	 */
	public class Point3D
	{
		/**The x or i or column location of the pixel*/
		public int x;
		/**The y or j or row location of the pixel*/
		public int y;
		/**The z or k or slice location of the pixel*/
		public int z;
		/**The x dimension of the pixel*/
		public double pixelWidth;
		/**The y dimension of the pixel*/
		public double pixelHeight;
		/**The z dimension of the pixel*/
		public double pixelDepth;
		/**The dimension unit e &period;g &period; cm or inch etc. of the pixel*/
		public String unit;
		/**The value at the pixel(x,y,z) location*/
		public double val;

		/**Default constructor*/
		public Point3D() {}
		/**
		 * @param x The x or horizontal or column position of a pixel
		 * @param y The y or vertical or row position of a pixel
		 * @param z The z or slice or depth position of a pixel
		 */
		public Point3D(int x, int y, int z) { this.x= x; this.y= y; this.z = z;}
		/**
		 * @param x The x or horizontal or column position of a pixel
		 * @param y The y or vertical or row position of a pixel
		 * @param z The z or slice or depth position of a pixel
		 * @param val The value of the pixel at the coordinates (x,y,z)
		 */
		public Point3D(int x, int y, int z, double val) { this.x= x; this.y= y; this.z = z; this.val= val; }
		/**Constructor*/
		/**
		 * @param x The x or horizontal or column position of a pixel
		 * @param y The y or vertical or row position of a pixel
		 * @param z The z or slice or depth position of a pixel
		 * @param pixelWidth The pixel dimension in the x direction
		 * @param pixelHeight The pixel dimension in the y direction
		 * @param pixelDepth The pixel dimension in the z direction
		 * @param val The value of the pixel at the coordinates (x,y,z)
		 */
		public Point3D(int x, int y, int z, double pixelWidth, double pixelHeight, double pixelDepth, double val)
		{
			this.x= x; this.y= y; this.z = z; 
			this.pixelWidth= pixelWidth; this.pixelHeight= pixelHeight; this.pixelDepth = pixelDepth;
			this.val= val;
		}
	}

	/**Class for holding histogram bin andvalue data
	 * @author John
	 *
	 */
	public class Histogram
	{
		/**The histogram bins data*/
		public double[] bin;
		/**The histogram value data*/
		public double[] count;
		/**Default constructor*/
		public Histogram() {}
		
		/**
		 * @param bin The bin value for the histogram
		 * @param count The count value in the corresponding bin
		 */
		public Histogram(double[] bin, double[] count)
		{
			super();
			/**The histogram bins data*/
			this.bin = bin;
			/**The histogram value data*/
			this.count = count;
		}
	}

	Random myRandom;

	//*******************************************************************************
	
	/**Required call, initializes the random generator using a seed
	 * drawn from the last 4 digits of the millisecond system clock.
	 */
	public void initRandom()
	{
		myRandom = new Random();
		long msec = System.currentTimeMillis();
		msec = msec % 10000;
		myRandom.setSeed(msec);		
	}

	/**Gets the coordinates and voxel values along a Bresenham line between p1 and p2 in an image
	 * @param image the image data as an array of 2D slices. Pass z=0 for a 2D image.
	 * @param imgWidth The width of the image in pixels
	 * @param imgHeight The height of the image in pixels
	 * @param imgDepth The depth of the image in pixels
	 * @param start A Point3D descriptor of the x,y,z point start of the line
	 * @param end A Point3D descriptor of the x,y,z point end of the line
	 * @return A Point3D Array containing the pixel coordinates and nearest neighbor values along the line
	 */	
	public  Point3D[] bresenhamLine(Object[] image, int imgWidth, int imgHeight, int imgDepth,Point3D start, Point3D end)
	{
		Point3D[]	theLine =  bresenhamLine(start,end);
		getLineValues(image,imgWidth,imgHeight,imgDepth,theLine);
		return theLine;
	}

	//*******************************************************************************

	/**Gets the coordinates along a Bresenham line between p1 and p2
	 * @param startPoint A Point3D descriptor of the x,y,z point start of the line
	 * @param endPoint A Point3D descriptor of the x,y,z point end of the line
	 * @return A Point3D Array containing the pixel coordinates along the line
	 */	
	public Point3D[] bresenhamLine(Point3D startPoint, Point3D endPoint)
	{
		// Bresenham code contributed by ishankhandelwals.
		// https://www.geeksforgeeks.org/bresenhams-algorithm-for-3-d-line-drawing/
		// Modified by LazzyIzzi to return the line coordinates
		ArrayList<Point3D>	theLine = new ArrayList<Point3D>();
		
		//make a local copy of the  input points
		Point3D start= new Point3D(startPoint.x,startPoint.y,startPoint.z);
		Point3D end = new Point3D(endPoint.x,endPoint.y,endPoint.z);

		theLine.add(new Point3D(start.x,start.y,start.z));

		int dx = Math.abs(end.x - start.x);
		int dy = Math.abs(end.y - start.y);
		int dz = Math.abs(end.z - start.z);
		int xs;
		int ys;
		int zs;
		if (end.x > start.x) {
			xs = 1;
		} else {
			xs = -1;
		}
		if (end.y > start.y) {
			ys = 1;
		} else {
			ys = -1;
		}
		if (end.z > start.z) {
			zs = 1;
		} else {
			zs = -1;
		}

		// Driving axis is X-axis"
		if (dx >= dy && dx >= dz) {
			int p1 = 2 * dy - dx;
			int p2 = 2 * dz - dx;
			while (start.x != end.x) {
				start.x += xs;
				if (p1 >= 0) {
					start.y += ys;
					p1 -= 2 * dx;
				}
				if (p2 >= 0) {
					start.z += zs;
					p2 -= 2 * dx;
				}
				p1 += 2 * dy;
				p2 += 2 * dz;

				theLine.add(new Point3D(start.x,start.y,start.z));
			}			
		}

		// Driving axis is Y-axis"
		else if (dy >= dx && dy >= dz)
		{
			int p1 = 2 * dx - dy;
			int p2 = 2 * dz - dy;
			while (start.y != end.y)
			{
				start.y += ys;
				if (p1 >= 0) {
					start.x += xs;
					p1 -= 2 * dy;
				}
				if (p2 >= 0)
				{
					start.z += zs;
					p2 -= 2 * dy;
				}
				p1 += 2 * dx;
				p2 += 2 * dz;
				theLine.add(new Point3D(start.x,start.y,start.z));
			}			
		}

		// Driving axis is Z-axis"
		else
		{
			int p1 = 2 * dy - dz;
			int p2 = 2 * dx - dz;
			while (start.z != end.z)
			{
				start.z += zs;
				if (p1 >= 0)
				{
					start.y += ys;
					p1 -= 2 * dz;
				}
				if (p2 >= 0)
				{
					start.x += xs;
					p2 -= 2 * dz;
				}
				p1 += 2 * dy;
				p2 += 2 * dx;
				theLine.add(new Point3D(start.x,start.y,start.z));
			}
		}
		Point3D[] lineArr = theLine.toArray(new Point3D[theLine.size()]);
		return lineArr;
	}

	//*********************************************************************************

	/**Computes the chord length distribution as a probability density, 
	 * <br>i &period;e &period; the likelihood of finding a chord in the interval between L +/- dL in image units
	 * This code measures the length of chords through the selected component along 
	 * <br>a random line of length = (min image dimension)/2. The chord collection is then binned by length.
	 * <br>Anisotropic pixels and voxels are supported. 
	 * @param image 1D Reference to a binarized 3D short image
	 * @param imgWidth The image width in pixels
	 * @param imgHeight The height width in pixels
	 * @param imgDepth The depth width in pixels
	 * @param pixelWidth the size of the voxel in the x, width, column, direction (e &period;g &period; cm , inch etc.)
	 * @param pixelHeight the size of the voxel in the y, imgHeight, row, direction (e &period;g &period; cm , inch etc.)
	 * @param pixelDepth the size of the voxel in the z, width, slice, direction (e &period;g &period; cm , inch etc.)
	 * @param nBins The number of bins in the histogram
	 * @param nSamples The number of samples for each line length
	 * @param valChoice use "Map 0" to measure the zero component and "Map !0" for the non-zero component
	 * @param showProgress true to display a progress bar
	 * @return Reference to a histogram of length (in user units ) vs length
	 */
	public Histogram chordLengthDistribution(Object[] image,int imgWidth,int imgHeight,int imgDepth,
			double pixelWidth, double pixelHeight, double pixelDepth,
			int nBins, int nSamples, String valChoice, boolean showProgress)
	{
		ProgressBars pBars = null;;
		if(showProgress)
		{	
			pBars = new ProgressBars("Chord Length Distribution");
			pBars.add("Chord Length", 400, 20, 0, nSamples);
			pBars.show();
		}

		Point3D	p1 = new Point3D();
		Point3D	p2 = new Point3D();
		float maxChord = Float.MIN_VALUE;
		float minChord = Float.MAX_VALUE;

		ArrayList<Float> allChords = new ArrayList<Float>();

		int minR = getMaxPixelSeparationDistance(imgWidth,imgHeight,imgDepth);
		Point3D[] theLine;
		float[] chordsLengths;

		double minScale = pixelWidth;
		if(minScale> pixelHeight) minScale=pixelHeight;
		if(minScale> pixelDepth) minScale=pixelDepth;
		
		double xScale = pixelWidth/minScale;
		double yScale = pixelHeight/minScale;
		double zScale = pixelDepth/minScale;

		for(int i=0;i<nSamples;i++)
		{
			randomPointPair(imgWidth,imgHeight,imgDepth,xScale,yScale,zScale,minR, p1,p2);

			theLine = bresenhamLine(image,imgWidth,imgHeight,imgDepth,p1,p2);

			chordsLengths = getChords(theLine,valChoice,pixelWidth,pixelHeight,pixelDepth);

			for(int j =0;j<chordsLengths.length;j++)
			{
				allChords.add(chordsLengths[j]);
				if(chordsLengths[j]>maxChord) maxChord = chordsLengths[j];
				if(chordsLengths[j]<minChord) minChord = chordsLengths[j];
			}
			if(showProgress) pBars.setValue("Chord Length",i);
		}							

		//linear scales
		int maxBin=nBins-1;
		int minBin=0;
		float m=(maxChord-minChord)/(maxBin-minBin);
		float b = minChord;

		double[] bin = new double[nBins];
		for(int j = 0;j<nBins;j++) bin[j] = m*j+b;

		double[] count = new double[nBins];
		for(int j =0;j<allChords.size();j++)
		{
			int binIndex = (int)((allChords.get(j)-b)/m);
			count[binIndex]++;
		}
		
		//Express as Probability density
		//i.e the probability a chord of length L
		//in the interval L1-L2
		for(int j =0;j<nBins;j++)
		{
			count[j]/=allChords.size();
		}

		Histogram hist= new Histogram(bin,count);		
		if(showProgress) pBars.close();
		return hist;

	}

	//*********************************************************************************

	/**Computes the chord length distribution as a probability density, 
	 * <br>i &period;e &period; the likelihood of finding a chord in the interval between L and dL in pixels or voxels
	 * This code measures the length of chords through the selected component along 
	 * <br>a random line of length = (min image dimension)/2. The chord collection is then binned by length.
	 * @param image 1D Reference to a binarized 3D short image
	 * @param imgWidth The image width in pixels
	 * @param imgHeight The height width in pixels
	 * @param imgDepth The depth width in pixels
	 * @param nBins The number of bins in the histogram
	 * @param nSamples The number of samples for each line length
	 * @param valChoice use "Map 0" to measure the zero component and "Map !0" for the non-zero component
	 * @param showProgress true to display a progress bar
	 * @return Reference to a histogram of length (in pixels) vs length
	 */
	public Histogram chordLengthDistribution(Object[] image,int imgWidth,int imgHeight,int imgDepth, int nBins, int nSamples, String valChoice, boolean showProgress)
	{
		return chordLengthDistribution(image, imgWidth, imgHeight, imgDepth,1.0, 1.0, 1.0, nBins, nSamples, valChoice, showProgress);
	}

	//*********************************************************************************

	/**Gets the coordinates and voxel values along a Bresenham line between p1 and p2 in an image and adds them to the Point3D val parameter.
	 * @param image the image data as an array of 2D slices. Pass z=0 for a 2D image.
	 * @param imgWidth The width of the image in pixels
	 * @param imgHeight The height of the image in pixels
	 * @param imgDepth The depth of the image in pixels
	 * @param line A Point3D descriptor of the x,y,z points of the line
	 */	
	public void getLineValues(Object[] image, int imgWidth, int imgHeight, int imgDepth, Point3D[] line)
	{

		if(image[0] instanceof byte[])
		{
			for(Point3D pt : line)
			{
				byte[] slice = (byte[])image[pt.z];
				pt.val= (double)slice[pt.x + pt.y*imgWidth];
			}
		}
		else if(image[0] instanceof short[])
		{
			for(Point3D pt : line)
			{
				short[] slice = (short[])image[pt.z];
				pt.val= (double)slice[pt.x + pt.y*imgWidth];
			}
		}
		else if(image[0] instanceof int[])
		{
			for(Point3D pt : line)
			{
				int[] slice = (int[])image[pt.z];
				pt.val= (double)slice[pt.x + pt.y*imgWidth];
			}
		}
		else if(image[0] instanceof float[])
		{
			for(Point3D pt : line)
			{
				float[] slice = (float[])image[pt.z];
				pt.val= (double)slice[pt.x + pt.y*imgWidth];
			}
		}
		else if(image[0] instanceof double[])
		{
			for(Point3D pt : line)
			{
				double[] slice = (double[])image[pt.z];
				pt.val= (double)slice[pt.x + pt.y*imgWidth];
			}
		}
	}

	/**The probe box is a cube in the image's physical dimensions, not pixel dimensions.
	 * <br>The probe box must be small enough to allow many samples each image dimension,
	 * <br>and larger than the typical pore size.
	 * @param imgWidth The image width in pixels
	 * @param imgHeight The height width in pixels
	 * @param imgDepth The depth width in pixels
	 * @param pixelWidth the size of the voxel in the x, width, column, direction (e &period;g &period; cm , inch etc.)
	 * @param pixelHeight the size of the voxel in the y, imgHeight, row, direction (e &period;g &period; cm , inch etc.)
	 * @param pixelDepth the size of the voxel in the z, width, slice, direction (e &period;g &period; cm , inch etc.)
	 * @return One fourth of the smallest image physical dimension
	 */
	public double getMaxProbeSize(int imgWidth,int imgHeight,int imgDepth, double pixelWidth,double pixelHeight,double pixelDepth)
	{
		double imgWidthDim =imgWidth*pixelWidth;
		double imgHeightDim =imgHeight*pixelHeight;
		double imgDepthDim =imgDepth*pixelDepth;
		double minDim;
		if(imgDepth==1)
		{
			minDim = imgWidthDim;
			if(minDim>imgHeightDim) minDim=imgHeightDim;			
		}
		else
		{
			minDim = imgWidthDim;
			if(minDim>imgHeightDim) minDim = imgHeightDim;
			if(minDim>imgDepthDim) minDim = imgDepthDim;
		}
		return minDim/4.0;		
	}
	//*********************************************************************************

	/**Returns the smallest image dimension divided by 2
	 * @param imgWidth the image width
	 * @param imgHeight the image height
	 * @param imgDepth the image depth
	 * @return The smallest image dimension divided by 2
	 * <br> For 2D images (imgDepth=1) only width and imgHeight are compared.
	 */
	public int getMaxPixelSeparationDistance(int imgWidth,int imgHeight,int imgDepth)
	{
		//find the smallest image dimension
		int minSize=Integer.MAX_VALUE;
		if(imgDepth==1) //single slice
		{			
			if(minSize>imgWidth) minSize=imgWidth;
			if(minSize>imgHeight) minSize = imgHeight;			
		}
		else
		{
			if(minSize>imgDepth) minSize=imgDepth;
			if(minSize>imgWidth) minSize=imgWidth;
			if(minSize>imgHeight) minSize = imgHeight;			
		}
		return minSize/2;		
	}

	//*********************************************************************************

	/**	Find the probability that a random line of length R in image units will lie completely in the chosen image space.
	 * Note 1. Distances are scanned in integer pixel steps.
	 * Note 2. The longest distance sampled is half the smallest image dimension.
	 * @param image  A reference to an 2D or 3D image segmented into zero and non-zero value regions
	 * @param imgWidth The image width in pixels
	 * @param imgHeight The image height in pixels
	 * @param imgDepth The image depth in pixels, pass 1 for 2D images.
	 * @param pixelWidth the size of the voxel in the x, width, column, direction (e &period;g &period; cm , inch etc.)
	 * @param pixelHeight the size of the voxel in the y, imgHeight, row, direction (e &period;g &period; cm , inch etc.)
	 * @param pixelDepth the size of the voxel in the z, width, slice, direction (e &period;g &period; cm , inch etc.), ignored for 2D images
	 * @param nSamples The number of samples at each distance step. More samples = less noise
	 * @param valChoice "Map 0" or "Map !0".
	 * @param showProgress Display progress Bar while running
	 * @return Reference to a Histogram of probability vs length in image units.
	 */
	public Histogram linealPathDistribution(Object[] image,int imgWidth,int imgHeight, int imgDepth,
			double pixelWidth, double pixelHeight, double pixelDepth,
			int nSamples, String valChoice, boolean showProgress)
	{
		ProgressBars pBars = null;;
		int maxR = getMaxPixelSeparationDistance(imgWidth,imgHeight,imgDepth);
		if(showProgress)
		{	
			pBars = new ProgressBars("Lineal Path Distribution");
			pBars.add("Length", 400, 20, 0, maxR);
			pBars.show();
		}
		Point3D	p1 = new Point3D();
		Point3D	p2 = new Point3D();
		Point3D p1a = new Point3D();
		Point3D p2a = new Point3D();

		double minScale = pixelWidth;
		if(minScale> pixelHeight) minScale=pixelHeight;
		if(minScale> pixelDepth) minScale=pixelDepth;
		
		double xScale = pixelWidth/minScale;
		double yScale = pixelHeight/minScale;
		double zScale = pixelDepth/minScale;
		
		int rPixels=0;
		double maxLength=0;
		double[] count=new double[maxR];
		switch(valChoice)
		{
		case "Map 0":
			for(rPixels = 0; rPixels < maxR; rPixels++)
			{
				for(int i = 0; i< nSamples; i++)
				{
					randomPointPair(imgWidth,imgHeight,imgDepth,xScale,yScale,zScale,rPixels, p1,p2);
					p1.val = getVoxelValue(image,imgWidth,imgHeight,imgDepth,p1);
					p2.val = getVoxelValue(image,imgWidth,imgHeight,imgDepth,p2);
					//bresenhamLine reassigns  p1 or p2, send copies
					p1a.x = p1.x; p1a.y = p1.y; p1a.z = p1.z;
					p2a.x = p2.x; p2a.y = p2.y; p2a.z = p2.z;
					Point3D[] theLine = bresenhamLine(image,imgWidth,imgHeight,imgDepth,p1a,p2a);
					double length = vectorLength(p1,p2,pixelWidth,pixelHeight,pixelDepth);
					if(length>maxLength) maxLength=length;

					boolean nonZeroFound = false;
					for(int j=0; j<theLine.length;j++)
					{
						if(theLine[j].val != 0)
						{
							nonZeroFound=true;
							break;
						}
					}
					//only zero values were found along the line
					if(nonZeroFound==false)
					{
						count[rPixels]++;
					}
				}
				if(showProgress) pBars.setValue("Length", rPixels);
			}
			break;
		case "Map !0":
			for(rPixels = 0; rPixels < maxR; rPixels++)
			{
				for(int i = 0; i< nSamples; i++)
				{
					randomPointPair(imgWidth,imgHeight,imgDepth,xScale,yScale,zScale,rPixels, p1,p2);
					p1.val = getVoxelValue(image,imgWidth,imgHeight,imgDepth,p1);
					p2.val = getVoxelValue(image,imgWidth,imgHeight,imgDepth,p2);
					//bresenhamLine reassigns  p1 or p2, send copies
					p1a.x = p1.x; p1a.y = p1.y; p1a.z = p1.z;
					p2a.x = p2.x; p2a.y = p2.y; p2a.z = p2.z;
					Point3D[] theLine = bresenhamLine(image,imgWidth,imgHeight,imgDepth,p1a,p2a);
					double length = vectorLength(p1,p2,pixelWidth,pixelHeight,pixelDepth);
					if(length>maxLength) maxLength=length;
					
					boolean zeroFound = false;
					for(int j=0; j<theLine.length;j++)
					{
						if(theLine[j].val==0)
						{
							zeroFound=true;
							break;
						}
					}
					if(zeroFound==false)
					{
						count[rPixels]++;
					}
				}
				if(showProgress) pBars.setValue("Length", rPixels);
			}
			break;
		}
		//double binWidth	= (maxLength)/(nBins-1);
		double binWidth	= (maxLength)/(maxR);
		double[] xAxis = new double[maxR];
		for(int n=0;n<maxR;n++) xAxis[n] = n*binWidth;
		
		//Show as Probability
		for(int i=0;i<maxR;i++) count[i]/=nSamples;
		if(showProgress) pBars.close();

		return new Histogram(xAxis,count);
	}

	//*********************************************************************************
	
	/**	Find the probability that a random line of length R in pixels will lie completely in the chosen image space.
	 * @param image 1D Reference to a binarized 2D short image
	 * @param imgWidth The image width in pixels
	 * @param imgHeight The height width in pixels
	 * @param imgDepth The depth in pixels, pass 1 for 2D images.
	 * @param nSamples The number of samples for each line length
	 * @param valChoice "Map 0" or "Map 255": Note that  the ImageJ 255 short value = -1
	 * @param showProgress Display progress Bar while running
	 * @return Reference to a histogram of length (in pixels) vs probability
	 */
	public Histogram linealPathDistribution(Object[] image,int imgWidth,int imgHeight, int imgDepth,
					int nSamples, String valChoice, boolean showProgress)
	{
		return linealPathDistribution(image, imgWidth, imgHeight, imgDepth, 1, 1, 1, nSamples, valChoice, showProgress);
	}

	//*********************************************************************************

	/**Experimental: Estimates the pore size distribution as a probability density using local measurements of the surface-to-volume ratio and computing the radius of an equivalent sphere.
	 * <br>i &period;e &period; the likelihood of finding a pore of radius R in the interval between R +/- dR in image units.
	 * <br>This code measures the surface and volume the selected component within 
	 * a randomly located cube of maximum length = (min image dimension)/2.
	 * the actual size of the cube is selected by the user. And should be larger than a typical pore.
	* The pore collection is then binned by radius.
	 * Anisotropic pixels and voxels are supported. 
	 * @param image 1D Reference to a binarized 2D short image
	 * @param imgWidth The image width in pixels
	 * @param imgHeight The height width in pixels
	 * @param imgDepth The depth width in pixels
	 * @param pixelWidth the size of the voxel in the x, width, column, direction (e &period;g &period; cm , inch etc.)
	 * @param pixelHeight the size of the voxel in the y, imgHeight, row, direction (e &period;g &period; cm , inch etc.)
	 * @param pixelDepth the size of the voxel in the z, width, slice, direction (e &period;g &period; cm , inch etc.)
	 * @param probeDim the edge dimension of the probe cube in image units, should be slightly larger than typical pore
	 * @param nBins The number of bins in the histogram
	 * @param nSamples The number of samples for each line length
	 * @param valChoice "Map 0" or "Map !0"
	 * @param showProgress Display progress Bar while running
	 * @return Reference to a histogram of estimated pore radius (in pixels) vs count
	 */
	public Histogram poreSizeDistribution(Object[] image,int imgWidth,int imgHeight, int imgDepth , double probeDim,
			double pixelWidth, double pixelHeight, double pixelDepth,
			int nBins, int nSamples, String valChoice, boolean showProgress)
	{

		ProgressBars pBars = null;;
		if(showProgress)
		{	
			pBars = new ProgressBars("Pore Size Distribution");
			pBars.add("Sampling", 400, 20, 0, nSamples);
			pBars.show();
		}
	
		//compute the box dimensions in pixels
		int xDimPixels = (int)(probeDim/pixelWidth);
		int yDimPixels = (int)(probeDim/pixelHeight);
		int zDimPixels = (int)(probeDim/pixelDepth);

		//compute the maximum box index in pixels
		int xMaxPixels = imgWidth - xDimPixels - 1; 
		int yMaxPixels = imgHeight - yDimPixels - 1; 
		int zMaxPixels = imgDepth - zDimPixels - 1; 

		double voxVal,testVal;
		double pixelUnitsVolume = pixelWidth*pixelHeight*pixelDepth;
		double boxUnitsVolume = pixelUnitsVolume*Math.pow(probeDim, 3.0);
		double xyArea = pixelWidth*pixelHeight;
		double xzArea = pixelWidth*pixelDepth;
		double yzArea = pixelHeight*pixelDepth;
		double pixelUnitsArea = pixelWidth*pixelHeight;
		double boxUnitsArea = pixelUnitsArea*Math.pow(probeDim, 2.0);
		double surface,volume,area,perimeter=0,maxR=0;
		
//		Random rand = new Random();
		Point3D boxTLF = new Point3D();
		Point3D imgPt = new Point3D();
		Point3D testPt = new Point3D();
		Point3D[] offsets = initializeOffsets();

		double radius1;
		ArrayList<Double> radius = new ArrayList<Double>();
		
		if(imgDepth>2) // Image must be 3D, 3 slices is minimal 3D
		{
			for(int sample=0;sample<nSamples;sample++)
			{
				//get a random point for the box top left front
				//avoiding the x,y,z = zero indices
				boxTLF.x = myRandom.nextInt(xMaxPixels);
				boxTLF.y = myRandom.nextInt(yMaxPixels);
				boxTLF.z = myRandom.nextInt(zMaxPixels);
				if(boxTLF.x==0) boxTLF.x++;
				if(boxTLF.y==0) boxTLF.y++;
				if(boxTLF.z==0) boxTLF.z++;
				//reset counters
				int zeroCnt = 0,topBtmCnt = 0,leftRightCnt = 0,fwdBackCnt = 0;
				for( imgPt.z=boxTLF.z; imgPt.z<(boxTLF.z + zDimPixels); imgPt.z++)
				{
					for(imgPt.y=boxTLF.y; imgPt.y<(boxTLF.y + yDimPixels); imgPt.y++)
					{
						for(imgPt.x=boxTLF.x; imgPt.x<(boxTLF.x + xDimPixels); imgPt.x++)
						{
							voxVal = getVoxelValue(image,imgWidth,imgHeight,imgDepth,imgPt);
							if(voxVal==0) zeroCnt++;
							else // it is a matrix pixel
							{//look around and count how many of its faces are touching the pore voxels
								for(int i= 0;i<2;i++) //top,btm
								{
									testPt.x = imgPt.x + offsets[i].x;
									testPt.y = imgPt.y + offsets[i].y;
									testPt.z = imgPt.z + offsets[i].z;
									testVal = getVoxelValue(image,imgWidth,imgHeight,imgDepth,testPt);
									if(testVal==0) topBtmCnt++;
								}							
								for(int i= 2;i<4;i++) //left,right
								{
									testPt.x = imgPt.x + offsets[i].x;
									testPt.y = imgPt.y + offsets[i].y;
									testPt.z = imgPt.z + offsets[i].z;
									testVal = getVoxelValue(image,imgWidth,imgHeight,imgDepth,testPt);
									if(testVal==0) leftRightCnt++;
								} 
								for(int i= 4;i<6;i++) //front, back
								{
									testPt.x = imgPt.x + offsets[i].x;
									testPt.y = imgPt.y + offsets[i].y;
									testPt.z = imgPt.z + offsets[i].z;
									testVal = getVoxelValue(image,imgWidth,imgHeight,imgDepth,testPt);
									if(testVal==0) fwdBackCnt++;
								}
							}
						}
					}
				}
				
				if(valChoice=="Map 0") volume = zeroCnt*pixelUnitsVolume;
				else volume = boxUnitsVolume - zeroCnt*pixelUnitsVolume;				
				surface = topBtmCnt*xzArea + leftRightCnt*yzArea + fwdBackCnt*xyArea;
				if(volume>0 && surface > 0)
				{
					radius1 = (4.488*volume/surface); 
					if(radius1>maxR) maxR=radius1;
					radius.add(radius1);
				}				
				if(showProgress) pBars.setValue("Sampling",sample);
			}
		}
		else
		{
			for(int sample=0;sample<nSamples;sample++)
			{
				//get a random point for the box top left front
				boxTLF.x = myRandom.nextInt(xMaxPixels);
				boxTLF.y = myRandom.nextInt(yMaxPixels);
				testPt.z = 0;
				//avoiding the image edges
				if(boxTLF.x==0) boxTLF.x++;
				if(boxTLF.y==0) boxTLF.y++;
				//reset counters
				int zeroCnt = 0;
				perimeter=0;
				for(imgPt.y=boxTLF.y; imgPt.y<(boxTLF.y + yDimPixels); imgPt.y++)
				{
					for(imgPt.x=boxTLF.x; imgPt.x<(boxTLF.x + xDimPixels); imgPt.x++)
					{
						voxVal = getVoxelValue(image,imgWidth,imgHeight,imgDepth,imgPt);
						if(voxVal==0) zeroCnt++;
						else
						{
							for(int i= 0;i<2;i++)//up,down
							{
								testPt.x = imgPt.x + offsets[i].x;
								testPt.y = imgPt.y + offsets[i].y;
								testVal = getVoxelValue(image,imgWidth,imgHeight,imgDepth,testPt);
								if(testVal==0) perimeter += pixelWidth;
							}							
							for(int i= 2;i<4;i++) //left,right
							{
								testPt.x = imgPt.x + offsets[i].x;
								testPt.y = imgPt.y + offsets[i].y;
								testVal = getVoxelValue(image,imgWidth,imgHeight,imgDepth,testPt);
								if(testVal==0) perimeter += pixelHeight;
							} 
						}
					}
				}
				if(valChoice=="Map 0") area = zeroCnt*pixelUnitsArea;
				else area = boxUnitsArea - zeroCnt*pixelUnitsArea;
				if(area>0 && perimeter > 0)
				{
					radius1 = (2.5196*area/perimeter); 
					if(radius1>maxR) maxR=radius1;
					radius.add(radius1);
				}				
				if(showProgress) pBars.setValue("Sampling",sample);
			}
		}

		//Bin the Radii
		double xAxis[] = new double[nBins];
		double count[] = new double[nBins];		
		double binWidth	= (maxR)/(nBins-1);
		for(int n=0;n<nBins;n++) xAxis[n] = n*binWidth;
		int index;
		for(int n=0;n<radius.size();n++)
		{
			index	= (int) (radius.get(n)/binWidth);
			count[index]++;
		}
		//normalize
		for(int n=0;n<count.length;n++) count[n]/=nSamples;

		if(showProgress) pBars.close();
			
		return new Histogram(xAxis,count);
	}

	//*********************************************************************************

	/**Function to return a random point p1 within a 2D or 3D space and 
	 * a second random point p2 separated by sepDist in image units 
	 * <br>The scale factor for the smallest pixel dimension should be 1 and the others are relative to it&period;
	 * <br>e &period;g &period; pixelWidth = 10um, pixelHeight=20um and pixelDepth = 30um, xScale = 1, yScale=2, zScale=3
	 * <br>The pixel units must be the same.
	 * <br> Use getVoxelValue to obtain the Point3D.val at the Point3D coordinates in an image.  
	 * @param imgWidth The width of the space in pixels
	 * @param imgHeight The height of the space in pixels
	 * @param imgDepth The depth of the space in pixels. Pass 1 for 2D spaces
	 * @param sepDist The desired distance of p2 from p1
	 * @param xScale the relative length in the x direction
	 * @param yScale the relative length in the y direction
	 * @param zScale the relative length in the z direction
	 * @param p1 Returns the first random pixel's coordinates
	 * @param p2 Returns the second random pixel's coordinates
	 */
	public void randomPointPair(int imgWidth, int imgHeight, int imgDepth, double xScale, double yScale, double zScale, double sepDist, Point3D p1, Point3D p2)
	{
		double	sinTheta,sinPhi,cosTheta,cosPhi,l;

		if(imgDepth==1)
		{
			//get a point anywhere in the image
			p1.x = (int) (myRandom.nextDouble()*imgWidth);
			p1.y = (int) (myRandom.nextDouble()*imgHeight);
			p1.z = 0;
			//get random sinTheta and cosTheta
			double theta = 2.0*Math.PI*myRandom.nextDouble();
			sinTheta	= Math.sin(theta);
			cosTheta	= Math.cos(theta);		
			
			//Get second point at sepDist
			p2.x = (int)(p1.x + sepDist/xScale * cosTheta);
			p2.y = (int)(p1.y + sepDist/yScale * sinTheta);
			
			//Change Sense if out of bounds
			if(p2.x <0)
			{
				p2.x = 2*p1.x - p2.x;
			}
			else if(p2.x >=imgWidth)
			{
				l = (imgWidth - p1.x) + (p2.x - imgWidth);
				p2.x = (int)(p1.x - l);		    	
			}

			if(p2.y <0)
			{
				p2.y = 2*p1.y - p2.y;
			}
			else if(p2.y >=imgHeight)
			{
				l = (imgHeight - p1.y) + (p2.y - imgHeight);
				p2.y = (int)(p1.y - l);		    	
			}

		}
		else //3D
		{
			//get a point anywhere in the image
			p1.x = (int) (myRandom.nextDouble()*imgWidth);
			p1.y = (int) (myRandom.nextDouble()*imgHeight);
			p1.z = (int) (myRandom.nextDouble()*imgDepth);

			//Calling sin and cos on theta an phi is not significantly slower
			//than using random directly and sorting out the quadrants.
			double theta = 2.0*Math.PI*myRandom.nextDouble();
			sinTheta	= Math.sin(theta);
			cosTheta	= Math.cos(theta);		
			double phi 	= 2.0*Math.PI*myRandom.nextDouble();
			sinPhi		= Math.sin(phi);
			cosPhi		= Math.cos(phi);

			//Get second point at scaled sepDist
			p2.x = (int) Math.round(p1.x + sepDist/xScale * sinPhi * cosTheta);
			p2.y = (int) Math.round(p1.y + sepDist/yScale * sinPhi * sinTheta);
			p2.z = (int) Math.round(p1.z + sepDist/zScale * cosPhi);
						
			//Change Sense if out of bounds
			if(p2.x <0)
			{
				p2.x = 2*p1.x - p2.x;
			}
			else if(p2.x >=imgWidth)
			{
				l = (imgWidth - p1.x) + (p2.x - imgWidth);
				p2.x = (int)(p1.x - l);		    	
			}

			if(p2.y <0)
			{
				p2.y = 2*p1.y - p2.y;
			}
			else if(p2.y >=imgHeight)
			{
				l = (imgHeight - p1.y) + (p2.y - imgHeight);
				p2.y = (int)(p1.y - l);		    	
			}

			if(p2.z <0)
			{
				p2.z = 2*p1.z - p2.z;
			}
			else if(p2.z >=imgDepth)
			{
				l = (imgHeight - p1.z) + (p2.z - imgHeight);
				p2.z = (int)(p1.z - l);		    	
			}
		}
	}

	//*********************************************************************************

	/**Function to return a random point p1 within a 2D or 3D space and
	 * a second random point p2 separated by sepDist in pixels  
	 * @param imgWidth The image imgWidth in pixels
	 * @param imgHeight The height imgWidth in pixels
	 * @param imgDepth The depth in pixels
	 * @param sepDist The length of the line between p1 and p2
	 * @param p1 Returns the first random pixel coordinates
	 * @param p2 Returns the second random pixel coordinates
	 */
	public void randomPointPair(int imgWidth, int imgHeight, int imgDepth, double sepDist, Point3D p1, Point3D p2)
	{
		randomPointPair(imgWidth,imgHeight,imgDepth,1,1,1,sepDist,p1,p2);
	}
	
	//*********************************************************************************

	/**Sets the coordinates and voxel values along a Bresenham line between p1 and p2 in an image
	 * @param image the image data as an array of 2D slices. Pass z=0 for a 2D image.
	 * @param imgWidth The width of the image in pixels
	 * @param imgHeight The height of the image in pixels
	 * @param imgDepth The depth of the image in pixels
	 * @param p1 A Point3D descriptor of the x,y,z point start of the line
	 * @param p2 A Point3D descriptor of the x,y,z point end of the line
	 * @param val the new value for the pixels between P1 and p2
	 */	
	public void setLine(Object[] image, int imgWidth, int imgHeight, int imgDepth, Point3D p1, Point3D p2, double val)
	{
		Point3D[] line = bresenhamLine(p1,p2);

		if(image[0] instanceof byte[])
		{
			for(Point3D pt : line)
			{
				byte[] slice = (byte[])image[pt.z];
				slice[pt.x + pt.y*imgWidth] = (byte) val;
			}
		}
		else if(image[0] instanceof short[])
		{
			for(Point3D pt : line)
			{
				short[] slice = (short[])image[pt.z];
				slice[pt.x + pt.y*imgWidth] = (short) val;
			}
		}
		else if(image[0] instanceof int[])
		{
			for(Point3D pt : line)
			{
				int[] slice = (int[])image[pt.z];
				slice[pt.x + pt.y*imgWidth] = (int) val;
			}
		}
		else if(image[0] instanceof float[])
		{
			for(Point3D pt : line)
			{
				float[] slice = (float[])image[pt.z];
				slice[pt.x + pt.y*imgWidth] = (float) val;
			}
		}
		else if(image[0] instanceof double[])
		{
			for(Point3D pt : line)
			{
				double[] slice = (double[])image[pt.z];
				slice[pt.x + pt.y*imgWidth] = val;
			}
		}
	}

	//*********************************************************************************

	/**	Find the probability that two points separated by a distance D will both lie in the selected image space.
	 * Note 1. Distances are scanned in integer pixel steps.
	 * Note 2. The longest distance sampled is half the smallest image dimension.
	 * @param image  A reference to an 2D or 3D image segmented into zero and non-zero value regions
	 * @param imgWidth The image width in pixels
	 * @param imgHeight The height in pixels
	 * @param imgDepth The depth in pixels, pass 1 for 2D images.
	 * @param pixelWidth the size of the voxel in the x, width, column, direction (e &period;g &period; cm , inch etc.)
	 * @param pixelHeight the size of the voxel in the y, imgHeight, row, direction (e &period;g &period; cm , inch etc.)
	 * @param pixelDepth the size of the voxel in the z, width, slice, direction (e &period;g &period; cm , inch etc.), pass 1 for 2D images
	 * @param nSamples The number of samples at each distance step. More samples = less noise
	 * @param valChoice "Map 0" or "Map !0".
	 * @param showProgress Display progress Bar while running
	 * @return Reference to a Histogram of probability vs length in image units.
	 */
	public Histogram twoPointFunction(Object[] image,int imgWidth,int imgHeight, int imgDepth,
			double pixelWidth, double pixelHeight, double pixelDepth,
			int nSamples, String valChoice, boolean showProgress)
	{
		ProgressBars pBars = null;;
		int maxR = getMaxPixelSeparationDistance(imgWidth,imgHeight,imgDepth);
		if(showProgress)
		{	
			pBars = new ProgressBars("Two Point Distribution");
			pBars.add("Length", 400, 20, 0, maxR);
			pBars.show();
		}
		Point3D	p1 = new Point3D();
		Point3D	p2 = new Point3D();

		double minScale = pixelWidth;
		if(minScale> pixelHeight) minScale=pixelHeight;
		if(minScale> pixelDepth) minScale=pixelDepth;
		
		double xScale = pixelWidth/minScale;
		double yScale = pixelHeight/minScale;
		double zScale = pixelDepth/minScale;
		
		int rPixels=0;
		double maxLength=0;
		double[] count=new double[maxR];
		switch(valChoice)
		{
		case "Map 0":
			for(rPixels = 0; rPixels < maxR; rPixels++)
			{
				for(int i = 0; i< nSamples; i++)
				{
					randomPointPair(imgWidth,imgHeight,imgDepth,xScale,yScale,zScale,rPixels, p1,p2);
					p1.val = getVoxelValue(image,imgWidth,imgHeight,imgDepth,p1);
					p2.val = getVoxelValue(image,imgWidth,imgHeight,imgDepth,p2);
					if(p1.val==0 && p2.val==0)
					{
						double length = vectorLength(p1,p2,pixelWidth,pixelHeight,pixelDepth);
						if(length>maxLength) maxLength=length;
						count[rPixels]++;
					}
				}
				if(showProgress) pBars.setValue("Length", rPixels);
			}
			break;
		case "Map !0":
			for(rPixels = 0; rPixels < maxR; rPixels++)
			{
				for(int i = 0; i< nSamples; i++)
				{
					randomPointPair(imgWidth,imgHeight,imgDepth,xScale,yScale,zScale,rPixels, p1,p2);
					p1.val = getVoxelValue(image,imgWidth,imgHeight,imgDepth,p1);
					p2.val = getVoxelValue(image,imgWidth,imgHeight,imgDepth,p2);
					if(p1.val!=0 && p2.val!=0)
					{
						double length = vectorLength(p1,p2,pixelWidth,pixelHeight,pixelDepth);
						if(length>maxLength) maxLength=length;
						count[rPixels]++;
					}
				}
				if(showProgress) pBars.setValue("Length", rPixels);
			}
			break;
		}
		//double binWidth	= (maxLength)/(nBins-1);
		double binWidth	= (maxLength)/(maxR);
		double[] xAxis = new double[maxR];
		for(int n=0;n<maxR;n++) xAxis[n] = n*binWidth;
		
		//Show as Probability
		for(int i=0;i<maxR;i++) count[i]/=nSamples;
		if(showProgress) pBars.close();

		return new Histogram(xAxis,count);
	}

	//*********************************************************************************

	/**Get the length of the vector p1-p2 in image units
	 * @param p1 A Point3D descriptor of the x,y,z point start of the line
	 * @param p2 A Point3D descriptor of the x,y,z point end of the line
	 * @param pixelWidth the size of the voxel in the x, width, column, direction (e &period;g &period; cm , inch etc.)
	 * @param pixelHeight the size of the voxel in the y, imgHeight, row, direction (e &period;g &period; cm , inch etc.)
	 * @param pixelDepth the size of the voxel in the z, width, slice, direction (e &period;g &period; cm , inch etc.), ignored for 2D images
	 * @return The length of the vector p1-p2 in image Units
	 */
	public double vectorLength(Point3D p1, Point3D p2, double pixelWidth, double pixelHeight, double pixelDepth)
	{		
		double dx = (p1.x - p2.x) * pixelWidth;
		double dy = (p1.y - p2.y) * pixelHeight;
		double dz = (p1.z - p2.z) * pixelDepth;
		return Math.sqrt(dx*dx + dy*dy + dz*dz);
	}

	//*********************************************************************************

	/**Returns an array containing the length of chords along a line through the selected component.
	 * @param theLine An array of  Point3D descriptors of the line
	 * @param valChoice "Map !0" get chords along the non-zero portions of line, "Map 0" get from zero component
	 * @param pixelWidth the size of the voxel in the x, width, column, direction (e &period;g &period; cm , inch etc.)
	 * @param pixelHeight the size of the voxel in the y, imgHeight, row, direction (e &period;g &period; cm , inch etc.)
	 * @param pixelDepth the size of the voxel in the z, width, slice, direction (e &period;g &period; cm , inch etc.), ignored for 2D images
	 * @return an array of chord length along the line.
	 */
	public float[] getChords(Point3D[] theLine, String valChoice, double pixelWidth, double pixelHeight, double pixelDepth)
	{
		//there cannot be more chords than points
		float[] chords = new float[theLine.length];
		int chordCnt=0;

		//there cannot be more edges than points
		Point3D[] edges = new Point3D[theLine.length];
		int edgeCnt=0;
		int i,startEdgeIndex;

		//binarize the line
		for(i=0; i< theLine.length;i++)
		{
			if(theLine[i].val != 0) theLine[i].val= 100;			
		}
		
		//Find the location of cells where the line changes value
		for(i=1; i< theLine.length;i++)
		{
			try
			{
				if(theLine[i].val - theLine[i - 1].val != 0)
				{
					edges[edgeCnt] = theLine[i];
					edgeCnt++;
				}
			}
			catch(Exception e)
			{
				System.out.println(e.toString());
			}
		}
		//truncate the edges array
		Point3D[] edges2 = Arrays.copyOf(edges, edgeCnt);


		double di,dj,dk;
		switch(valChoice)
		{
		case "Map 0":
			//If the first cell is 0, skip to the second edge.
			//The first edge steps out of the chordVal space,
			//the second edge steps into it.
			if(theLine[0].val  == 0)  startEdgeIndex = 1; else startEdgeIndex = 0;
			//Find the chord lengths
			chordCnt = 0;
			for( i = startEdgeIndex;i< edges2.length - 1;i+=2)
			{
				di = (edges2[i + 1].x - edges2[i].x) * pixelWidth;
				dj = (edges2[i + 1].y - edges2[i].y) * pixelHeight;
				dk = (edges2[i + 1].z - edges2[i].z) * pixelDepth;
				chords[chordCnt] = (float)Math.sqrt(di*di + dj*dj + dk*dk);
				chordCnt ++;
			}			
			break;

		case "Map !0":
			//If the first cell is !0, skip to the second edge.
			//The first edge steps out of the chordVal space,
			//the second edge steps into it.
			if(theLine[0].val  != 0)  startEdgeIndex = 1; else startEdgeIndex = 0;
			//Find the chord lengths
			chordCnt = 0;
			for( i = startEdgeIndex;i< edges2.length - 1;i+=2)
			{
				di = (edges2[i + 1].x - edges2[i].x) * pixelWidth;
				dj = (edges2[i + 1].y - edges2[i].y) * pixelHeight;
				dk = (edges2[i + 1].z - edges2[i].z) * pixelDepth;
				chords[chordCnt] = (float)Math.sqrt(di*di + dj*dj + dk*dk);
				chordCnt ++;
			}			
			break;
		}

		return Arrays.copyOf(chords, chordCnt);
	}

	//*********************************************************************************

	private double getVoxelValue(Object[] image, int imgWidth, int imgHeight, int imgDepth, Point3D pt)
	{
		if (pt.x<0 && pt.x>imgWidth && pt.y<0 && pt.y>imgHeight && pt.z<0 && pt.z>imgDepth)
		{
			throw new IndexOutOfBoundsException();			
		}
		double value=Double.NaN;

		if(image[0] instanceof byte[])
		{
			byte[] slice = (byte[])image[pt.z];
			value = (double)slice[pt.x + pt.y*imgWidth];
		}
		else if(image[0] instanceof short[])
		{
			short[] slice = (short[])image[pt.z];
			value = slice[pt.x + pt.y*imgWidth];
		}
		else if(image[0] instanceof int[])
		{
			int[] slice = (int[])image[pt.z];
			value = (double)slice[pt.x + pt.y*imgWidth];
		}
		else if(image[0] instanceof float[])
		{
			float[] slice = (float[])image[pt.z];
			value = (double)slice[pt.x + pt.y*imgWidth];
		}
		else if(image[0] instanceof double[])
		{
			double[] slice = (double[])image[pt.z];
			value = (double)slice[pt.x + pt.y*imgWidth];
		}
		return value;
	}

	//*********************************************************************************

	private Point3D[] initializeOffsets()
	{// Creates an array of offsets to the 26 voxels touching the home voxel

		// offsetList is a list of offsets to i,j,k that give the location of nearest
		// neighbors ranging from 6 through 18 to 26 connected.
		// to vary the connectivity we simply go deeper into the list.
		//Point3D[] offset = new Point3D[26];
		Point3D[] offset = new Point3D[6];

		// left = -1;		up = -1;		fwd = -1
		// right = +1		down = +1		back = +1

		// initialize the offsets list  6 face touching neighbors
		offset[0] = new Point3D(0,-1,0);	//0, up, 0
		offset[1] = new Point3D(0,1,0);		//0, down, 0
		offset[2] = new Point3D(-1,0,0);	//left, 0, 0
		offset[3] = new Point3D(1,0,0);		//right, 0, 0
		offset[4] = new Point3D(0,0,-1);	//0, 0, fwd
		offset[5] = new Point3D(0,0,1);		//0, 0, back

		/*	// initialize the offsets list  12 edge touching neighbors
	offset[6] = new Point3D(-1,-1,0);	//left, up, 0
	offset[7] = new Point3D(1,-1,0);	//right, up, 0
	offset[8] = new Point3D(-1,1,0);	//left, down, 0
	offset[9] = new Point3D(1,1,0);		//right, down,0
	offset[10] = new Point3D(0,-1,-1);	//0, up, fwd
	offset[11] = new Point3D(0,-1,1);	//0, up, back
	offset[12] = new Point3D(0,1,-1);	//0, down, fwd
	offset[13] = new Point3D(0,1,1);	//0, down, back
	offset[14] = new Point3D(-1,0,-1);	//left, 0, fwd
	offset[15] = new Point3D(1,0,-1);	//right, 0, fwd
	offset[16] = new Point3D(-1,0,1);	//left, 0, back
	offset[17] = new Point3D(1,0,1);	//right, 0, back

	// initialize the offsets list  8 corner touching neighbors
	offset[18] = new Point3D(-1,-1,-1);	//left up fwd
	offset[19] = new Point3D(1,-1,-1);	//right up fwd 
	offset[20] = new Point3D(-1,1,-1);	//left down fwd 
	offset[21] = new Point3D(1,1,-1);	//right down fwd 
	offset[22] = new Point3D(-1,-1,1);	//left up back 
	offset[23] = new Point3D(1,-1,1);	//right up back 
	offset[24] = new Point3D(-1,1,1);	//left down back 
	offset[25] = new Point3D(1,1,1);	//right down back 
		 */
		return offset;
	}
}
	