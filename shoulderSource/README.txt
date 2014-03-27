Explanation of pipeline: 

* main cpp file: prostateTRUSwithSnake.cpp
* data file: Series181.png

* aim: to extract contour of prostate from ultrasound images
  - initial plan was to extract partial contours (i.e. open paths)
      - the code is still in "prostateTRUS.cpp" source code, but to run this, you need to change the CMakeLists file to include this source code instead
      - the parameter list for command line 
          - prostateTRUS.exe (filepath) 12 0.08 2 0.04 0.04 0.08 1 5
      - documentation is in the word document "parameter_testing.docx" from the start up to and including "partial conclusion"
  - but it seems like the blurring improves the full contour, so I'm actually looking at active contour (i.e. closed contours) now

* This code has a combination of SimpleITK, ITK and an external polar transform commands
  - SimpleITK is used to conveniently carry out simple procedures like reading the input image, writing the output image and showing images in imageJ.
  - IJMacros, itkCartesianToPolarTransform.cpp, itkCartesianToPolarTransform.txx, itkPolarToCartesianTransform.cpp, itkPolarToCartesianTransform.txx are all for implementing polar coordinates.

* The initial part is to blur the image by converting the cartesian coordinates to polar coordinates:
  - I take the origin at either the probe center or the prostate center to smooth along different directions.
  - For polar transform to sample correctly, origin has to be first shifted using ChangeInformationImageFilter.
  - Overall blur is performed with MedianImageFilter to remove speckle and noise. 
  - Subsequent blurring is performed with RecursiveGaussianImageFilter.
  - The sequence of blurring is as follows:
      - Centered at probe, blur along theta,
      - Centered at prostate, blur along theta,
      - Centered at probe, blur along radius. 

* Then I apply a sequence of filters to extract the contour
  - CurvatureFlowImageFilter to smooth out jaggy edges
  - GradientMagnitudeRecursiveGaussianImageFilter to obtain feature input
  - SigmoidImageFilter
  - FastMarchingImageFilter with three seed points
      - (350, 400), (310, 500) and (300, 600)
  - IntensityWindowImageFilter to get rid of extreme points
  - BinaryThresholdImageFilter

* Running the code:
  - This code currently shows the image and image information for
    - all iamges for blurring,
    - speed image for fast marching,
    - fast marching results
    - binary thresholding results.
  - If you want to see any other intermmediate images, just uncomment the short code following every filter/procedure. It will display the image in ImageJ. 
  - Type this in command line to run "prostateTRUS.exe (filepath) 12 0.04 0.04 0.08 100 0.02 0 1 0.6 220".
  - Current version does not write file, but only displays images in ImageJ. If you want to save the files, uncomment the last bit at the end about writing.

* Decision making
  - Lots of tests were carried out to search for the best filter sequence and parameters. They are documented in the word document named "parameter_testing.docx" in this folder.

* Evaluation
  - It's not very accurate because the obvious contours are not closed and there is "leakage" of the fast marching map
  - It is also not generic enough to be applied to all images without redoing the parameter tuning
  - But the blurring allows better visualization of what is the actual contour because it gets rid of the artefacts and noise

* Other files in this folder
  - ppt for class presentation.
  - prostateTRUS.cpp is an older version that searches for partial contours without involving active contours.
