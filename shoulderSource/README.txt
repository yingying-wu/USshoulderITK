README

* main cpp file: usshoulder_codeForTesting.cpp
* data file: shoulder_4_2_2014_86.png 

* project:
  - 3D reconstruction of ultrasound (US) images

* aim: 
  - to equalize the contrast between bright and dark region of image taken by a slightly malfunctioning US probe 
  - If we get a better probe, we won't need to do this, except for the first cropping to get the scan region

* approach:
  - The input image is the screen capture on the US imaging machine, so the first step is to crop the image to the region of the scan i.e. ignoring patient labels and other information about the scan that is on the machine screen
  - The bright and dark segments are cropped out to another image of the same size as the region of the scan because they will be later added together. The cropped image are padded with pizel of value 0 so they can be directly added without affecting each other. 
  - The bright and dark segments are processed separately with the SigmoidImageFilter to get approximately equal contrast (judged by eye). SigmoidImageFilter takes in a beta value (center of sigmoid curve) and an alpha value (spread of curve)
  - The two segments, after contrast treatment, are added together using the AddImageFilter
  - The merged image is then windowed to give pixel values between 0 - 255
  - We then blur the boundary between the two segments using a neighborhood iterator and a GaussianOperator along a thin rectangle fitted about the boundaries

* cropping:
  - boundaries are manually located. the top left corner and the bottom right corner are noted. From experience, even input from the same machine can have changing boundaries...
  - the RegionOfInterestImageFilter is used instead of the ExtractImageFilter because the ExtractImageFilter does not move the origin of the image and thus doesn't allow images to be added easily

* This code has a combination of SimpleITK, ITK 
  - SimpleITK is used to conveniently carry out simple procedures like showing images in imageJ where you can create histogram of the pixel values of a specific region
  - IJMacros was left in from using an external polar coordinates library for prostate image processing. Didn't take it out because the last time I did, cmake complained and I didn't have time to resolve that

* Running the code:
  - image file: hard coded
  - image boundaries: hard coded
  - input: usshoulder.exe beta_bright alpha_bright beta_dark alpha_dark region("bright", "dark", "all")
  - use "region" to tell program which image to show. bright/dark: show bright/dark region only after contrast. all: show combined final image after windowing

* Evaluation
  - There is still a dark boundary between the dark and bright regions, but because we are reconstructing the volume by taking maximum, the dark region will hopefully be filled in by other image slices

* Other files in this folder
  - four different images were tested and the parameters used for the images are recorded at the beginning of the cpp file
