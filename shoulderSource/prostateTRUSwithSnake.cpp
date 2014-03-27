// Spring 2013 Medical Image Analysis
// Project name: Extracting Partial Contours of Prostate on 2D Ultrasound Images 
// Created by: Ying Ying Wu

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "SimpleITK.h"

#include "itkVector.h"

#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"

#include "itkCartesianToPolarTransform.h"
#include "itkPolarToCartesianTransform.h"
#include "vnl/vnl_math.h"

#include "itkMedianImageFilter.h"
#include "itkChangeInformationImageFilter.h"
#include "itkResampleImageFilter.h"
#include "itkNearestNeighborInterpolateImageFunction.h"
#include "itkRecursiveGaussianImageFilter.h"

#include "itkCurvatureFlowImageFilter.h"
#include "itkGradientMagnitudeRecursiveGaussianImageFilter.h"
#include "itkStatisticsImageFilter.h"
#include "itkSigmoidImageFilter.h"
#include "itkFastMarchingImageFilter.h"
#include "itkIntensityWindowingImageFilter.h"
#include "itkBinaryThresholdImageFilter.h"

// define mathematical constant
const double PI = 3.14159265359;

// using SimpleITK for simple procedures like reading in images and showing images in imageJ
namespace sitk = itk::simple;

// define a few global typedefs that are repeatedly used to define images and classes
typedef itk::Image< float, 2 > InternalITKImageType;

// make theta resolution in polar coordinates half a degree
const double thetaResolution = PI/360;


int main( int argc, char * argv[] )
{
  // make sure enough arguments are passed into exe
  if( argc < 12 )
    {
    std::cerr << "Usage: " << std::endl;
    std::cerr << argv[0] << "prostateTRUS.exe \
							filepath \
							radius(MedianImageFilter) \
							sigma(RecursiveGaussianImageFilter) \
							sigma(RecursiveGaussianImageFilter) \
							sigma(RecursiveGaussianImageFilter) \
							numberOfIterations(CurvatureFlowImageFilter) \
							sigma(GradientMagnitudeRecursiveGaussianImageFilter) \
							beta(SigmoidImageFilter) \
							alpha(SigmoidImageFilter) \
							normalizationFactor(FastMarchingImageFilter) \
							upperThreshold(BinaryThresholdImageFilter)"<< std::endl;
    return EXIT_FAILURE;
    }

  // declare typedefs here
  typedef itk::MedianImageFilter< InternalITKImageType, InternalITKImageType >  MedianFilterType;
  typedef itk::ChangeInformationImageFilter<InternalITKImageType> ChangeOriginFilterType;
  typedef itk::ResampleImageFilter<InternalITKImageType,InternalITKImageType> ResampleFilterType;
  typedef itk::PolarToCartesianTransform< double, 2 > PolarTransformType;
  typedef itk::NearestNeighborInterpolateImageFunction< InternalITKImageType, double > InterpolatorType;
  typedef itk::RecursiveGaussianImageFilter< InternalITKImageType, InternalITKImageType > GaussianFilterType;
  typedef itk::CartesianToPolarTransform< double, 2 > CartesianTransformType; 
  typedef itk::BinaryThresholdImageFilter< InternalITKImageType, InternalITKImageType > ThresholdingFilterType;
  typedef itk::GradientMagnitudeRecursiveGaussianImageFilter< InternalITKImageType, InternalITKImageType > GradientFilterType;
  typedef itk::StatisticsImageFilter< InternalITKImageType > StatisticsFilterType;
  typedef itk::SigmoidImageFilter< InternalITKImageType, InternalITKImageType > SigmoidFilterType;
  typedef itk::FastMarchingImageFilter< InternalITKImageType, InternalITKImageType >FastMarchingFilterType;
  typedef FastMarchingFilterType::NodeContainer NodeContainer;
  typedef FastMarchingFilterType::NodeType NodeType;
  typedef itk::IntensityWindowingImageFilter< InternalITKImageType, InternalITKImageType >  IntensityFilterType;
  typedef itk::CurvatureFlowImageFilter< InternalITKImageType, InternalITKImageType > FilterType;

  //
  // ImageFileReader
  //
  sitk::ImageFileReader reader;
  reader.SetFileName( argv[1] );
  sitk::Image sitkImageIn = reader.Execute();
  if ( sitkImageIn.GetDimension() != 2 )
  {
	std::cerr << "Image dimensions must match!"<<std::endl;
	return EXIT_FAILURE; 
  }
  
  //
  // convert sitk image to itk image
  //
  sitk::CastImageFilter caster;
  caster.SetOutputPixelType( sitk::sitkFloat32 );
  sitkImageIn = caster.Execute( sitkImageIn );
  
  InternalITKImageType::Pointer itkImage;
  itkImage = dynamic_cast <InternalITKImageType*> ( sitkImageIn.GetITKBase() );
  InternalITKImageType::SpacingValueType inputImageSpacing = itkImage->GetSpacing()[0];
  InternalITKImageType::SizeValueType inputImageSizeX = itkImage->GetLargestPossibleRegion().GetSize()[0];
  InternalITKImageType::SizeValueType inputImageSizeY = itkImage->GetLargestPossibleRegion().GetSize()[1];
  
  /*
  printf("itkImage size %d %d\n", itkImage->GetLargestPossibleRegion().GetSize()[0], itkImage->GetLargestPossibleRegion().GetSize()[1] );
  printf("itkImage spacing %f %f\n", itkImage->GetSpacing()[0], itkImage->GetSpacing()[1] );
  printf("itkImage origin %f %f\n", itkImage->GetOrigin()[0], itkImage->GetOrigin()[1] );
    
  // Show image
  sitk::Show( sitk::Image( itkImage ) );
  */
  
  //---------------------
  // BLUR
  //---------------------
  
  //
  // MedianImageFilter 
  // blurs input image to remove speckles and artifacts
  //
  MedianFilterType::Pointer medianFilter = MedianFilterType::New();  
  medianFilter->SetInput( itkImage ); 
  medianFilter->SetRadius( atoi( argv[2] ) );
  medianFilter->Update();
  
  /*
  InternalITKImageType::Pointer medianImage = medianFilter->GetOutput();
  printf("medianImage size %d %d\n", medianImage->GetLargestPossibleRegion().GetSize()[0], medianImage->GetLargestPossibleRegion().GetSize()[1] );
  printf("medianImage spacing %f %f\n", medianImage->GetSpacing()[0], medianImage->GetSpacing()[1] );
  printf("medianImage origin %f %f\n", medianImage->GetOrigin()[0], medianImage->GetOrigin()[1] );
   
  // Show image
  sitk::Show( sitk::Image( medianFilter->GetOutput() ) );
  */

  //
  // ChangeInformationImageFilter 
  // shift origin to center of probe (right-most column, middle row) 
  //
  ChangeOriginFilterType::Pointer probeOriginFilter = ChangeOriginFilterType::New();
  
  // value where the output image origin should start i.e. start value of output image sampling
  InternalITKImageType::PointType probeOrigin;
  probeOrigin[0] = -( inputImageSizeX - inputImageSpacing );
  probeOrigin[1] = -( floor( inputImageSizeY / 2 ) - inputImageSpacing );
  probeOriginFilter->ChangeOriginOn();
  probeOriginFilter->SetOutputOrigin( probeOrigin );

  probeOriginFilter->SetInput( medianFilter->GetOutput() );
  probeOriginFilter->Update();
  
  /*
  InternalITKImageType::Pointer probeOriginImage = probeOriginFilter->GetOutput();
  printf("probeOriginImage size %d %d\n", probeOriginImage->GetLargestPossibleRegion().GetSize()[0], probeOriginImage->GetLargestPossibleRegion().GetSize()[1] );
  printf("probeOriginImage spacing %f %f\n", probeOriginImage->GetSpacing()[0], probeOriginImage->GetSpacing()[1] );
  printf("probeOriginImage origin %f %f\n", probeOriginImage->GetOrigin()[0], probeOriginImage->GetOrigin()[1] );
  
  // Show image
  sitk::Show( sitk::Image( probeOriginFilter->GetOutput() ) );
  */

  //
  // ResampleImageFilter - PolarToCartesianTransform
  // converts to polar coordinates, naming is based on output coordinates to input coordinates
  // everything in spatial coordinates i.e. mm
  //
  ResampleFilterType::Pointer polarResampleFilter = ResampleFilterType::New();

  PolarTransformType::Pointer polarTransform = PolarTransformType::New();
  polarResampleFilter->SetTransform( polarTransform );
  
  InterpolatorType::Pointer polarInterpolator = InterpolatorType::New();
  polarResampleFilter->SetInterpolator( polarInterpolator );

  polarResampleFilter->SetDefaultPixelValue( 0 );
  
  // make r the same spatial resolution as input image
  double polarSpacing[ 2 ] = { thetaResolution, inputImageSpacing }; // (theta, r) pixel spacing in millimeters
  polarResampleFilter->SetOutputSpacing( polarSpacing );
  
  // value where the output image origin should start i.e. start value of output image sampling
  // theta range from pi/2 to 3*pi/2 since origin is on rightmost column
  double polarOrigin[ 2 ] = { PI/2, 0.0 }; // (theta, r) space coordinate of origin
  polarResampleFilter->SetOutputOrigin( polarOrigin );

  InternalITKImageType::SizeType polarSize = { 360, 512 }; // (theta, r) number of pixels
  polarResampleFilter->SetSize( polarSize );
  
  polarResampleFilter->SetInput( probeOriginFilter->GetOutput() );
  polarResampleFilter->Update();
  
  /*
  InternalITKImageType::Pointer polarImage = polarResampleFilter->GetOutput();
  printf("polarImage size %d %d\n", polarImage->GetLargestPossibleRegion().GetSize()[0], polarImage->GetLargestPossibleRegion().GetSize()[1] );
  printf("polarImage spacing %f %f\n", polarImage->GetSpacing()[0], polarImage->GetSpacing()[1] );
  printf("polarImage origin %f %f\n", polarImage->GetOrigin()[0], polarImage->GetOrigin()[1] );
  
  // Show image
  sitk::Show( sitk::Image( polarResampleFilter->GetOutput() ) );
  */
  
  //
  // RecursiveGaussianImageFilter
  // smooth along theta to remove artifacts radiating out from probe
  //  
  GaussianFilterType::Pointer gaussianFilter = GaussianFilterType::New();

  gaussianFilter->SetSigma( atof(argv[3]) ); 
  gaussianFilter->SetDirection( 0 ); // (theta, r)
  gaussianFilter->SetNormalizeAcrossScale( true );

  gaussianFilter->SetInput( polarResampleFilter->GetOutput() );
  gaussianFilter->Update();
  
  /*
  InternalITKImageType::Pointer gaussianImage = gaussianFilter->GetOutput();
  printf("gaussianImage size %d %d\n", gaussianImage->GetLargestPossibleRegion().GetSize()[0], gaussianImage->GetLargestPossibleRegion().GetSize()[1] );
  printf("gaussianImage spacing %f %f\n", gaussianImage->GetSpacing()[0], gaussianImage->GetSpacing()[1] );
  printf("gaussianImage origin %f %f\n", gaussianImage->GetOrigin()[0], gaussianImage->GetOrigin()[1] );
  
  // Show image
  sitk::Show( sitk::Image( gaussianFilter->GetOutput() ) );
  */

  //
  // ResampleImageFilter - CartesianToPolarTransform
  //  
  ResampleFilterType::Pointer cartesianResampleFilter = ResampleFilterType::New();
 
  CartesianTransformType::Pointer cartesianTransformBlurred = CartesianTransformType::New();
  cartesianResampleFilter->SetTransform( cartesianTransformBlurred );
  
  InterpolatorType::Pointer cartesianInterpolatorBlurred = InterpolatorType::New();
  cartesianResampleFilter->SetInterpolator( cartesianInterpolatorBlurred );
  
  cartesianResampleFilter->SetDefaultPixelValue( 0 );
  
  // make it back to same spacing as input image
  double cartesianSpacing[ 2 ] = { inputImageSpacing, inputImageSpacing }; // (x, y) pixel spacing in millimeters
  cartesianResampleFilter->SetOutputSpacing( cartesianSpacing );
  
  cartesianResampleFilter->SetOutputOrigin( probeOrigin );
  
  InternalITKImageType::SizeType cartesianSize = { inputImageSizeX, inputImageSizeY }; // (x, y) number of pixels 
  cartesianResampleFilter->SetSize( cartesianSize );
  
  cartesianResampleFilter->SetInput( gaussianFilter->GetOutput() );
  cartesianResampleFilter->Update();
  
  /*
  InternalITKImageType::Pointer cartesianImage = cartesianResampleFilter->GetOutput();
  printf("cartesianImage size %d %d\n", cartesianImage->GetLargestPossibleRegion().GetSize()[0], cartesianImage->GetLargestPossibleRegion().GetSize()[1] );
  printf("cartesianImage spacing %f %f\n", cartesianImage->GetSpacing()[0], cartesianImage->GetSpacing()[1] );
  printf("cartesianImage origin %f %f\n", cartesianImage->GetOrigin()[0], cartesianImage->GetOrigin()[1] );
  
  // Show image
  sitk::Show( sitk::Image( cartesianResampleFilter->GetOutput() ) );
  */
  
  //
  // ChangeInformationImageFilter
  // shift to center of prostate (middle column, middle row)
  //  
  ChangeOriginFilterType::Pointer prostateOriginFilter = ChangeOriginFilterType::New();
  
  // value where the output image origin should start i.e. start value of output image sampling
  InternalITKImageType::PointType prostateOrigin;
  prostateOrigin[0] = -( floor( inputImageSizeX / 2 ) - inputImageSpacing ); 
  prostateOrigin[1] = -( floor( inputImageSizeY / 2 ) - inputImageSpacing );
  prostateOriginFilter->ChangeOriginOn();
  prostateOriginFilter->SetOutputOrigin( prostateOrigin );
  printf("origin %f %f\n", prostateOrigin[0], prostateOrigin[1] );

  prostateOriginFilter->SetInput( cartesianResampleFilter->GetOutput() );
  prostateOriginFilter->Update();  

  /*
  InternalITKImageType::Pointer prostateOriginImage = prostateOriginFilter->GetOutput();
  printf("prostateOriginImage size %d %d\n", prostateOriginImage->GetLargestPossibleRegion().GetSize()[0], prostateOriginImage->GetLargestPossibleRegion().GetSize()[1] );
  printf("prostateOriginImage spacing %f %f\n", prostateOriginImage->GetSpacing()[0], prostateOriginImage->GetSpacing()[1] );
  printf("prostateOriginImage origin %f %f\n", prostateOriginImage->GetOrigin()[0], prostateOriginImage->GetOrigin()[1] );
  
  // Show image
  sitk::Show( sitk::Image( prostateOriginFilter->GetOutput() ) );
  */
  
  //
  // ResampleImageFilter - PolarToCartesianTransform2
  // converts to polar coordinates, naming is based on output coordinates to input coordinates
  // everything in spatial coordinates i.e. mm
  //
  ResampleFilterType::Pointer polarResampleFilter2 = ResampleFilterType::New();
  
  PolarTransformType::Pointer polarTransform2 = PolarTransformType::New();
  polarResampleFilter2->SetTransform( polarTransform2 );
  
  InterpolatorType::Pointer polarInterpolator2 = InterpolatorType::New();
  polarResampleFilter2->SetInterpolator( polarInterpolator2 );
    
  polarResampleFilter2->SetDefaultPixelValue( 0 );
  polarResampleFilter2->SetOutputSpacing( polarSpacing );
  
  // value where the output image origin should start i.e. start value of output image sampling
  // theta range from 0 to 2*pi since origin is in the center
  double polarOrigin2[ 2 ] = { 0.0, 0.0 }; // (theta, r) space coordinate of origin
  polarResampleFilter2->SetOutputOrigin( polarOrigin2 );

  InternalITKImageType::SizeType polarSize2 = { 720, 512 }; // (theta, r) number of pixels
  polarResampleFilter2->SetSize( polarSize2 );
  
  polarResampleFilter2->SetInput( prostateOriginFilter->GetOutput() );
  polarResampleFilter2->Update();

  /*
  InternalITKImageType::Pointer polarImage2 = polarResampleFilter2->GetOutput();
  printf("polarImage2 size %d %d\n", polarImage2->GetLargestPossibleRegion().GetSize()[0], polarImage2->GetLargestPossibleRegion().GetSize()[1] );
  printf("polarImage2 spacing %f %f\n", polarImage2->GetSpacing()[0], polarImage2->GetSpacing()[1] );
  printf("polarImage2 origin %f %f\n", polarImage2->GetOrigin()[0], polarImage2->GetOrigin()[1] );
  
  // Show image
  sitk::Show( sitk::Image( polarResampleFilter2->GetOutput() ) );
  */

  //
  // RecursiveGaussianImageFilter
  // smooth along theta to enhance prostate boundary
  //  
  GaussianFilterType::Pointer gaussianFilter2 = GaussianFilterType::New();

  gaussianFilter2->SetSigma( atof(argv[4]) );
  gaussianFilter2->SetDirection( 0 ); //(theta, r)
  gaussianFilter2->SetNormalizeAcrossScale( true );

  gaussianFilter2->SetInput( polarResampleFilter2->GetOutput() );
  gaussianFilter2->Update();  

  /*
  InternalITKImageType::Pointer gaussianImage2 = gaussianFilter2->GetOutput();
  printf("gaussianImage2 size %d %d\n", gaussianImage2->GetLargestPossibleRegion().GetSize()[0], gaussianImage2->GetLargestPossibleRegion().GetSize()[1] );
  printf("gaussianImage2 spacing %f %f\n", gaussianImage2->GetSpacing()[0], gaussianImage2->GetSpacing()[1] );
  printf("gaussianImage2 origin %f %f\n", gaussianImage2->GetOrigin()[0], gaussianImage2->GetOrigin()[1] );
  
  // Show image
  sitk::Show( sitk::Image( gaussianFilter2->GetOutput() ) );
  */
  
  //
  // ResampleImageFilter - CartesianToPolarTransform 
  //  
  ResampleFilterType::Pointer cartesianResampleFilter2 = ResampleFilterType::New();
 
  CartesianTransformType::Pointer cartesianTransform2 = CartesianTransformType::New();
  cartesianResampleFilter2->SetTransform( cartesianTransform2 );
  
  InterpolatorType::Pointer cartesianInterpolator2 = InterpolatorType::New();
  cartesianResampleFilter2->SetInterpolator( cartesianInterpolator2 );
  
  cartesianResampleFilter2->SetDefaultPixelValue( 0 );
  cartesianResampleFilter2->SetOutputSpacing( cartesianSpacing );
  cartesianResampleFilter2->SetOutputOrigin( prostateOrigin );
  cartesianResampleFilter2->SetSize( cartesianSize );
  
  cartesianResampleFilter2->SetInput( gaussianFilter2->GetOutput() );
  cartesianResampleFilter2->Update();  

  /*
  InternalITKImageType::Pointer cartesianImage2 = cartesianResampleFilter2->GetOutput();
  printf("cartesianImage2 size %d %d\n", cartesianImage2->GetLargestPossibleRegion().GetSize()[0], cartesianImage2->GetLargestPossibleRegion().GetSize()[1] );
  printf("cartesianImage2 spacing %f %f\n", cartesianImage2->GetSpacing()[0], cartesianImage2->GetSpacing()[1] );
  printf("cartesianImage2 origin %f %f\n", cartesianImage2->GetOrigin()[0], cartesianImage2->GetOrigin()[1] );
  
  // Show image
  sitk::Show( sitk::Image( cartesianResampleFilter2->GetOutput() ) );
  */
  
  //
  // ChangeInformationImageFilter
  // shift to center of probe (right-most column, middle row)
  //
  ChangeOriginFilterType::Pointer probeOriginFilter2 = ChangeOriginFilterType::New();
  
  // value where the output image origin should start i.e. start value of output image sampling
  probeOriginFilter2->ChangeOriginOn();
  probeOriginFilter2->SetOutputOrigin( probeOrigin );

  probeOriginFilter2->SetInput( cartesianResampleFilter2->GetOutput() );
  probeOriginFilter2->Update();
  
  /*
  InternalITKImageType::Pointer probeOriginImage2 = probeOriginFilter2->GetOutput();
  printf("probeOriginImage2 size %d %d\n", probeOriginImage2->GetLargestPossibleRegion().GetSize()[0], probeOriginImage2->GetLargestPossibleRegion().GetSize()[1] );
  printf("probeOriginImage2 spacing %f %f\n", probeOriginImage2->GetSpacing()[0], probeOriginImage2->GetSpacing()[1] );
  printf("probeOriginImage2 origin %f %f\n", probeOriginImage2->GetOrigin()[0], probeOriginImage2->GetOrigin()[1] );
  
  // Show image
  sitk::Show( sitk::Image( probeOriginFilter2->GetOutput() ) );
  */

  //
  // ResampleImageFilter - PolarToCartesianTransform
  // converts to polar coordinates, naming is based on output coordinates to input coordinates
  // everything in spatial coordinates i.e. mm
  //
  ResampleFilterType::Pointer polarResampleFilter3 = ResampleFilterType::New();

  PolarTransformType::Pointer polarTransform3 = PolarTransformType::New();
  polarResampleFilter3->SetTransform( polarTransform3 );
  
  InterpolatorType::Pointer polarInterpolator3 = InterpolatorType::New();
  polarResampleFilter3->SetInterpolator( polarInterpolator3 );

  polarResampleFilter3->SetDefaultPixelValue( 0 );
  polarResampleFilter3->SetOutputSpacing( polarSpacing );
  polarResampleFilter3->SetOutputOrigin( polarOrigin );
  polarResampleFilter3->SetSize( polarSize );
  
  polarResampleFilter3->SetInput( probeOriginFilter2->GetOutput() );
  polarResampleFilter3->Update();
  
  /*
  InternalITKImageType::Pointer polarImage3 = polarResampleFilter3->GetOutput();
  printf("polarImage3 size %d %d\n", polarImage3->GetLargestPossibleRegion().GetSize()[0], polarImage3->GetLargestPossibleRegion().GetSize()[1] );
  printf("polarImage3 spacing %f %f\n", polarImage3->GetSpacing()[0], polarImage3->GetSpacing()[1] );
  printf("polarImage3 origin %f %f\n", polarImage3->GetOrigin()[0], polarImage3->GetOrigin()[1] );
  
  // Show image
  sitk::Show( sitk::Image( polarResampleFilter3->GetOutput() ) );
  */

  //
  // RecursiveGaussianImageFilter
  // smooth along radius to further enhance short curves
  //
  GaussianFilterType::Pointer gaussianFilter3 = GaussianFilterType::New();

  gaussianFilter3->SetSigma( atof(argv[5]) ); 
  gaussianFilter3->SetDirection( 1 ); // (theta, r)
  gaussianFilter3->SetNormalizeAcrossScale( true );

  gaussianFilter3->SetInput( polarResampleFilter3->GetOutput() );
  gaussianFilter3->Update();
  
  /*
  InternalITKImageType::Pointer gaussianImage3 = gaussianFilter3->GetOutput();
  printf("gaussianImage3 size %d %d\n", gaussianImage3->GetLargestPossibleRegion().GetSize()[0], gaussianImage3->GetLargestPossibleRegion().GetSize()[1] );
  printf("gaussianImage3 spacing %f %f\n", gaussianImage3->GetSpacing()[0], gaussianImage3->GetSpacing()[1] );
  printf("gaussianImage3 origin %f %f\n", gaussianImage3->GetOrigin()[0], gaussianImage3->GetOrigin()[1] );
  
  // Show image
  sitk::Show( sitk::Image( gaussianFilter3->GetOutput() ) );
  */

  //
  // ResampleImageFilter - CartesianToPolarTransform 
  //
  ResampleFilterType::Pointer cartesianResampleFilter3 = ResampleFilterType::New();
 
  CartesianTransformType::Pointer cartesianTransform3 = CartesianTransformType::New();
  cartesianResampleFilter3->SetTransform( cartesianTransform3 );
  
  InterpolatorType::Pointer cartesianInterpolator3 = InterpolatorType::New();
  cartesianResampleFilter3->SetInterpolator( cartesianInterpolator3 );
  
  cartesianResampleFilter3->SetDefaultPixelValue( 0 );
  cartesianResampleFilter3->SetOutputSpacing( cartesianSpacing );
  cartesianResampleFilter3->SetOutputOrigin( probeOrigin );
  cartesianResampleFilter3->SetSize( cartesianSize );
  
  cartesianResampleFilter3->SetInput( gaussianFilter3->GetOutput() );
  cartesianResampleFilter3->Update();
  
  /*
  InternalITKImageType::Pointer cartesianImage3 = cartesianResampleFilter3->GetOutput();
  printf("cartesianImage3 size %d %d\n", cartesianImage3->GetLargestPossibleRegion().GetSize()[0], cartesianImage3->GetLargestPossibleRegion().GetSize()[1] );
  printf("cartesianImage3 spacing %f %f\n", cartesianImage3->GetSpacing()[0], cartesianImage3->GetSpacing()[1] );
  printf("cartesianImage3 origin %f %f\n", cartesianImage3->GetOrigin()[0], cartesianImage3->GetOrigin()[1] );
  
  // Show image
  sitk::Show( sitk::Image( cartesianResampleFilter3->GetOutput() ) );
  */

  //
  // ChangeInformationImageFilter 
  // shift to lower left corner (0,0)
  //
  ChangeOriginFilterType::Pointer cartesianOriginFilter = ChangeOriginFilterType::New();
  
  // value where the output image origin should start i.e. start value of output image sampling
  InternalITKImageType::PointType finalImageOrigin;
  finalImageOrigin[0] = 0;
  finalImageOrigin[1] = 0;
  cartesianOriginFilter->ChangeOriginOn();
  cartesianOriginFilter->SetOutputOrigin( finalImageOrigin );
  printf("origin %f %f\n", finalImageOrigin[0], finalImageOrigin[1] );

  cartesianOriginFilter->SetInput( cartesianResampleFilter3->GetOutput() );
  cartesianOriginFilter->Update();
  
  /*
  InternalITKImageType::Pointer cartesianOriginImage = cartesianOriginFilter->GetOutput();
  printf("cartesianOriginImage size %d %d\n", cartesianOriginImage->GetLargestPossibleRegion().GetSize()[0], cartesianOriginImage->GetLargestPossibleRegion().GetSize()[1] );
  printf("cartesianOriginImage spacing %f %f\n", cartesianOriginImage->GetSpacing()[0], cartesianOriginImage->GetSpacing()[1] );
  printf("cartesianOriginImage origin %f %f\n", cartesianOriginImage->GetOrigin()[0], cartesianOriginImage->GetOrigin()[1] );
  
  // Show image
  sitk::Show( sitk::Image( cartesianOriginFilter->GetOutput() ) );
  */

  //---------------------
  // SNAKE
  //---------------------
  
  //
  // CurvatureFlowImageFilter
  //
  FilterType::Pointer curvatureFlowFilter = FilterType::New();

  curvatureFlowFilter->SetNumberOfIterations( atof( argv[6] ) );
  curvatureFlowFilter->SetTimeStep( 0.125 );

  curvatureFlowFilter->SetInput( cartesianOriginFilter->GetOutput() );
  curvatureFlowFilter->Update();


  InternalITKImageType::Pointer curvatureFlowImage = curvatureFlowFilter->GetOutput();
  printf("curvatureFlowImage size %d %d\n", curvatureFlowImage->GetLargestPossibleRegion().GetSize()[0], curvatureFlowImage->GetLargestPossibleRegion().GetSize()[1] );
  printf("curvatureFlowImage spacing %f %f\n", curvatureFlowImage->GetSpacing()[0], curvatureFlowImage->GetSpacing()[1] );
  printf("curvatureFlowImage origin %f %f\n", curvatureFlowImage->GetOrigin()[0], curvatureFlowImage->GetOrigin()[1] );
  
  // Show image
  sitk::Show( sitk::Image( curvatureFlowFilter->GetOutput() ) );


  // 
  // GradientMagniudeRecursiveGaussianImageFilter
  //
  GradientFilterType::Pointer gradientMagnitude = GradientFilterType::New();
  
  gradientMagnitude->SetSigma( atof( argv[7] ) );
  
  gradientMagnitude->SetInput( curvatureFlowFilter->GetOutput() );
  gradientMagnitude->Update();
  
  /*
  InternalITKImageType::Pointer gradientMagImage = gradientMagnitude->GetOutput();
  printf("gradientMagImage size %d %d\n", gradientMagImage->GetLargestPossibleRegion().GetSize()[0], gradientMagImage->GetLargestPossibleRegion().GetSize()[1] );
  printf("gradientMagImage spacing %f %f\n", gradientMagImage->GetSpacing()[0], gradientMagImage->GetSpacing()[1] );
  printf("gradientMagImage origin %f %f\n", gradientMagImage->GetOrigin()[0], gradientMagImage->GetOrigin()[1] );
  
  // Show image
  sitk::Show( sitk::Image( gradientMagnitude->GetOutput() ) );
  */
  
  //
  // StatisticsImageFilter
  //
  StatisticsFilterType::Pointer statisticsFilter = StatisticsFilterType::New();
  statisticsFilter->SetInput( gradientMagnitude->GetOutput() );
  statisticsFilter->Update();
  
  printf("Max: %f\nMin: %f\nMean: %f\nSD: %f\n\n", 
			statisticsFilter->GetMaximum(), statisticsFilter->GetMinimum(), 
			statisticsFilter->GetMean(), statisticsFilter->GetSigma());
  
  //
  // SigmoidImageFilter
  //
  SigmoidFilterType::Pointer sigmoid = SigmoidFilterType::New();
  
  sigmoid->SetOutputMinimum( 0.0 );
  sigmoid->SetOutputMaximum( 1.0 );
  sigmoid->SetBeta( atof( argv[8] ) );
  sigmoid->SetAlpha( - atof( argv[9] ) );
  
  sigmoid->SetInput( gradientMagnitude->GetOutput() );
  sigmoid->Update();
  
  
  InternalITKImageType::Pointer sigmoidImage = sigmoid->GetOutput();
  printf("sigmoidImage size %d %d\n", sigmoidImage->GetLargestPossibleRegion().GetSize()[0], sigmoidImage->GetLargestPossibleRegion().GetSize()[1] );
  printf("sigmoidImage spacing %f %f\n", sigmoidImage->GetSpacing()[0], sigmoidImage->GetSpacing()[1] );
  printf("sigmoidImage origin %f %f\n", sigmoidImage->GetOrigin()[0], sigmoidImage->GetOrigin()[1] );
  
  // Show image
  sitk::Show( sitk::Image( sigmoid->GetOutput() ) );
   

  // 
  // FastMarchingImageFilter
  //
  
  FastMarchingFilterType::Pointer fastMarching = FastMarchingFilterType::New();
  NodeContainer::Pointer seeds = NodeContainer::New();
  
  NodeType node;
  node.SetValue( 0.0 );
  InternalITKImageType::IndexType seedPosition = {350, 400};
  node.SetIndex( seedPosition );
  
  NodeType node2;
  node2.SetValue( 0.0 );
  InternalITKImageType::IndexType seedPosition2 = {350, 600};
  node2.SetIndex( seedPosition2 );

  NodeType node3;
  node3.SetValue( 0.0 );
  InternalITKImageType::IndexType seedPosition3 = {300, 500};
  node3.SetIndex( seedPosition3 );

  seeds->Initialize();
  seeds->InsertElement( 0, node );
  seeds->InsertElement( 1, node2 );
  seeds->InsertElement( 2, node3 );
  
  fastMarching->SetTrialPoints( seeds );
  fastMarching->SetNormalizationFactor( atof( argv[10] ) );
  fastMarching->SetStoppingValue( 500000 );
  fastMarching->SetOutputSize( cartesianSize );
  
  fastMarching->SetInput( sigmoid->GetOutput() );
  fastMarching->Update();
  
  /*
  InternalITKImageType::Pointer fastMarchingImage = fastMarching->GetOutput();
  printf("fastMarchingImage size %d %d\n", fastMarchingImage->GetLargestPossibleRegion().GetSize()[0], fastMarchingImage->GetLargestPossibleRegion().GetSize()[1] );
  printf("fastMarchingImage spacing %f %f\n", fastMarchingImage->GetSpacing()[0], fastMarchingImage->GetSpacing()[1] );
  printf("fastMarchingImage origin %f %f\n", fastMarchingImage->GetOrigin()[0], fastMarchingImage->GetOrigin()[1] );
  
  // Show image
  sitk::Show( sitk::Image( fastMarching->GetOutput() ) );
  */

  //
  // IntensityWindowImageFilter
  //
  IntensityFilterType::Pointer intensityWindowing = IntensityFilterType::New();

  intensityWindowing->SetWindowMinimum( 0 );
  intensityWindowing->SetWindowMaximum( 1000 );

  intensityWindowing->SetOutputMinimum(   0.0 );
  intensityWindowing->SetOutputMaximum( 1000 ); // floats but in the range of chars.

  intensityWindowing->SetInput( fastMarching->GetOutput() );

  
  InternalITKImageType::Pointer intensityWindowingImage = intensityWindowing->GetOutput();
  printf("intensityWindowingImage size %d %d\n", intensityWindowingImage->GetLargestPossibleRegion().GetSize()[0], intensityWindowingImage->GetLargestPossibleRegion().GetSize()[1] );
  printf("intensityWindowingImage spacing %f %f\n", intensityWindowingImage->GetSpacing()[0], intensityWindowingImage->GetSpacing()[1] );
  printf("intensityWindowingImage origin %f %f\n", intensityWindowingImage->GetOrigin()[0], intensityWindowingImage->GetOrigin()[1] );
  
  // Show image
  sitk::Show( sitk::Image( intensityWindowing->GetOutput() ) );
  

  //
  // StatisticsImageFilter
  //
  StatisticsFilterType::Pointer statisticsFilter2 = StatisticsFilterType::New();
  statisticsFilter2->SetInput( intensityWindowing->GetOutput() );
  statisticsFilter2->Update();
  
  printf("Max: %f\nMin: %f\nMean: %f\nSD: %f\n\n", 
			statisticsFilter2->GetMaximum(), statisticsFilter2->GetMinimum(), 
			statisticsFilter2->GetMean(), statisticsFilter2->GetSigma());
    
  //
  // BinaryThresholdImageFilter
  //
  ThresholdingFilterType::Pointer thresholder = ThresholdingFilterType::New();
  
  thresholder->SetLowerThreshold( 0.0 );
  thresholder->SetUpperThreshold( atof( argv[11] ) );
  thresholder->SetOutsideValue( 0 );
  thresholder->SetInsideValue( 255 );
  
  thresholder->SetInput( fastMarching->GetOutput() );
  thresholder->Update();

  
  InternalITKImageType::Pointer thresholderImage = thresholder->GetOutput();
  printf("thresholderImage size %d %d\n", thresholderImage->GetLargestPossibleRegion().GetSize()[0], thresholderImage->GetLargestPossibleRegion().GetSize()[1] );
  printf("thresholderImage spacing %f %f\n", thresholderImage->GetSpacing()[0], thresholderImage->GetSpacing()[1] );
  printf("thresholderImage origin %f %f\n", thresholderImage->GetOrigin()[0], thresholderImage->GetOrigin()[1] );
  
  // Show image
  sitk::Show( sitk::Image( thresholder->GetOutput() ) );
  

  /*
  sitk::ImageFileWriter writer;
  writer.SetFileName( "C:\\Users\\yingying\\Desktop\\test181.png" );
  writer.Execute( thresholder->GetOutput() );
  */
  /*
  itkImage = dynamic_cast <InternalITKImageType*> ( thresholder->GetOutput().GetITKBase() );
  outputOrigin = itkImage->GetOrigin();
  printf("origin %f %f", outputOrigin[0], outputOrigin[1] );
  */

  return EXIT_SUCCESS;
}






