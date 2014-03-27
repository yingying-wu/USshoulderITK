/*
parameters
radius of MedianImageFilter = 6
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "SimpleITK.h"

#include "itkVector.h"

#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkImageKernelOperator.h"

#include "itkRescaleIntensityImageFilter.h"
#include "itkMedianImageFilter.h"

#include "itkCartesianToPolarTransform.h"
#include "itkPolarToCartesianTransform.h"
#include "vnl/vnl_math.h"

#include "itkChangeInformationImageFilter.h"
#include "itkResampleImageFilter.h"
#include "itkNearestNeighborInterpolateImageFunction.h"
#include "itkCropImageFilter.h"

#include "itkRecursiveGaussianImageFilter.h"
#include "itkDerivativeImageFilter.h"
#include "itkAbsImageFilter.h"

#include "itkStatisticsImageFilter.h"
#include "itkBinaryThresholdImageFilter.h"
#include "itkIntensityWindowingImageFilter.h"
#include "itkSigmoidImageFilter.h"

#include "itkNeighborhoodOperator.h"
#include "itkBinaryClosingByReconstructionImageFilter.h"
#include "itkBinaryOpeningByReconstructionImageFilter.h"
#include "itkBinaryBallStructuringElement.h"
#include "itkAddImageFilter.h"

#include "itkFastMarchingImageFilter.h"

const double PI = 3.14159265359;
const double thetaResolution = PI/360;
namespace sitk = itk::simple;
typedef itk::Image< float, 2 > InternalITKImageType;
typedef itk::ImageKernelOperator< float, 2, itk::NeighborhoodAllocator< float > > ITKKernelType;


int main( int argc, char * argv[] )
{
  if( argc < 10 )
    {
    std::cerr << "Usage: " << std::endl;
	//**UPDATE
    std::cerr << argv[0] << " prostateTRUS.exe \
							filepath \
							radius(MedianImageFilter) \
							sigma(RecursiveGaussianImageFilter) \
							numberOfSD(BinaryThresholdFilter) \
							sigma(RecursiveGaussianImageFilter) \
							sigma(RecursiveGaussianImageFilter) \
							sigma(RecursiveGaussianImageFilter) \
							numberOfSD(BinaryThresholdFilter) \
							radius(BinaryClosingByReconstructionImageFilter  and BinaryOpeningByReconstructionImageFilter)"<< std::endl;
    return EXIT_FAILURE;
    }

  // declare typedefs here
  typedef itk::MedianImageFilter< InternalITKImageType, InternalITKImageType >  MedianFilterType;
  typedef itk::ChangeInformationImageFilter<InternalITKImageType> ChangeOriginFilterType;
  typedef itk::ResampleImageFilter<InternalITKImageType,InternalITKImageType> ResampleFilterType;
  typedef itk::PolarToCartesianTransform< double, 2 > PolarTransformType;
  typedef itk::NearestNeighborInterpolateImageFunction< InternalITKImageType, double > InterpolatorType;
  typedef itk::CropImageFilter< InternalITKImageType, InternalITKImageType > CropImageFilterType;
  typedef itk::RecursiveGaussianImageFilter< InternalITKImageType, InternalITKImageType > GaussianFilterType;
  typedef itk::DerivativeImageFilter< InternalITKImageType, InternalITKImageType > FilterType;
  typedef itk::CartesianToPolarTransform< double, 2 > CartesianTransformType;
  typedef itk::AbsImageFilter< InternalITKImageType, InternalITKImageType > AbsFilterType;
  typedef itk::StatisticsImageFilter< InternalITKImageType > StatisticsImageFilterType;
  typedef itk::BinaryThresholdImageFilter< InternalITKImageType, InternalITKImageType > BinaryThresholdImageFilterType;
  typedef itk::AddImageFilter< InternalITKImageType, InternalITKImageType, InternalITKImageType > AddFilterType;
  
  typedef itk::BinaryBallStructuringElement<InternalITKImageType::PixelType, InternalITKImageType::ImageDimension> StructuringElementType;
  typedef itk::BinaryClosingByReconstructionImageFilter<InternalITKImageType, StructuringElementType> ClosingFilterType;
  typedef itk::BinaryOpeningByReconstructionImageFilter<InternalITKImageType, StructuringElementType> OpeningFilterType;
  //typedef itk::IntensityWindowingImageFilter< InternalITKImageType,  InternalITKImageType > IntensityWindowingImageFilterType;
  //typedef itk::SigmoidImageFilter< InternalITKImageType,  InternalITKImageType > SigmoidImageFilterType;
  //typedef itk::FastMarchingImageFilter< InternalITKImageType,  InternalITKImageType > FastMarchingImageFilterType;
  

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
  // sitk to itk image
  //
  sitk::CastImageFilter caster;
  caster.SetOutputPixelType( sitk::sitkFloat32 );
  sitkImageIn = caster.Execute( sitkImageIn );
  
  InternalITKImageType::Pointer itkImage;
  itkImage = dynamic_cast <InternalITKImageType*> ( sitkImageIn.GetITKBase() );
  InternalITKImageType::SpacingValueType inputImageSpacing = itkImage->GetSpacing()[0];
  InternalITKImageType::SizeValueType inputImageSizeX = itkImage->GetLargestPossibleRegion().GetSize()[0];
  InternalITKImageType::SizeValueType inputImageSizeY = itkImage->GetLargestPossibleRegion().GetSize()[1];
  
  // print original image
  /*
  printf("itkImage size %d %d\n", itkImage->GetLargestPossibleRegion().GetSize()[0], itkImage->GetLargestPossibleRegion().GetSize()[1] );
  printf("itkImage spacing %f %f\n", itkImage->GetSpacing()[0], itkImage->GetSpacing()[1] );
  printf("itkImage origin %f %f\n", itkImage->GetOrigin()[0], itkImage->GetOrigin()[1] );
    
  // Show image
  sitk::Show( sitk::Image( itkImage ) );
  */
  
  //
  // MedianImageFilter 
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

  //--------------------------
  // Get long edge of prostate
  //--------------------------
  
  //
  // ChangeInformationImageFilter - shift to right-most column, middle row i.e. center of probe
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
  // everything in spatial coordinates i.e. mm
  //
  ResampleFilterType::Pointer polarResampleFilter = ResampleFilterType::New();

  PolarTransformType::Pointer polarTransform = PolarTransformType::New();
  polarResampleFilter->SetTransform( polarTransform );
  
  InterpolatorType::Pointer polarInterpolator = InterpolatorType::New();
  polarResampleFilter->SetInterpolator( polarInterpolator );

  polarResampleFilter->SetDefaultPixelValue( 0 );
  
  // make theta resolution half a degree
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
  // RecursiveGaussianImageFilter - Smooth THETA
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
  // DerivativeImageFilter - find gradient along RADIUS
  //
  FilterType::Pointer derivativeFilter = FilterType::New();

  derivativeFilter->SetOrder( 1 );
  derivativeFilter->SetDirection( 1 ); // (theta,r)

  derivativeFilter->SetInput( gaussianFilter->GetOutput() );  
  derivativeFilter->Update();
  
  /*
  InternalITKImageType::Pointer derivativeImage = derivativeFilter->GetOutput();
  printf("derivativeImage size %d %d\n", derivativeImage->GetLargestPossibleRegion().GetSize()[0], derivativeImage->GetLargestPossibleRegion().GetSize()[1] );
  printf("derivativeImage spacing %f %f\n", derivativeImage->GetSpacing()[0], derivativeImage->GetSpacing()[1] );
  printf("derivativeImage origin %f %f\n", derivativeImage->GetOrigin()[0], derivativeImage->GetOrigin()[1] );
  
  // Show image
  sitk::Show( sitk::Image( derivativeFilter->GetOutput() ) );
  */

  //
  // ResampleImageFilter - CartesianToPolarTransform 
  //
  ResampleFilterType::Pointer cartesianResampleFilter = ResampleFilterType::New();
 
  CartesianTransformType::Pointer cartesianTransform = CartesianTransformType::New();
  cartesianResampleFilter->SetTransform( cartesianTransform );
  
  InterpolatorType::Pointer cartesianInterpolator = InterpolatorType::New();
  cartesianResampleFilter->SetInterpolator( cartesianInterpolator );
  
  cartesianResampleFilter->SetDefaultPixelValue( 0 );
  
  // make it back to same spacing as input image
  double cartesianSpacing[ 2 ] = { inputImageSpacing, inputImageSpacing }; // (x, y) pixel spacing in millimeters
  cartesianResampleFilter->SetOutputSpacing( cartesianSpacing );
  
  double cartesianOrigin[ 2 ]; // (x, y) space coordinate of origin 
  cartesianOrigin[0] = probeOrigin[0]; 
  cartesianOrigin[1] = probeOrigin[1]; 
  cartesianResampleFilter->SetOutputOrigin( cartesianOrigin );
  
  InternalITKImageType::SizeType cartesianSize = { inputImageSizeX, inputImageSizeY }; // (x, y) number of pixels 
  cartesianResampleFilter->SetSize( cartesianSize );
  
  cartesianResampleFilter->SetInput( derivativeFilter->GetOutput() );
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
  // ResampleImageFilter - CartesianToPolarTransform after blur along theta
  //
  ResampleFilterType::Pointer cartesianResampleFilterBlurred = ResampleFilterType::New();
 
  CartesianTransformType::Pointer cartesianTransformBlurred = CartesianTransformType::New();
  cartesianResampleFilterBlurred->SetTransform( cartesianTransformBlurred );
  
  InterpolatorType::Pointer cartesianInterpolatorBlurred = InterpolatorType::New();
  cartesianResampleFilterBlurred->SetInterpolator( cartesianInterpolatorBlurred );
  
  cartesianResampleFilterBlurred->SetDefaultPixelValue( 0 );
  
  // use same setting as the transformation after derivative
  cartesianResampleFilterBlurred->SetOutputSpacing( cartesianSpacing );  
  cartesianResampleFilterBlurred->SetOutputOrigin( cartesianOrigin );
  cartesianResampleFilterBlurred->SetSize( cartesianSize );
  
  cartesianResampleFilterBlurred->SetInput( gaussianFilter->GetOutput() );
  cartesianResampleFilterBlurred->Update();
  
  
  InternalITKImageType::Pointer cartesianImageBlurred = cartesianResampleFilterBlurred->GetOutput();
  printf("cartesianImageBlurred size %d %d\n", cartesianImageBlurred->GetLargestPossibleRegion().GetSize()[0], cartesianImageBlurred->GetLargestPossibleRegion().GetSize()[1] );
  printf("cartesianImageBlurred spacing %f %f\n", cartesianImageBlurred->GetSpacing()[0], cartesianImageBlurred->GetSpacing()[1] );
  printf("cartesianImageBlurred origin %f %f\n", cartesianImageBlurred->GetOrigin()[0], cartesianImageBlurred->GetOrigin()[1] );
  
  // Show image
  sitk::Show( sitk::Image( cartesianResampleFilterBlurred->GetOutput() ) );
  

  //
  // ChangeInformationImageFilter - shift to lower left corner (0,0)
  //
  ChangeOriginFilterType::Pointer cartesianOriginFilter = ChangeOriginFilterType::New();
  
  // value where the output image origin should start i.e. start value of output image sampling
  InternalITKImageType::PointType finalImageOrigin;
  finalImageOrigin[0] = 0;
  finalImageOrigin[1] = 0;
  cartesianOriginFilter->ChangeOriginOn();
  cartesianOriginFilter->SetOutputOrigin( finalImageOrigin );
  printf("origin %f %f\n", finalImageOrigin[0], finalImageOrigin[1] );

  cartesianOriginFilter->SetInput( cartesianResampleFilter->GetOutput() );
  cartesianOriginFilter->Update();
  
  /*
  InternalITKImageType::Pointer cartesianOriginImage = cartesianOriginFilter->GetOutput();
  printf("cartesianOriginImage size %d %d\n", cartesianOriginImage->GetLargestPossibleRegion().GetSize()[0], cartesianOriginImage->GetLargestPossibleRegion().GetSize()[1] );
  printf("cartesianOriginImage spacing %f %f\n", cartesianOriginImage->GetSpacing()[0], cartesianOriginImage->GetSpacing()[1] );
  printf("cartesianOriginImage origin %f %f\n", cartesianOriginImage->GetOrigin()[0], cartesianOriginImage->GetOrigin()[1] );
  
  // Show image
  sitk::Show( sitk::Image( cartesianOriginFilter->GetOutput() ) );
  */

  //
  // AbsImageFilter
  //
  AbsFilterType::Pointer absFilter = AbsFilterType::New();
  absFilter->SetInput( cartesianOriginFilter->GetOutput() );
  absFilter->Update();
  
  /*
  InternalITKImageType::Pointer absImage = absFilter->GetOutput();
  printf("absImage size %d %d\n", absImage->GetLargestPossibleRegion().GetSize()[0], absImage->GetLargestPossibleRegion().GetSize()[1] );
  printf("absImage spacing %f %f\n", absImage->GetSpacing()[0], absImage->GetSpacing()[1] );
  printf("absImage origin %f %f\n", absImage->GetOrigin()[0], absImage->GetOrigin()[1] );
  
  // Show image
  sitk::Show( sitk::Image( absFilter->GetOutput() ) );
  */

  //
  // StatisticsImageFilter
  //
  StatisticsImageFilterType::Pointer statisticsFilter = StatisticsImageFilterType::New();
  statisticsFilter->SetInput( absFilter->GetOutput() );
  statisticsFilter->Update();
  
  printf("Max: %f\nMin: %f\nMean: %f\nSD: %f\n\n", 
			statisticsFilter->GetMaximum(), statisticsFilter->GetMinimum(), 
			statisticsFilter->GetMean(), statisticsFilter->GetSigma());
  
  // 
  // BinaryThresholdImageFilter
  //
  BinaryThresholdImageFilterType::Pointer binaryFilter = BinaryThresholdImageFilterType::New();

  binaryFilter->SetInsideValue( 255 );
  binaryFilter->SetOutsideValue( 0 );
  binaryFilter->SetLowerThreshold( statisticsFilter->GetMean() + atof(argv[4]) * statisticsFilter->GetSigma() );
  binaryFilter->SetUpperThreshold( statisticsFilter->GetMaximum() );

  binaryFilter->SetInput( absFilter->GetOutput() );
  binaryFilter->Update();
  
  
  InternalITKImageType::Pointer binaryImage = binaryFilter->GetOutput();
  printf("binaryImage size %d %d\n", binaryImage->GetLargestPossibleRegion().GetSize()[0], binaryImage->GetLargestPossibleRegion().GetSize()[1] );
  printf("binaryImage spacing %f %f\n", binaryImage->GetSpacing()[0], binaryImage->GetSpacing()[1] );
  printf("binaryImage origin %f %f\n", binaryImage->GetOrigin()[0], binaryImage->GetOrigin()[1] );
  
  // Show image
  sitk::Show( sitk::Image( binaryFilter->GetOutput() ) );
  

  //---------------------------
  // Get short edge of prostate after long edge smoothing!!
  //---------------------------

  //
  // RecursiveGaussianImageFilter - Smooth THETA
  //
  
  GaussianFilterType::Pointer gaussianFilter11 = GaussianFilterType::New();

  gaussianFilter11->SetSigma( atof(argv[5]) ); 
  gaussianFilter11->SetDirection( 0 ); // (theta, r)
  gaussianFilter11->SetNormalizeAcrossScale( true );

  gaussianFilter11->SetInput( polarResampleFilter->GetOutput() );
  gaussianFilter11->Update();
  

  /*
  InternalITKImageType::Pointer gaussianImage11 = gaussianFilter11->GetOutput();
  printf("gaussianImage11 size %d %d\n", gaussianImage11->GetLargestPossibleRegion().GetSize()[0], gaussianImage11->GetLargestPossibleRegion().GetSize()[1] );
  printf("gaussianImage11 spacing %f %f\n", gaussianImage11->GetSpacing()[0], gaussianImage11->GetSpacing()[1] );
  printf("gaussianImage11 origin %f %f\n", gaussianImage11->GetOrigin()[0], gaussianImage11->GetOrigin()[1] );
  
  // Show image
  sitk::Show( sitk::Image( gaussianFilter11->GetOutput() ) );
  */

  //
  // ResampleImageFilter - CartesianToPolarTransform after blur along theta
  //
  
  ResampleFilterType::Pointer cartesianResampleFilterBlurred11 = ResampleFilterType::New();
 
  CartesianTransformType::Pointer cartesianTransformBlurred11 = CartesianTransformType::New();
  cartesianResampleFilterBlurred11->SetTransform( cartesianTransformBlurred11 );
  
  InterpolatorType::Pointer cartesianInterpolatorBlurred11 = InterpolatorType::New();
  cartesianResampleFilterBlurred11->SetInterpolator( cartesianInterpolatorBlurred11 );
  
  cartesianResampleFilterBlurred11->SetDefaultPixelValue( 0 );
  
  // use same setting as the transformation after derivative
  cartesianResampleFilterBlurred11->SetOutputSpacing( cartesianSpacing );  
  cartesianResampleFilterBlurred11->SetOutputOrigin( cartesianOrigin );
  cartesianResampleFilterBlurred11->SetSize( cartesianSize );
  
  cartesianResampleFilterBlurred11->SetInput( gaussianFilter11->GetOutput() );
  cartesianResampleFilterBlurred11->Update();
  

  /*
  InternalITKImageType::Pointer cartesianImageBlurred11 = cartesianResampleFilterBlurred11->GetOutput();
  printf("cartesianImageBlurred11 size %d %d\n", cartesianImageBlurred11->GetLargestPossibleRegion().GetSize()[0], cartesianImageBlurred11->GetLargestPossibleRegion().GetSize()[1] );
  printf("cartesianImageBlurred11 spacing %f %f\n", cartesianImageBlurred11->GetSpacing()[0], cartesianImageBlurred11->GetSpacing()[1] );
  printf("cartesianImageBlurred11 origin %f %f\n", cartesianImageBlurred11->GetOrigin()[0], cartesianImageBlurred11->GetOrigin()[1] );
  
  // Show image
  sitk::Show( sitk::Image( cartesianResampleFilterBlurred11->GetOutput() ) );
  */

  //
  // ChangeInformationImageFilter - shift to middle column, middle row i.e. center of prostate
  //
  
  ChangeOriginFilterType::Pointer prostateOriginFilter11 = ChangeOriginFilterType::New();
  
  // value where the output image origin should start i.e. start value of output image sampling
  InternalITKImageType::PointType prostateOrigin11;
  prostateOrigin11[0] = -( floor( inputImageSizeX / 2 ) - inputImageSpacing ); //-( floor( 3 * inputImageSizeX / 5 ) - inputImageSpacing );
  prostateOrigin11[1] = -( floor( inputImageSizeY / 2 ) - inputImageSpacing );
  prostateOriginFilter11->ChangeOriginOn();
  prostateOriginFilter11->SetOutputOrigin( prostateOrigin11 );
  printf("origin %f %f\n", prostateOrigin11[0], prostateOrigin11[1] );

  prostateOriginFilter11->SetInput( cartesianResampleFilterBlurred11->GetOutput() );
  prostateOriginFilter11->Update();
  

  /*
  InternalITKImageType::Pointer prostateOriginImage11 = prostateOriginFilter11->GetOutput();
  printf("prostateOriginImage11 size %d %d\n", prostateOriginImage11->GetLargestPossibleRegion().GetSize()[0], prostateOriginImage11->GetLargestPossibleRegion().GetSize()[1] );
  printf("prostateOriginImage11 spacing %f %f\n", prostateOriginImage11->GetSpacing()[0], prostateOriginImage11->GetSpacing()[1] );
  printf("prostateOriginImage11 origin %f %f\n", prostateOriginImage11->GetOrigin()[0], prostateOriginImage11->GetOrigin()[1] );
  
  // Show image
  sitk::Show( sitk::Image( prostateOriginFilter11->GetOutput() ) );
  */
  
  //
  // ResampleImageFilter - PolarToCartesianTransform2
  // everything in spatial coordinates i.e. mm
  //

  ResampleFilterType::Pointer polarResampleFilter12 = ResampleFilterType::New();
  
  PolarTransformType::Pointer polarTransform12 = PolarTransformType::New();
  polarResampleFilter12->SetTransform( polarTransform12 );
  
  InterpolatorType::Pointer polarInterpolator12 = InterpolatorType::New();
  polarResampleFilter12->SetInterpolator( polarInterpolator12 );
    
  polarResampleFilter12->SetDefaultPixelValue( 0 );

  // make theta resolution half a degree
  // make r the same spatial resolution as input image
  double polarSpacing12[ 2 ] = { thetaResolution, inputImageSpacing }; // (theta, r) pixel spacing in millimeters
  polarResampleFilter12->SetOutputSpacing( polarSpacing12 );
  
  // value where the output image origin should start i.e. start value of output image sampling
  // theta range from 0 to 2*pi since origin is in the center
  double polarOrigin12[ 2 ] = { 0.0, 0.0 }; // (theta, r) space coordinate of origin
  polarResampleFilter12->SetOutputOrigin( polarOrigin12 );

  InternalITKImageType::SizeType polarSize12 = { 720, 512 }; // (theta, r) number of pixels
  polarResampleFilter12->SetSize( polarSize12 );
  
  polarResampleFilter12->SetInput( prostateOriginFilter11->GetOutput() );
  polarResampleFilter12->Update();
  

  /*
  InternalITKImageType::Pointer polarImage12 = polarResampleFilter12->GetOutput();
  printf("polarImage12 size %d %d\n", polarImage12->GetLargestPossibleRegion().GetSize()[0], polarImage12->GetLargestPossibleRegion().GetSize()[1] );
  printf("polarImage12 spacing %f %f\n", polarImage12->GetSpacing()[0], polarImage12->GetSpacing()[1] );
  printf("polarImage12 origin %f %f\n", polarImage12->GetOrigin()[0], polarImage12->GetOrigin()[1] );
  
  // Show image
  sitk::Show( sitk::Image( polarResampleFilter12->GetOutput() ) );
  */

  //
  // RecursiveGaussianImageFilter - Smooth THETA
  //
  
  GaussianFilterType::Pointer gaussianFilter12 = GaussianFilterType::New();

  gaussianFilter12->SetSigma( atof(argv[6]) );
  gaussianFilter12->SetDirection( 0 ); //(theta, r)
  gaussianFilter12->SetNormalizeAcrossScale( true );

  gaussianFilter12->SetInput( polarResampleFilter12->GetOutput() );
  gaussianFilter12->Update();
  

  /*
  InternalITKImageType::Pointer gaussianImage12 = gaussianFilter12->GetOutput();
  printf("gaussianImage12 size %d %d\n", gaussianImage12->GetLargestPossibleRegion().GetSize()[0], gaussianImage12->GetLargestPossibleRegion().GetSize()[1] );
  printf("gaussianImage12 spacing %f %f\n", gaussianImage12->GetSpacing()[0], gaussianImage12->GetSpacing()[1] );
  printf("gaussianImage12 origin %f %f\n", gaussianImage12->GetOrigin()[0], gaussianImage12->GetOrigin()[1] );
  
  // Show image
  sitk::Show( sitk::Image( gaussianFilter12->GetOutput() ) );
  */
  
  // 
  // DerivativeImageFilter - find gradient along RADIUS
  //
  
  FilterType::Pointer derivativeFilter12 = FilterType::New();

  derivativeFilter12->SetOrder( 1 );
  derivativeFilter12->SetDirection( 1 ); //(theta,r)

  derivativeFilter12->SetInput( gaussianFilter12->GetOutput() );  
  derivativeFilter12->Update();
  

  /*
  InternalITKImageType::Pointer derivativeImage12 = derivativeFilter12->GetOutput();
  printf("derivativeImage12 size %d %d\n", derivativeImage12->GetLargestPossibleRegion().GetSize()[0], derivativeImage12->GetLargestPossibleRegion().GetSize()[1] );
  printf("derivativeImage12 spacing %f %f\n", derivativeImage12->GetSpacing()[0], derivativeImage12->GetSpacing()[1] );
  printf("derivativeImage12 origin %f %f\n", derivativeImage12->GetOrigin()[0], derivativeImage12->GetOrigin()[1] );
  
  // Show image
  sitk::Show( sitk::Image( derivativeFilter12->GetOutput() ) );
  */
  
  //
  // ResampleImageFilter - CartesianToPolarTransform 
  //
  
  ResampleFilterType::Pointer cartesianResampleFilter12 = ResampleFilterType::New();
 
  CartesianTransformType::Pointer cartesianTransform12 = CartesianTransformType::New();
  cartesianResampleFilter12->SetTransform( cartesianTransform12 );
  
  InterpolatorType::Pointer cartesianInterpolator12 = InterpolatorType::New();
  cartesianResampleFilter12->SetInterpolator( cartesianInterpolator12 );
  
  cartesianResampleFilter12->SetDefaultPixelValue( 0 );
  
  // make it back to same spacing as input image
  double cartesianSpacing12[ 2 ] = { inputImageSpacing, inputImageSpacing }; // (x, y) pixel spacing in millimeters
  cartesianResampleFilter12->SetOutputSpacing( cartesianSpacing12 );
  
  double cartesianOrigin12[ 2 ]; // (x, y) space coordinate of origin 
  cartesianOrigin12[0] = prostateOrigin11[0]; 
  cartesianOrigin12[1] = prostateOrigin11[1]; 
  cartesianResampleFilter12->SetOutputOrigin( cartesianOrigin12 );

  InternalITKImageType::SizeType cartesianSize12 = { inputImageSizeX, inputImageSizeY }; // (x, y) number of pixels 
  cartesianResampleFilter12->SetSize( cartesianSize12 );
  
  cartesianResampleFilter12->SetInput( derivativeFilter12->GetOutput() );
  cartesianResampleFilter12->Update();
  

  /*
  InternalITKImageType::Pointer cartesianImage12 = cartesianResampleFilter12->GetOutput();
  printf("cartesianImage12 size %d %d\n", cartesianImage12->GetLargestPossibleRegion().GetSize()[0], cartesianImage12->GetLargestPossibleRegion().GetSize()[1] );
  printf("cartesianImage12 spacing %f %f\n", cartesianImage12->GetSpacing()[0], cartesianImage12->GetSpacing()[1] );
  printf("cartesianImage12 origin %f %f\n", cartesianImage12->GetOrigin()[0], cartesianImage12->GetOrigin()[1] );
  
  // Show image
  sitk::Show( sitk::Image( cartesianResampleFilter12->GetOutput() ) );
  */
  
  
  //
  // ResampleImageFilter - CartesianToPolarTransform after blur along theta
  //
  
  ResampleFilterType::Pointer cartesianResampleFilterBlurred12 = ResampleFilterType::New();
 
  CartesianTransformType::Pointer cartesianTransformBlurred12 = CartesianTransformType::New();
  cartesianResampleFilterBlurred12->SetTransform( cartesianTransformBlurred12 );
  
  InterpolatorType::Pointer cartesianInterpolatorBlurred12 = InterpolatorType::New();
  cartesianResampleFilterBlurred12->SetInterpolator( cartesianInterpolatorBlurred12 );
  
  cartesianResampleFilterBlurred12->SetDefaultPixelValue( 0 );
  
  // use same setting as the transformation after derivative
  cartesianResampleFilterBlurred12->SetOutputSpacing( cartesianSpacing12 );
  cartesianResampleFilterBlurred12->SetOutputOrigin( cartesianOrigin12 );
  cartesianResampleFilterBlurred12->SetSize( cartesianSize12 );
  
  cartesianResampleFilterBlurred12->SetInput( gaussianFilter12->GetOutput() );
  cartesianResampleFilterBlurred12->Update();
  

  /*
  InternalITKImageType::Pointer cartesianImageBlurred12 = cartesianResampleFilterBlurred12->GetOutput();
  printf("cartesianImageBlurred12 size %d %d\n", cartesianImageBlurred12->GetLargestPossibleRegion().GetSize()[0], cartesianImageBlurred12->GetLargestPossibleRegion().GetSize()[1] );
  printf("cartesianImageBlurred12 spacing %f %f\n", cartesianImageBlurred12->GetSpacing()[0], cartesianImageBlurred12->GetSpacing()[1] );
  printf("cartesianImageBlurred12 origin %f %f\n", cartesianImageBlurred12->GetOrigin()[0], cartesianImageBlurred12->GetOrigin()[1] );
  
  // Show image
  sitk::Show( sitk::Image( cartesianResampleFilterBlurred12->GetOutput() ) );
  */

  //---------------------------
  // Get short edge of prostate - centered at probe
  //---------------------------

  //
  // ChangeInformationImageFilter - shift to right-most column, middle row i.e. center of probe
  //
  ChangeOriginFilterType::Pointer probeOriginFilter41 = ChangeOriginFilterType::New();
  
  // value where the output image origin should start i.e. start value of output image sampling
  InternalITKImageType::PointType probeOrigin41;
  probeOrigin41[0] = -( inputImageSizeX - inputImageSpacing );
  probeOrigin41[1] = -( floor( inputImageSizeY / 2 ) - inputImageSpacing );
  probeOriginFilter41->ChangeOriginOn();
  probeOriginFilter41->SetOutputOrigin( probeOrigin41 );

  probeOriginFilter41->SetInput( cartesianResampleFilterBlurred12->GetOutput() );
  probeOriginFilter41->Update();
  
  /*
  InternalITKImageType::Pointer probeOriginImage41 = probeOriginFilter41->GetOutput();
  printf("probeOriginImage41 size %d %d\n", probeOriginImage41->GetLargestPossibleRegion().GetSize()[0], probeOriginImage41->GetLargestPossibleRegion().GetSize()[1] );
  printf("probeOriginImage41 spacing %f %f\n", probeOriginImage41->GetSpacing()[0], probeOriginImage41->GetSpacing()[1] );
  printf("probeOriginImage41 origin %f %f\n", probeOriginImage41->GetOrigin()[0], probeOriginImage41->GetOrigin()[1] );
  
  // Show image
  sitk::Show( sitk::Image( probeOriginFilter41->GetOutput() ) );
  */

  //
  // ResampleImageFilter - PolarToCartesianTransform
  // everything in spatial coordinates i.e. mm
  //
  ResampleFilterType::Pointer polarResampleFilter41 = ResampleFilterType::New();

  PolarTransformType::Pointer polarTransform41 = PolarTransformType::New();
  polarResampleFilter41->SetTransform( polarTransform41 );
  
  InterpolatorType::Pointer polarInterpolator41 = InterpolatorType::New();
  polarResampleFilter41->SetInterpolator( polarInterpolator41 );

  polarResampleFilter41->SetDefaultPixelValue( 0 );
  
  // make theta resolution half a degree
  // make r the same spatial resolution as input image
  double polarSpacing41[ 2 ] = { thetaResolution, inputImageSpacing }; // (theta, r) pixel spacing in millimeters
  polarResampleFilter41->SetOutputSpacing( polarSpacing41 );
  
  // value where the output image origin should start i.e. start value of output image sampling
  // theta range from pi/2 to 3*pi/2 since origin is on rightmost column
  double polarOrigin41[ 2 ] = { PI/2, 0.0 }; // (theta, r) space coordinate of origin
  polarResampleFilter41->SetOutputOrigin( polarOrigin41 );

  InternalITKImageType::SizeType polarSize41 = { 360, 512 }; // (theta, r) number of pixels
  polarResampleFilter41->SetSize( polarSize41 );
  
  polarResampleFilter41->SetInput( probeOriginFilter41->GetOutput() );
  polarResampleFilter41->Update();
  
  /*
  InternalITKImageType::Pointer polarImage41 = polarResampleFilter41->GetOutput();
  printf("polarImage41 size %d %d\n", polarImage41->GetLargestPossibleRegion().GetSize()[0], polarImage41->GetLargestPossibleRegion().GetSize()[1] );
  printf("polarImage41 spacing %f %f\n", polarImage41->GetSpacing()[0], polarImage41->GetSpacing()[1] );
  printf("polarImage41 origin %f %f\n", polarImage41->GetOrigin()[0], polarImage41->GetOrigin()[1] );
  
  // Show image
  sitk::Show( sitk::Image( polarResampleFilter41->GetOutput() ) );
  */

  //
  // RecursiveGaussianImageFilter - Smooth RADIUS
  //
  GaussianFilterType::Pointer gaussianFilter31 = GaussianFilterType::New();

  gaussianFilter31->SetSigma( atof(argv[7]) ); 
  gaussianFilter31->SetDirection( 1 ); // (theta, r)
  gaussianFilter31->SetNormalizeAcrossScale( true );

  gaussianFilter31->SetInput( polarResampleFilter41->GetOutput() );
  gaussianFilter31->Update();
  
  /*
  InternalITKImageType::Pointer gaussianImage31 = gaussianFilter31->GetOutput();
  printf("gaussianImage31 size %d %d\n", gaussianImage31->GetLargestPossibleRegion().GetSize()[0], gaussianImage31->GetLargestPossibleRegion().GetSize()[1] );
  printf("gaussianImage31 spacing %f %f\n", gaussianImage31->GetSpacing()[0], gaussianImage31->GetSpacing()[1] );
  printf("gaussianImage31 origin %f %f\n", gaussianImage31->GetOrigin()[0], gaussianImage31->GetOrigin()[1] );
  
  // Show image
  sitk::Show( sitk::Image( gaussianFilter31->GetOutput() ) );
  */

  // 
  // DerivativeImageFilter - find gradient along THETA
  //
  FilterType::Pointer derivativeFilter31 = FilterType::New();

  derivativeFilter31->SetOrder( 1 );
  derivativeFilter31->SetDirection( 0 ); // (theta,r)

  derivativeFilter31->SetInput( gaussianFilter31->GetOutput() );  
  derivativeFilter31->Update();
  
  /*
  InternalITKImageType::Pointer derivativeImage31 = derivativeFilter31->GetOutput();
  printf("derivativeImage31 size %d %d\n", derivativeImage31->GetLargestPossibleRegion().GetSize()[0], derivativeImage31->GetLargestPossibleRegion().GetSize()[1] );
  printf("derivativeImage31 spacing %f %f\n", derivativeImage31->GetSpacing()[0], derivativeImage31->GetSpacing()[1] );
  printf("derivativeImage31 origin %f %f\n", derivativeImage31->GetOrigin()[0], derivativeImage31->GetOrigin()[1] );
  
  // Show image
  sitk::Show( sitk::Image( derivativeFilter31->GetOutput() ) );
  */

  //
  // ResampleImageFilter - CartesianToPolarTransform 
  //
  ResampleFilterType::Pointer cartesianResampleFilter31 = ResampleFilterType::New();
 
  CartesianTransformType::Pointer cartesianTransform31 = CartesianTransformType::New();
  cartesianResampleFilter31->SetTransform( cartesianTransform31 );
  
  InterpolatorType::Pointer cartesianInterpolator31 = InterpolatorType::New();
  cartesianResampleFilter31->SetInterpolator( cartesianInterpolator31 );
  
  cartesianResampleFilter31->SetDefaultPixelValue( 0 );
  
  // make it back to same spacing as input image
  double cartesianSpacing31[ 2 ] = { inputImageSpacing, inputImageSpacing }; // (x, y) pixel spacing in millimeters
  cartesianResampleFilter31->SetOutputSpacing( cartesianSpacing31 );
  
  double cartesianOrigin31[ 2 ]; // (x, y) space coordinate of origin 
  cartesianOrigin31[0] = probeOrigin[0]; 
  cartesianOrigin31[1] = probeOrigin[1]; 
  cartesianResampleFilter31->SetOutputOrigin( cartesianOrigin31 );
  
  InternalITKImageType::SizeType cartesianSize31 = { inputImageSizeX, inputImageSizeY }; // (x, y) number of pixels 
  cartesianResampleFilter31->SetSize( cartesianSize31 );
  
  cartesianResampleFilter31->SetInput( derivativeFilter31->GetOutput() );
  cartesianResampleFilter31->Update();
  
  /*
  InternalITKImageType::Pointer cartesianImage31 = cartesianResampleFilter31->GetOutput();
  printf("cartesianImage31 size %d %d\n", cartesianImage31->GetLargestPossibleRegion().GetSize()[0], cartesianImage31->GetLargestPossibleRegion().GetSize()[1] );
  printf("cartesianImage31 spacing %f %f\n", cartesianImage31->GetSpacing()[0], cartesianImage31->GetSpacing()[1] );
  printf("cartesianImage31 origin %f %f\n", cartesianImage31->GetOrigin()[0], cartesianImage31->GetOrigin()[1] );
  
  // Show image
  sitk::Show( sitk::Image( cartesianResampleFilter31->GetOutput() ) );
  */

  //
  // ResampleImageFilter - CartesianToPolarTransform after blur along theta
  //
  ResampleFilterType::Pointer cartesianResampleFilterBlurred31 = ResampleFilterType::New();
 
  CartesianTransformType::Pointer cartesianTransformBlurred31 = CartesianTransformType::New();
  cartesianResampleFilterBlurred31->SetTransform( cartesianTransformBlurred31 );
  
  InterpolatorType::Pointer cartesianInterpolatorBlurred31 = InterpolatorType::New();
  cartesianResampleFilterBlurred31->SetInterpolator( cartesianInterpolatorBlurred31 );
  
  cartesianResampleFilterBlurred31->SetDefaultPixelValue( 0 );
  
  // use same setting as the transformation after derivative
  cartesianResampleFilterBlurred31->SetOutputSpacing( cartesianSpacing31 );  
  cartesianResampleFilterBlurred31->SetOutputOrigin( cartesianOrigin31 );
  cartesianResampleFilterBlurred31->SetSize( cartesianSize31 );
  
  cartesianResampleFilterBlurred31->SetInput( gaussianFilter31->GetOutput() );
  cartesianResampleFilterBlurred31->Update();
  
  
  InternalITKImageType::Pointer cartesianImageBlurred31 = cartesianResampleFilterBlurred31->GetOutput();
  printf("cartesianImageBlurred31 size %d %d\n", cartesianImageBlurred31->GetLargestPossibleRegion().GetSize()[0], cartesianImageBlurred31->GetLargestPossibleRegion().GetSize()[1] );
  printf("cartesianImageBlurred31 spacing %f %f\n", cartesianImageBlurred31->GetSpacing()[0], cartesianImageBlurred31->GetSpacing()[1] );
  printf("cartesianImageBlurred31 origin %f %f\n", cartesianImageBlurred31->GetOrigin()[0], cartesianImageBlurred31->GetOrigin()[1] );
  
  // Show image
  sitk::Show( sitk::Image( cartesianResampleFilterBlurred31->GetOutput() ) );
  

  //
  // ChangeInformationImageFilter - shift to lower left corner (0,0)
  //
  ChangeOriginFilterType::Pointer cartesianOriginFilter31 = ChangeOriginFilterType::New();
  
  // value where the output image origin should start i.e. start value of output image sampling
  InternalITKImageType::PointType finalImageOrigin31;
  finalImageOrigin31[0] = 0;
  finalImageOrigin31[1] = 0;
  cartesianOriginFilter31->ChangeOriginOn();
  cartesianOriginFilter31->SetOutputOrigin( finalImageOrigin31 );
  printf("origin %f %f\n", finalImageOrigin31[0], finalImageOrigin31[1] );

  cartesianOriginFilter31->SetInput( cartesianResampleFilter31->GetOutput() );
  cartesianOriginFilter31->Update();
  
  /*
  InternalITKImageType::Pointer cartesianOriginImage31 = cartesianOriginFilter31->GetOutput();
  printf("cartesianOriginImage31 size %d %d\n", cartesianOriginImage31->GetLargestPossibleRegion().GetSize()[0], cartesianOriginImage31->GetLargestPossibleRegion().GetSize()[1] );
  printf("cartesianOriginImage31 spacing %f %f\n", cartesianOriginImage31->GetSpacing()[0], cartesianOriginImage31->GetSpacing()[1] );
  printf("cartesianOriginImage31 origin %f %f\n", cartesianOriginImage31->GetOrigin()[0], cartesianOriginImage31->GetOrigin()[1] );
  
  // Show image
  sitk::Show( sitk::Image( cartesianOriginFilter31->GetOutput() ) );
  */

  //
  // AbsImageFilter
  //
  AbsFilterType::Pointer absFilter31 = AbsFilterType::New();
  absFilter31->SetInput( cartesianOriginFilter31->GetOutput() );
  absFilter31->Update();
  
  /*
  InternalITKImageType::Pointer absImage31 = absFilter31->GetOutput();
  printf("absImage31 size %d %d\n", absImage31->GetLargestPossibleRegion().GetSize()[0], absImage31->GetLargestPossibleRegion().GetSize()[1] );
  printf("absImage31 spacing %f %f\n", absImage31->GetSpacing()[0], absImage31->GetSpacing()[1] );
  printf("absImage31 origin %f %f\n", absImage31->GetOrigin()[0], absImage31->GetOrigin()[1] );
  
  // Show image
  sitk::Show( sitk::Image( absFilter31->GetOutput() ) );
  */

  //
  // StatisticsImageFilter
  //
  StatisticsImageFilterType::Pointer statisticsFilter31 = StatisticsImageFilterType::New();
  statisticsFilter31->SetInput( absFilter31->GetOutput() );
  statisticsFilter31->Update();
  
  printf("Max: %f\nMin: %f\nMean: %f\nSD: %f\n\n", 
			statisticsFilter31->GetMaximum(), statisticsFilter31->GetMinimum(), 
			statisticsFilter31->GetMean(), statisticsFilter31->GetSigma());
  
  // 
  // BinaryThresholdImageFilter
  //
  BinaryThresholdImageFilterType::Pointer binaryFilter31 = BinaryThresholdImageFilterType::New();

  binaryFilter31->SetInsideValue( 255 );
  binaryFilter31->SetOutsideValue( 0 );
  binaryFilter31->SetLowerThreshold( statisticsFilter31->GetMean() + atof(argv[8]) * statisticsFilter31->GetSigma() );
  binaryFilter31->SetUpperThreshold( statisticsFilter31->GetMaximum() );

  binaryFilter31->SetInput( absFilter31->GetOutput() );
  binaryFilter31->Update();

  
  InternalITKImageType::Pointer binaryImage31 = binaryFilter31->GetOutput();
  printf("binaryImage31 size %d %d\n", binaryImage31->GetLargestPossibleRegion().GetSize()[0], binaryImage31->GetLargestPossibleRegion().GetSize()[1] );
  printf("binaryImage31 spacing %f %f\n", binaryImage31->GetSpacing()[0], binaryImage31->GetSpacing()[1] );
  printf("binaryImage31 origin %f %f\n", binaryImage31->GetOrigin()[0], binaryImage31->GetOrigin()[1] );
  
  // Show image
  sitk::Show( sitk::Image( binaryFilter31->GetOutput() ) );

  
  //
  // BinaryClosingByReconstructionImageFilter
  //
  StructuringElementType structuringElement;
  structuringElement.SetRadius(atoi( argv[9] ));
  structuringElement.CreateStructuringElement();
   
  ClosingFilterType::Pointer closingFilter = ClosingFilterType::New();
  
  closingFilter->SetKernel( structuringElement );
  closingFilter->SetForegroundValue( 255 );

  closingFilter->SetInput( binaryFilter31->GetOutput() );
  closingFilter->Update();
  
  
  InternalITKImageType::Pointer closingImage = closingFilter->GetOutput();
  printf("closingImage size %d %d\n", closingImage->GetLargestPossibleRegion().GetSize()[0], closingImage->GetLargestPossibleRegion().GetSize()[1] );
  printf("closingImage spacing %f %f\n", closingImage->GetSpacing()[0], closingImage->GetSpacing()[1] );
  printf("closingImage origin %f %f\n", closingImage->GetOrigin()[0], closingImage->GetOrigin()[1] );
  
  // Show image
  sitk::Show( sitk::Image( closingFilter->GetOutput() ) );
  
  //
  // BinaryOpeningByReconstructionImageFilter
  //
  OpeningFilterType::Pointer openingFilter = OpeningFilterType::New();
  
  openingFilter->SetKernel( structuringElement );
  openingFilter->SetForegroundValue( 255 );

  openingFilter->SetInput( closingFilter->GetOutput() );
  openingFilter->Update();
  
  
  InternalITKImageType::Pointer openingImage = openingFilter->GetOutput();
  printf("openingImage size %d %d\n", openingImage->GetLargestPossibleRegion().GetSize()[0], openingImage->GetLargestPossibleRegion().GetSize()[1] );
  printf("openingImage spacing %f %f\n", openingImage->GetSpacing()[0], openingImage->GetSpacing()[1] );
  printf("openingImage origin %f %f\n", openingImage->GetOrigin()[0], openingImage->GetOrigin()[1] );
  
  // Show image
  sitk::Show( sitk::Image( openingFilter->GetOutput() ) );
  
  //
  // AddImageFilter - add derivatives along radius and theta
  //
  
  AddFilterType::Pointer addFilter41 = AddFilterType::New();

  addFilter41->SetInput1( binaryFilter->GetOutput() );
  addFilter41->SetInput2( openingFilter->GetOutput() );

  addFilter41->Update();
  

  
  InternalITKImageType::Pointer addImage41 = addFilter41->GetOutput();
  printf("addImage41 size %d %d\n", addImage41->GetLargestPossibleRegion().GetSize()[0], addImage41->GetLargestPossibleRegion().GetSize()[1] );
  printf("addImage41 spacing %f %f\n", addImage41->GetSpacing()[0], addImage41->GetSpacing()[1] );
  printf("addImage41 origin %f %f\n", addImage41->GetOrigin()[0], addImage41->GetOrigin()[1] );
  
  // Show image
  sitk::Show( sitk::Image( addFilter41->GetOutput() ) );
  

  /*
  sitk::ImageFileWriter writer;
  writer.SetFileName( "C:\\Users\\yingying\\Desktop\\test181.png" );
  writer.Execute( sitkImageOut );
  */
  /*
  itkImage = dynamic_cast <InternalITKImageType*> ( sitkImageOut.GetITKBase() );
  outputOrigin = itkImage->GetOrigin();
  printf("origin %f %f", outputOrigin[0], outputOrigin[1] );
  */
  return EXIT_SUCCESS;
}






