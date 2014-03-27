
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "SimpleITK.h"

#include "itkVector.h"

#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"

#include "vnl/vnl_math.h"

#include "itkMedianImageFilter.h"
#include "itkChangeInformationImageFilter.h"
#include "itkRecursiveGaussianImageFilter.h"

#include "itkStatisticsImageFilter.h"
#include "itkSigmoidImageFilter.h"
#include "itkIntensityWindowingImageFilter.h"

#include "itkResampleImageFilter.h"
#include "itkTranslationTransform.h"
#include "itkLinearInterpolateImageFunction.h"
#include "itkCropImageFilter.h"
#include "itkAddImageFilter.h"

// define mathematical constant
const double PI = 3.14159265359;

// using SimpleITK for simple procedures like reading in images and showing images in imageJ
namespace sitk = itk::simple;

// define a few global typedefs that are repeatedly used to define images and classes
typedef itk::Image< float, 2 > SliceImageType;

// make theta resolution in polar coordinates half a degree
const double thetaResolution = PI/360;


int main( int argc, char * argv[] )
{
  // make sure enough arguments are passed into exe
  if( argc < 2 )
  {
    std::cerr << "Usage: " << std::endl;
    std::cerr << argv[0] << "USshoulder.exe (radius of MedianFilter) (parameters)"<< std::endl;
    return EXIT_FAILURE;
  }

							
  //
  // ImageFileReader
  //
  std::string filename = "c:\\miia\\yingyinw\\USshoulder\\shoulderSource\\yy_3_26_2014_36.png";
  std::string filterType = argv[1];

  typedef itk::ImageFileReader<SliceImageType> SliceReader;
  SliceReader::Pointer reader = SliceReader::New();
  reader->SetFileName(filename);
  reader->Update();
  
  SliceImageType::SizeType inputImageSize;
  inputImageSize[0] = reader->GetOutput()->GetLargestPossibleRegion().GetSize()[0];
  inputImageSize[1] = reader->GetOutput()->GetLargestPossibleRegion().GetSize()[1];
  
  printf("itkImage size %d %d\n", inputImageSize[0], inputImageSize[1]);
  printf("itkImage spacing %f %f\n", reader->GetOutput()->GetSpacing()[0], reader->GetOutput()->GetSpacing()[1] );
  printf("itkImage origin %f %f\n", reader->GetOutput()->GetOrigin()[0], reader->GetOutput()->GetOrigin()[1] );
  /*  
  // Show image
  sitk::Show( sitk::Image( reader->GetOutput() ) );
  */
  
  
  
  if (filterType.compare("crop") == 0)
  {
	int x1 = atoi(argv[2]);
	int x2 = atoi(argv[3]);
	int y1 = atoi(argv[4]);
	int y2 = atoi(argv[5]);

	printf("\nCrop corners x1:%d, x2:%d, y1:%d, y2:%d\n", x1, x2, y1, y2);

	// cropsize value is the amoutn removed on one side
	SliceImageType::SizeType cropSize;
	cropSize[0] = (inputImageSize[0]-(x2-x1))/2;
	cropSize[1] = (inputImageSize[1]-(y2-y1))/2;

	printf("\nCrop size lower:%d, upper:%d\n", cropSize[0], cropSize[1]);

	//
	// ResampleImageFilter 
	//
	typedef itk::ResampleImageFilter<SliceImageType,SliceImageType> ResampleFilterType;
	ResampleFilterType::Pointer ResampleFilter = ResampleFilterType::New();
	
	// translate
	typedef itk::TranslationTransform< double, 2 >  TransformType;
	TransformType::Pointer cropShift = TransformType::New();
	TransformType::OutputVectorType translation;
	translation[0] = -int((cropSize[0])-x1);
	translation[1] = -int((cropSize[1])-y1);
	cropShift->Translate( translation );
	
	// interpolator
	typedef itk::LinearInterpolateImageFunction<SliceImageType, double> InterpolatorType;
	InterpolatorType::Pointer interpolator = InterpolatorType::New();
			
	typedef itk::ResampleImageFilter<SliceImageType, SliceImageType> ResampleFilterType;
	ResampleFilterType::Pointer resampler = ResampleFilterType::New();
	resampler->SetTransform(cropShift);	
	resampler->SetInterpolator(interpolator);	
	resampler->SetSize(inputImageSize);
	resampler->SetDefaultPixelValue(0);
	resampler->SetInput(reader->GetOutput());
	resampler->Update();	
	
	printf("\nImage shifted\n");
	printf("Image translate by x:%f y:%f\n\n", translation[0], translation[1]);
	
	// Show image
	SliceImageType::Pointer shiftedImage = resampler->GetOutput();
	printf("shiftedImage size %d %d\n", shiftedImage->GetLargestPossibleRegion().GetSize()[0], shiftedImage->GetLargestPossibleRegion().GetSize()[1] );
	printf("shiftedImage spacing %f %f\n", shiftedImage->GetSpacing()[0], shiftedImage->GetSpacing()[1] );
	printf("shiftedImage origin %f %f\n", shiftedImage->GetOrigin()[0], shiftedImage->GetOrigin()[1] );
	
	sitk::Show( sitk::Image( resampler->GetOutput() ) );

	
	//
	// CropImageFilter
	//
	typedef itk::CropImageFilter <SliceImageType, SliceImageType> CropImageFilterType;
	CropImageFilterType::Pointer cropFilterTop = CropImageFilterType::New();

	cropFilterTop->SetBoundaryCropSize(cropSize);
	cropFilterTop->SetInput(resampler->GetOutput());

	printf("\nImage cropped\n");
	printf("Image cropped by x:%d y:%d on each side\n\n", cropSize[0]/2, cropSize[1]/2);

	// Show image
	SliceImageType::Pointer croppedImage = cropFilterTop->GetOutput();
	printf("croppedImage size %d %d\n", croppedImage->GetLargestPossibleRegion().GetSize()[0], croppedImage->GetLargestPossibleRegion().GetSize()[1] );
	printf("croppedImage spacing %f %f\n", croppedImage->GetSpacing()[0], croppedImage->GetSpacing()[1] );
	printf("croppedImage origin %f %f\n", croppedImage->GetOrigin()[0], croppedImage->GetOrigin()[1] );

	sitk::Show( sitk::Image( cropFilterTop->GetOutput() ) );
	// 
	// AddImageFilter
	//	
	
	

  }
  

  
	
  
  //
  // MedianImageFilter 
  // blurs input image to remove speckles and artifacts
  //
  if (filterType.compare("median") == 0)
  {
	typedef itk::MedianImageFilter< SliceImageType, SliceImageType >  MedianFilterType;
	MedianFilterType::Pointer medianFilter = MedianFilterType::New();  
	medianFilter->SetInput( reader->GetOutput() ); 
	medianFilter->SetRadius( atof( argv[2] ) );
	medianFilter->Update();

	SliceImageType::Pointer medianImage = medianFilter->GetOutput();
	printf("medianImage size %d %d\n", medianImage->GetLargestPossibleRegion().GetSize()[0], medianImage->GetLargestPossibleRegion().GetSize()[1] );
	printf("medianImage spacing %f %f\n", medianImage->GetSpacing()[0], medianImage->GetSpacing()[1] );
	printf("medianImage origin %f %f\n", medianImage->GetOrigin()[0], medianImage->GetOrigin()[1] );

	// Show image
	sitk::Show( sitk::Image( medianFilter->GetOutput() ) );
  }

  
  //
  // RecursiveGaussianImageFilter
  //  
  if (filterType.compare("gaussian") == 0)
  {
	typedef itk::RecursiveGaussianImageFilter< SliceImageType, SliceImageType > GaussianFilterType;
	GaussianFilterType::Pointer gaussianFilter = GaussianFilterType::New();

	gaussianFilter->SetSigma( atof(argv[2]) ); 
	gaussianFilter->SetDirection( 0 ); 
	gaussianFilter->SetNormalizeAcrossScale( true );

	gaussianFilter->SetInput( reader->GetOutput() );
	gaussianFilter->Update();
  
	SliceImageType::Pointer gaussianImage = gaussianFilter->GetOutput();
	printf("gaussianImage size %d %d\n", gaussianImage->GetLargestPossibleRegion().GetSize()[0], gaussianImage->GetLargestPossibleRegion().GetSize()[1] );
	printf("gaussianImage spacing %f %f\n", gaussianImage->GetSpacing()[0], gaussianImage->GetSpacing()[1] );
	printf("gaussianImage origin %f %f\n", gaussianImage->GetOrigin()[0], gaussianImage->GetOrigin()[1] );
  
	// Show image
	sitk::Show( sitk::Image( gaussianFilter->GetOutput() ) );
  }
  
  
  
  //
  // StatisticsImageFilter
  //
  /*
  typedef itk::StatisticsImageFilter< SliceImageType > StatisticsFilterType;
  StatisticsFilterType::Pointer statisticsFilter = StatisticsFilterType::New();
  statisticsFilter->SetInput( medianFilter->GetOutput() );
  statisticsFilter->Update();
  
  printf("Max: %f\nMin: %f\nMean: %f\nSD: %f\n\n", 
			statisticsFilter->GetMaximum(), statisticsFilter->GetMinimum(), 
			statisticsFilter->GetMean(), statisticsFilter->GetSigma());
  */			
			
  //
  // SigmoidImageFilter
  //
  
  if (filterType.compare("sigmoid") == 0)
  {
	typedef itk::SigmoidImageFilter< SliceImageType, SliceImageType > SigmoidFilterType;
	SigmoidFilterType::Pointer sigmoid = SigmoidFilterType::New();
	sigmoid->SetOutputMinimum( 0.0 );
	sigmoid->SetOutputMaximum( 1.0 );
	sigmoid->SetBeta( atof( argv[2] ) );
	sigmoid->SetAlpha( - atof( argv[3] ) );
  
	sigmoid->SetInput( reader->GetOutput() );
	sigmoid->Update();
  
  
	SliceImageType::Pointer sigmoidImage = sigmoid->GetOutput();
	printf("sigmoidImage size %d %d\n", sigmoidImage->GetLargestPossibleRegion().GetSize()[0], sigmoidImage->GetLargestPossibleRegion().GetSize()[1] );
	printf("sigmoidImage spacing %f %f\n", sigmoidImage->GetSpacing()[0], sigmoidImage->GetSpacing()[1] );
	printf("sigmoidImage origin %f %f\n", sigmoidImage->GetOrigin()[0], sigmoidImage->GetOrigin()[1] );
  
	// Show image
	sitk::Show( sitk::Image( sigmoid->GetOutput() ) );
  }
  
  //
  // IntensityWindowImageFilter
  //
  /*
  typedef itk::IntensityWindowingImageFilter< SliceImageType, SliceImageType >  IntensityFilterType;
  IntensityFilterType::Pointer intensityWindowing = IntensityFilterType::New();

  intensityWindowing->SetWindowMinimum( 0 );
  intensityWindowing->SetWindowMaximum( 1000 );

  intensityWindowing->SetOutputMinimum(   0.0 );
  intensityWindowing->SetOutputMaximum( 1000 ); // floats but in the range of chars.

  intensityWindowing->SetInput( fastMarching->GetOutput() );

  
  SliceImageType::Pointer intensityWindowingImage = intensityWindowing->GetOutput();
  printf("intensityWindowingImage size %d %d\n", intensityWindowingImage->GetLargestPossibleRegion().GetSize()[0], intensityWindowingImage->GetLargestPossibleRegion().GetSize()[1] );
  printf("intensityWindowingImage spacing %f %f\n", intensityWindowingImage->GetSpacing()[0], intensityWindowingImage->GetSpacing()[1] );
  printf("intensityWindowingImage origin %f %f\n", intensityWindowingImage->GetOrigin()[0], intensityWindowingImage->GetOrigin()[1] );
  
  // Show image
  sitk::Show( sitk::Image( intensityWindowing->GetOutput() ) );
  */

  /*
  sitk::ImageFileWriter writer;
  writer.SetFileName( "C:\\Users\\yingying\\Desktop\\test181.png" );
  writer.Execute( thresholder->GetOutput() );
  */
  /*
  itkImage = dynamic_cast <SliceImageType*> ( thresholder->GetOutput().GetITKBase() );
  outputOrigin = itkImage->GetOrigin();
  printf("origin %f %f", outputOrigin[0], outputOrigin[1] );
  */

  return EXIT_SUCCESS;
}






