
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
#include "itkConstantPadImageFilter.h"
#include "itkImageDuplicator.h"
#include "itkImageRegionIterator.h"

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

							
	//-----------------------
	// Read and process image
	//-----------------------

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
    
	//sitk::Show( sitk::Image( reader->GetOutput() ) );
  
	// hard code screen borders
	int screenX1 = 192;
	int screenX2 = 588;
	int screenY1 = 31;
	int screenY2 = 411;

	// cropsize value is the amount removed on one side
	SliceImageType::SizeType screenCropSize;
	screenCropSize[0] = (inputImageSize[0]-(screenX2-screenX1))/2;
	screenCropSize[1] = (inputImageSize[1]-(screenY2-screenY1))/2;

	//
	// ResampleImageFilter 
	//
	
	// translate
	typedef itk::TranslationTransform< double, 2 >  TransformType;
	TransformType::Pointer screenCropShift = TransformType::New();
	TransformType::OutputVectorType screenTranslation;
	screenTranslation[0] = -int((screenCropSize[0])-screenX1);
	screenTranslation[1] = -int((screenCropSize[1])-screenY1);
	screenCropShift->Translate( screenTranslation );
	
	// interpolator
	typedef itk::LinearInterpolateImageFunction<SliceImageType, double> InterpolatorType;
	InterpolatorType::Pointer interpolator = InterpolatorType::New();
			
	typedef itk::ResampleImageFilter<SliceImageType, SliceImageType> ResampleFilterType;
	ResampleFilterType::Pointer screenResampler = ResampleFilterType::New();
	screenResampler->SetTransform(screenCropShift);	
	screenResampler->SetInterpolator(interpolator);	
	screenResampler->SetSize(inputImageSize);
	screenResampler->SetDefaultPixelValue(0);
	screenResampler->SetInput(reader->GetOutput());
	screenResampler->Update();	
	
	//
	// CropImageFilter
	//
	typedef itk::CropImageFilter <SliceImageType, SliceImageType> CropImageFilterType;
	CropImageFilterType::Pointer screenCropFilter = CropImageFilterType::New();

	screenCropFilter->SetBoundaryCropSize(screenCropSize);
	screenCropFilter->SetInput(screenResampler->GetOutput());
 
	// for easy reference
	SliceImageType::Pointer inputImage = screenCropFilter->GetOutput();
	printf("\nImage cropped to ultrasound data\n");
	//sitk::Show( sitk::Image( inputImage ) );
	inputImageSize[0] = inputImage->GetLargestPossibleRegion().GetSize()[0];
	inputImageSize[1] = inputImage->GetLargestPossibleRegion().GetSize()[1];
	printf("itkImage size %d %d\n", inputImageSize[0], inputImageSize[1]);


	//-------------------------
	// Create image section top
	//-------------------------

	int topX1 = 115;
	int topX2 = 285;
	int topY1 = 1;
	int topY2 = 181;
	
	// cropsize value is the amoutn removed on one side
	SliceImageType::SizeType cropSize1;
	cropSize1[0] = (inputImageSize[0]-(topX2-topX1))/2;
	cropSize1[1] = (inputImageSize[1]-(topY2-topY1))/2;
	//printf("\n crop x:%d y:%d\n", cropSize1[0], cropSize1[1]);
	
	//
	// ResampleImageFilter 
	//

	// translate
	TransformType::Pointer cropShift1 = TransformType::New();
	TransformType::OutputVectorType translation1;
	translation1[0] = -int((cropSize1[0])-topX1);
	translation1[1] = -int((cropSize1[1])-topY1);
	cropShift1->Translate( translation1 );
	//printf("\n translation x:%f y:%f\n", translation1[0], translation1[1]);

	InterpolatorType::Pointer interpolator1 = InterpolatorType::New();
				
	ResampleFilterType::Pointer resampler1 = ResampleFilterType::New();
	resampler1->SetTransform(cropShift1);	
	resampler1->SetInterpolator(interpolator1);	
	resampler1->SetSize(inputImageSize);
	SliceImageType::PointType newOrigin;
	newOrigin[0] = screenX1-screenTranslation[0];
	newOrigin[1] = screenY1-screenTranslation[1];
	resampler1->SetOutputOrigin(newOrigin);
	resampler1->SetDefaultPixelValue(0);
	resampler1->SetInput(screenCropFilter->GetOutput());
	resampler1->Update();	
	//sitk::Show( sitk::Image( resampler1->GetOutput() ) );
	
	//
	// CropImageFilter
	//
	CropImageFilterType::Pointer cropFilter1 = CropImageFilterType::New();

	cropFilter1->SetBoundaryCropSize(cropSize1);
	cropFilter1->SetInput(resampler1->GetOutput());

	//
	// Pad to full image size
	//
	typedef itk::ConstantPadImageFilter <SliceImageType, SliceImageType> ConstantPadImageFilterType;
  
	ConstantPadImageFilterType::Pointer padFilter1 = ConstantPadImageFilterType::New();
	padFilter1->SetInput(cropFilter1->GetOutput());
	padFilter1->SetPadLowerBound(cropSize1);
	padFilter1->SetPadUpperBound(cropSize1);
	padFilter1->SetConstant(0);
	padFilter1->Update();

	SliceImageType::Pointer topSection = padFilter1->GetOutput();
	sitk::Show( sitk::Image( topSection ) );
	printf("\nTop section cropped out.\n");

	
	//-------------------------
	// Create image section bot
	//-------------------------

	int botX1 = 98;
	int botX2 = 303;
	int botY1 = 181;
	int botY2 = 379;
	
	// cropsize value is the amoutn removed on one side
	SliceImageType::SizeType cropSize2;
	cropSize2[0] = (inputImageSize[0]-(botX2-botX1))/2;
	cropSize2[1] = (inputImageSize[1]-(botY2-botY1))/2;
	
	//
	// ResampleImageFilter 
	//

	// translate	
	TransformType::Pointer cropShift2 = TransformType::New();
	TransformType::OutputVectorType translation2;
	translation2[0] = -int((cropSize2[0])-botX1);
	translation2[1] = -int((cropSize2[1])-botY1);
	cropShift2->Translate( translation2 );
	printf("\n translation x:%f y:%f\n", translation2[0], translation2[1]);
				
	ResampleFilterType::Pointer resampler2 = ResampleFilterType::New();
	resampler2->SetTransform(cropShift2);	
	resampler2->SetInterpolator(interpolator);	
	resampler2->SetSize(inputImageSize);
	resampler2->SetOutputOrigin(newOrigin);
	resampler2->SetDefaultPixelValue(0);
	resampler2->SetInput(inputImage);
	resampler2->Update();	
	
	//
	// CropImageFilter
	//
	CropImageFilterType::Pointer cropFilter2 = CropImageFilterType::New();
	cropFilter2->SetBoundaryCropSize(cropSize2);
	cropFilter2->SetInput(resampler2->GetOutput());

	//
	// Pad to full image size
	//  
	ConstantPadImageFilterType::Pointer padFilter2 = ConstantPadImageFilterType::New();
	padFilter2->SetInput(cropFilter2->GetOutput());
	padFilter2->SetPadLowerBound(cropSize2);
	padFilter2->SetPadUpperBound(cropSize2);
	padFilter2->SetConstant(0);
	padFilter2->Update();

	SliceImageType::Pointer botSection = padFilter2->GetOutput();
	//sitk::Show( sitk::Image( botSection ) );
	printf("\nBotton section cropped out.\n");

	
	//---------------------------
	// Create image section sides
	//---------------------------
	 
	typedef itk::ImageDuplicator< SliceImageType > DuplicatorType;
	DuplicatorType::Pointer duplicator = DuplicatorType::New();
	duplicator->SetInputImage(inputImage);
	duplicator->Update();
	SliceImageType::Pointer sideSection = duplicator->GetOutput();

	//
	// Set crop region to black 0
	//
	SliceImageType::RegionType region;
	SliceImageType::IndexType start;
	start[0] = topX1;
	start[1] = topY1;
 
	SliceImageType::SizeType size;
	size[0] = topX2-topX1;
	size[1] = topY2-topY1;
 
	region.SetSize(size);
	region.SetIndex(start);

	itk::ImageRegionIterator<SliceImageType> imageIterator1(sideSection,region);
 
	while(!imageIterator1.IsAtEnd())
    {
		SliceImageType::PixelType pixel = 0;		
		imageIterator1.Set(pixel); 
		++imageIterator1;
	}

	//sitk::Show( sitk::Image(sideSection) );
	
	start[0] = botX1;
	start[1] = botY1;
 	size[0] = botX2-botX1;
	size[1] = botY2-botY1;
	region.SetSize(size);
	region.SetIndex(start);

	itk::ImageRegionIterator<SliceImageType> imageIterator2(sideSection,region);
 
	while(!imageIterator2.IsAtEnd())
    {
		SliceImageType::PixelType pixel = 0;		
		imageIterator2.Set(pixel); 
		++imageIterator1;
	}

	//sitk::Show( sitk::Image(sideSection) );	

	printf("\nSides cropped out.\n");



  

  
	
  
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






