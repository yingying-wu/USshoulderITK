
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

#include "itkImageRegionIterator.h"

// define mathematical constant
const double PI = 3.14159265359;

// using SimpleITK for simple procedures like reading in images and showing images in imageJ
namespace sitk = itk::simple;

// define a few global typedefs that are repeatedly used to define images and classes
typedef itk::Image< float, 2 > SliceImageType;

// make theta resolution in polar coordinates half a degree
const double thetaResolution = PI/360;


int main(int argc, char * argv[] )
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
	
	SliceImageType::Pointer inputImage = SliceImageType::New();
	inputImage = reader->GetOutput();
	SliceImageType::SizeType inputImageSize;
	inputImageSize[0] = inputImage->GetLargestPossibleRegion().GetSize()[0];
	inputImageSize[1] = inputImage->GetLargestPossibleRegion().GetSize()[1];
  
	printf("Original image size %d %d\n", inputImageSize[0], inputImageSize[1]);
	printf("Image spacing %f %f\n", reader->GetOutput()->GetSpacing()[0], reader->GetOutput()->GetSpacing()[1] );
    		
	//sitk::Show( sitk::Image( inputImage ) );

	// interpolator - the same used for all 
	typedef itk::LinearInterpolateImageFunction<SliceImageType, double> InterpolatorType;
	InterpolatorType::Pointer interpolator = InterpolatorType::New();
	  
	// hard code screen borders
	int screenX1 = 192;
	int screenX2 = 588;
	int screenY1 = 31;
	int screenY2 = 411;

	// cropsize value is the amount removed on one side
	SliceImageType::SizeType screenCropSize;
	screenCropSize[0] = (inputImageSize[0]-(screenX2-screenX1))/2;
	screenCropSize[1] = (inputImageSize[1]-(screenY2-screenY1))/2;
	//printf("Screen crop %d %d\n", screenCropSize[0], screenCropSize[1]);

	//
	// ResampleImageFilter 
	//
	
	typedef itk::TranslationTransform< double, 2 >  TransformType;
	TransformType::Pointer screenCropShift = TransformType::New();
	TransformType::OutputVectorType screenTranslation;
	screenTranslation[0] = -int( (screenCropSize[0])-screenX1 );
	screenTranslation[1] = -int( (screenCropSize[1])-screenY1 ); 
	screenCropShift->Translate( screenTranslation );
	//printf("Screen translation %f %f\n", screenTranslation[0], screenTranslation[1]);
				
	typedef itk::ResampleImageFilter<SliceImageType, SliceImageType> ResampleFilterType;
	ResampleFilterType::Pointer screenResampler = ResampleFilterType::New();
	screenResampler->SetTransform( screenCropShift );	
	screenResampler->SetInterpolator( interpolator );	
	screenResampler->SetSize( inputImageSize );
	screenResampler->SetDefaultPixelValue( 0 );
	screenResampler->SetInput( inputImage );
	screenResampler->Update();	
	//sitk::Show( sitk::Image( screenResampler->GetOutput() ) );
	
	//
	// CropImageFilter
	//
	typedef itk::CropImageFilter <SliceImageType, SliceImageType> CropImageFilterType;
	CropImageFilterType::Pointer screenCropFilter = CropImageFilterType::New();

	screenCropFilter->SetBoundaryCropSize( screenCropSize );
	screenCropFilter->SetInput( screenResampler->GetOutput() );
	screenCropFilter->Update();

	// change pointer to cropped image. update image size
	SliceImageType::Pointer croppedImage = SliceImageType::New();
	croppedImage = screenCropFilter->GetOutput();
	//sitk::Show( sitk::Image( croppedImage ) );
	SliceImageType::SizeType croppedImageSize;
	croppedImageSize[0] = croppedImage->GetLargestPossibleRegion().GetSize()[0];
	croppedImageSize[1] = croppedImage->GetLargestPossibleRegion().GetSize()[1];
	//printf("Cropped image size %d %d\n", croppedImageSize[0], croppedImageSize[1]);


	printf("\nImage cropped to ultrasound data\n");
	// origin offset after cropping
	SliceImageType::PointType newOrigin;
	newOrigin[0] = screenX1-screenTranslation[0];
	newOrigin[1] = screenY1-screenTranslation[1];

	//-------------------------
	// Create image section top
	//-------------------------


	int topX1 = 115;
	int topX2 = 285;
	int topY1 = 1;
	int topY2 = 181;
	
	// cropsize value is the amount removed on one side
	SliceImageType::SizeType cropSize1;
	cropSize1[0] = ( croppedImageSize[0]-(topX2-topX1) )/2;
	cropSize1[1] = ( croppedImageSize[1]-(topY2-topY1) )/2;
	//printf("\nTop section crop %d %d\n", cropSize1[0], cropSize1[1]);
	
	//
	// ResampleImageFilter 
	//
	
	// translate
	TransformType::Pointer cropShift1 = TransformType::New();
	TransformType::OutputVectorType translation1;
	translation1[0] = -int( (cropSize1[0])-topX1 );
	translation1[1] = -int( (cropSize1[1])-topY1 );
	cropShift1->Translate( translation1 );
	//printf("\nTop section translation x:%f y:%f\n", translation1[0], translation1[1]);
	
				
	ResampleFilterType::Pointer resampler1 = ResampleFilterType::New();
	resampler1->SetTransform( cropShift1 );	
	resampler1->SetInterpolator( interpolator );	
	resampler1->SetSize( croppedImageSize );
	resampler1->SetOutputOrigin( newOrigin );
	resampler1->SetDefaultPixelValue( 0 );
	resampler1->SetInput( croppedImage );
	resampler1->Update();	
	//printf("resampler1 %d %d", resampler1->GetOutput()->GetLargestPossibleRegion().GetSize()[0], resampler1->GetOutput()->GetLargestPossibleRegion().GetSize()[1]);
	//sitk::Show( sitk::Image( resampler1->GetOutput() ) );
	

	//
	// CropImageFilter
	//
	CropImageFilterType::Pointer cropFilter1 = CropImageFilterType::New();
	cropFilter1->SetBoundaryCropSize( cropSize1 );
	cropFilter1->SetInput( resampler1->GetOutput() );
	//sitk::Show( sitk::Image( cropFilter1->GetOutput() ));
		
	TransformType::Pointer cropShiftInv1 = TransformType::New();
	TransformType::OutputVectorType translationInv1;
	translationInv1[0] = -translation1[0];
	translationInv1[1] = -translation1[1];
	cropShiftInv1->Translate( translationInv1 );

	ResampleFilterType::Pointer resamplerInvert1 = ResampleFilterType::New();
	resamplerInvert1->SetTransform( cropShiftInv1 );	
	resamplerInvert1->SetInterpolator( interpolator );	
	resamplerInvert1->SetSize( croppedImageSize );
	resamplerInvert1->SetOutputOrigin( newOrigin );
	resamplerInvert1->SetDefaultPixelValue( 0 );
	resamplerInvert1->SetInput( cropFilter1->GetOutput() );
	resamplerInvert1->Update();	
	//sitk::Show( sitk::Image( resamplerInvert1->GetOutput() ));
	

	SliceImageType::Pointer topSection = SliceImageType::New();
	topSection = resamplerInvert1->GetOutput();
	//sitk::Show( sitk::Image( topSection ) );
	printf("\nTop section cropped\n");

	
	//-------------------------
	// Create image section bot
	//-------------------------
	
	/*DuplicatorType::Pointer duplicatorBot = DuplicatorType::New();
	duplicatorBot->SetInputImage( croppedImage );
	duplicatorBot->Update();
	SliceImageType::Pointer botCopy = duplicatorBot->GetOutput();
	printf("\nbotCopy %f %f\n", botCopy->GetLargestPossibleRegion().GetSize()[0], botCopy->GetLargestPossibleRegion().GetSize()[1]);
	printf("\nDuplicated cropped image\n");
	sitk::Show( sitk::Image( botCopy ));
	*/
	int botX1 = 98;
	int botX2 = 303;
	int botY1 = 181;
	int botY2 = 380;
	
	// cropsize value is the amount removed on one side
	SliceImageType::SizeType cropSize2;
	cropSize2[0] = ( croppedImageSize[0]-(botX2-botX1) )/2;
	cropSize2[1] = ( croppedImageSize[1]-(botY2-botY1) )/2;
	
	//
	// ResampleImageFilter 
	//

	// translate	
	TransformType::Pointer cropShift2 = TransformType::New();
	TransformType::OutputVectorType translation2;
	translation2[0] = -int( (cropSize2[0])-botX1 );
	translation2[1] = -int( (cropSize2[1])-botY1 );
	cropShift2->Translate( translation2 );
	//printf("\nBottom section translation x:%f y:%f\n", translation2[0], translation2[1]);
				

	ResampleFilterType::Pointer resampler2 = ResampleFilterType::New();
	resampler2->SetTransform( cropShift2 );	
	resampler2->SetInterpolator( interpolator );	
	resampler2->SetSize( croppedImageSize );
	resampler2->SetOutputOrigin( newOrigin );
	resampler2->SetDefaultPixelValue( 0 );
	resampler2->SetInput( croppedImage );
	resampler2->Update();	
	

	//
	// CropImageFilter
	//
	CropImageFilterType::Pointer cropFilter2 = CropImageFilterType::New();
	cropFilter2->SetBoundaryCropSize( cropSize2 );
	cropFilter2->SetInput( resampler2->GetOutput() );

	TransformType::Pointer cropShiftInv2 = TransformType::New();
	TransformType::OutputVectorType translationInv2;
	translationInv2[0] = -translation2[0];
	translationInv2[1] = -translation2[1];
	cropShiftInv2->Translate( translationInv2 );


	ResampleFilterType::Pointer resamplerInvert2 = ResampleFilterType::New();
	resamplerInvert2->SetTransform( cropShiftInv2 );	
	resamplerInvert2->SetInterpolator( interpolator );	
	resamplerInvert2->SetSize( croppedImageSize );
	resamplerInvert2->SetOutputOrigin( newOrigin );
	resamplerInvert2->SetDefaultPixelValue( 0 );
	resamplerInvert2->SetInput( cropFilter2->GetOutput() );
	resamplerInvert2->Update();	
	//sitk::Show( sitk::Image( resamplerInvert2->GetOutput() ));
	

	SliceImageType::Pointer botSection = SliceImageType::New();
	botSection = resamplerInvert2->GetOutput();
	//sitk::Show( sitk::Image( botSection ) );
	printf("\nBotton section cropped\n");

	//
	// Add images
	//
	typedef itk::AddImageFilter < SliceImageType, SliceImageType > AddFilterType;
	AddFilterType::Pointer addFilter = AddFilterType::New();

	addFilter->SetInput1( botSection );
	addFilter->SetInput2( topSection );

	SliceImageType::Pointer midSection = SliceImageType::New();
	midSection = addFilter->GetOutput();
	sitk::Show( sitk::Image( midSection ) );
	printf("\nTop and bottom merged\n");

	
	//---------------------------
	// Create image section sides
	//---------------------------

	/*	
	//
	// Set crop region to black 0
	//
	SliceImageType::RegionTSype region;
	SliceImageType::IndexType start;
	start[0] = topX1+newOrigin[0];
	start[1] = topY1+newOrigin[1];
 
	SliceImageType::SizeType size;
	size[0] = topX2-topX1;
	size[1] = topY2-topY1;
 
	region.SetSize( size );
	region.SetIndex( start );

	itk::ImageRegionIterator<SliceImageType> imageIterator1(croppedImage,region);
 
	while(!imageIterator1.IsAtEnd())
    {
		SliceImageType::PixelType pixel = 0;		
		imageIterator1.Set( pixel ); 
		++imageIterator1;
	}

	printf("\nDone iterating\n");
	sitk::Show( sitk::Image(croppedImage) );
	

	start[0] = botX1+newOrigin[0];
	start[1] = botY1+newOrigin[1];
 	size[0] = botX2-botX1;
	size[1] = botY2-botY1;
	region.SetSize( size );
	region.SetIndex( start );

	itk::ImageRegionIterator<SliceImageType> imageIterator2(sideSection,region);
 
	while(!imageIterator2.IsAtEnd())
    {
		SliceImageType::PixelType pixel = 0;		
		imageIterator2.Set( pixel ); 
		++imageIterator1;
	}

	sitk::Show( sitk::Image(sideSection) );	

	printf("\nSides cropped out.\n");
	*/


	//
	// StatisticsImageFilter
	//
	typedef itk::StatisticsImageFilter< SliceImageType > StatisticsFilterType;
	StatisticsFilterType::Pointer midStatisticsFilter = StatisticsFilterType::New();
	midStatisticsFilter->SetInput( midSection );
	midStatisticsFilter->Update();
  
	printf("Max: %f\nMin: %f\nMean: %f\nSD: %f\n\n", 
			midStatisticsFilter->GetMaximum(), midStatisticsFilter->GetMinimum(), 
			midStatisticsFilter->GetMean(), midStatisticsFilter->GetSigma());
		

  

  
	
  
	//
	// MedianImageFilter 
	// blurs input image to remove speckles and artifacts
	//
	if (filterType.compare("median") == 0)
	{
	typedef itk::MedianImageFilter< SliceImageType, SliceImageType >  MedianFilterType;
	MedianFilterType::Pointer medianFilter = MedianFilterType::New();  
	medianFilter->SetInput( midSection ); 
	medianFilter->SetRadius( atof( argv[2] ) );
	medianFilter->Update();

	sitk::Show( sitk::Image( medianFilter->GetOutput() ) );
	}

  
	//
	// RecursiveGaussianImageFilter
	//  
	if (filterType.compare("gaussian") == 0)
	{
	typedef itk::RecursiveGaussianImageFilter< SliceImageType, SliceImageType > GaussianFilterType;
	GaussianFilterType::Pointer gaussianFilterX = GaussianFilterType::New();

	gaussianFilterX->SetSigma( atof(argv[2]) ); 
	gaussianFilterX->SetDirection( 0 ); 
	gaussianFilterX->SetNormalizeAcrossScale( true );

	gaussianFilterX->SetInput( midSection );
	gaussianFilterX->Update();

	GaussianFilterType::Pointer gaussianFilterY = GaussianFilterType::New();

	gaussianFilterY->SetSigma( atof(argv[2]) ); 
	gaussianFilterY->SetDirection( 1 ); 
	gaussianFilterY->SetNormalizeAcrossScale( true );

	gaussianFilterY->SetInput( gaussianFilterX->GetOutput() );
	gaussianFilterY->Update();

	printf("\nGaussian blur %f\n", atof(argv[2]) );
	sitk::Show( sitk::Image( gaussianFilterY->GetOutput() ) );
	}
  

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
	sigmoid->SetAlpha( -atof( argv[3] ) );
  
	sigmoid->SetInput( midSection );
	sigmoid->Update();
  
	sitk::Show( sitk::Image( sigmoid->GetOutput() ) );

	StatisticsFilterType::Pointer winStatisticsFilter = StatisticsFilterType::New();
	winStatisticsFilter->SetInput( sigmoid->GetOutput() );
	winStatisticsFilter->Update();
  
	printf("\nMidSection after windowing \nMax: %f\nMin: %f\nMean: %f\nSD: %f\n\n", 
			winStatisticsFilter->GetMaximum(), winStatisticsFilter->GetMinimum(), 
			winStatisticsFilter->GetMean(), winStatisticsFilter->GetSigma());

	}
  
	//
	// IntensityWindowImageFilter
	//
	if (filterType.compare("window") == 0)
	{
	typedef itk::IntensityWindowingImageFilter< SliceImageType, SliceImageType >  IntensityFilterType;
	IntensityFilterType::Pointer intensityWindowing = IntensityFilterType::New();

	intensityWindowing->SetWindowMinimum( atof( argv[2] ) );
	intensityWindowing->SetWindowMaximum( atof( argv[3] ) );

	intensityWindowing->SetOutputMinimum( 0 );
	intensityWindowing->SetOutputMaximum( 255 ); // floats but in the range of chars.

	intensityWindowing->SetInput( midSection );

	sitk::Show( sitk::Image( intensityWindowing->GetOutput() ) );
	}

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






