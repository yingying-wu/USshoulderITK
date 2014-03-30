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

#include "vnl/vnl_math.h"

#include "itkMedianImageFilter.h"
#include "itkRecursiveGaussianImageFilter.h"
#include "itkStatisticsImageFilter.h"
#include "itkIntensityWindowingImageFilter.h"
#include "itkSigmoidImageFilter.h"

#include "itkResampleImageFilter.h"
#include "itkTranslationTransform.h"
#include "itkLinearInterpolateImageFunction.h"
#include "itkAddImageFilter.h"
#include "itkImageRegionIterator.h"
#include "itkRegionOfInterestImageFilter.h"
#include "itkImageDuplicator.h"

const double PI = 3.14159265359;
namespace sitk = itk::simple;
typedef itk::Image< float, 2 > SliceImageType;


int main( int argc, char * argv[] )
{
  if( argc < 2 )
    {
    std::cerr << "Usage: " << std::endl;
    std::cerr << argv[0] << " USshoulder.exe"<< std::endl;
    return EXIT_FAILURE;
    }

  // declare typedefs here
  	typedef itk::ImageFileReader<SliceImageType> SliceReader;	
	typedef itk::ImageFileWriter< SliceImageType > WriterType;

	typedef itk::MedianImageFilter< SliceImageType, SliceImageType >  MedianFilterType;
	typedef itk::RecursiveGaussianImageFilter< SliceImageType, SliceImageType > GaussianFilterType;
	typedef itk::SigmoidImageFilter< SliceImageType, SliceImageType > SigmoidFilterType;
	typedef itk::IntensityWindowingImageFilter< SliceImageType, SliceImageType >  IntensityFilterType;

	typedef itk::LinearInterpolateImageFunction<SliceImageType, double> InterpolatorType;
	typedef itk::ResampleImageFilter<SliceImageType, SliceImageType> ResampleFilterType;
	typedef itk::TranslationTransform< double, 2 >  TransformType;
	typedef itk::AddImageFilter < SliceImageType, SliceImageType > AddFilterType;
	typedef itk::RegionOfInterestImageFilter< SliceImageType, SliceImageType > ROIType;
	typedef itk::ImageDuplicator< SliceImageType > DuplicatorType;


	//-----------------------
	// Read and process image
	//-----------------------

	std::string filename = "c:\\miia\\yingyinw\\USshoulder\\shoulderSource\\yy_3_26_2014_36.png";
	std::string filterType = argv[1];
  
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
	//sitk::Show( sitk::Image( reader->GetOutput() ) );

	
	//---------------------------
	// Crop Original Screen Image
	
	// hard code screen borders
	int screenX1 = 192;
	int screenX2 = 588;
	int screenY1 = 31;
	int screenY2 = 411;

	SliceImageType::IndexType screenStart;
	screenStart[0] = screenX1;
	screenStart[1] = screenY1; 
	SliceImageType::SizeType screenSize;
	screenSize[0] = screenX2 - screenX1;
	screenSize[1] = screenY2 - screenY1; 
	SliceImageType::RegionType desiredRegionScreen( screenStart, screenSize );

	ROIType::Pointer ROIscreen = ROIType::New();
	ROIscreen->SetRegionOfInterest( desiredRegionScreen );
	ROIscreen->SetInput( inputImage );
	ROIscreen->Update();

	SliceImageType::Pointer croppedImage = SliceImageType::New();
	croppedImage = ROIscreen->GetOutput();

	//sitk::Show( sitk::Image( croppedImage ) );
	printf("\nExtracted image from screenshot\n");
	printf("Extracted image size %d %d\n", croppedImage->GetLargestPossibleRegion().GetSize()[0], croppedImage->GetLargestPossibleRegion().GetSize()[1]);

	SliceImageType::SizeType croppedImageSize;
	croppedImageSize[0] = croppedImage->GetLargestPossibleRegion().GetSize()[0];
	croppedImageSize[1] = croppedImage->GetLargestPossibleRegion().GetSize()[1];
	//printf("Cropped image size %d %d\n", croppedImageSize[0], croppedImageSize[1]);


	//--------------
	// Bright region
	//--------------

	//----------------------
	// Crop Top Bright Image 

	// hard code screen borders
	int topX1 = 115;
	int topX2 = 285;
	int topY1 = 1;
	int topY2 = 181;

	SliceImageType::IndexType topStart;
	topStart[0] = topX1;
	topStart[1] = topY1; 
	SliceImageType::SizeType topSize;
	topSize[0] = topX2 - topX1;
	topSize[1] = topY2 - topY1; 
	SliceImageType::RegionType desiredRegionTop( topStart, topSize );

	ROIType::Pointer ROITop = ROIType::New();
	ROITop->SetRegionOfInterest( desiredRegionTop );
	ROITop->SetInput( croppedImage );
	ROITop->Update();

	SliceImageType::Pointer topImage = SliceImageType::New();
	topImage = ROITop->GetOutput();

	//sitk::Show( sitk::Image( topImage ) );
	printf("\nTop image\n");
	printf("Top image size %d %d\n", topImage->GetLargestPossibleRegion().GetSize()[0], topImage->GetLargestPossibleRegion().GetSize()[1]);

	// translate	
	TransformType::Pointer topCropShift = TransformType::New();
	TransformType::OutputVectorType topTranslation;
	topTranslation[0] = 0;//-int( topX1 );
	topTranslation[1] = 0;//-int( topX2 );
	topCropShift->Translate( topTranslation );
	//printf("\nBottom section translation x:%f y:%f\n", translation2[0], translation2[1]);		

	// interpolator
	typedef itk::LinearInterpolateImageFunction<SliceImageType, double> InterpolatorType;
	InterpolatorType::Pointer interpolator = InterpolatorType::New();

	SliceImageType::PointType topOrigin;
	topOrigin[0] = screenX1;//-int( topTranslation[0] );
	topOrigin[1] = screenY1;//-int( topTranslation[1] );

	ResampleFilterType::Pointer resamplerTop = ResampleFilterType::New();
	resamplerTop->SetTransform( topCropShift );	
	resamplerTop->SetInterpolator( interpolator );	
	resamplerTop->SetSize( croppedImageSize );
	resamplerTop->SetOutputOrigin( topOrigin );
	resamplerTop->SetDefaultPixelValue( 0 );
	resamplerTop->SetInput( topImage );
	resamplerTop->Update();	

	SliceImageType::Pointer topImageShifted = SliceImageType::New();
	topImageShifted = resamplerTop->GetOutput();

	//sitk::Show( sitk::Image( topImageShifted ) );
	printf("\nTop image shifted\n");
	printf("Top image shifted size %d %d\n", topImageShifted->GetLargestPossibleRegion().GetSize()[0], topImageShifted->GetLargestPossibleRegion().GetSize()[1]);
		
	//-------------------------
	// Crop Bottom Bright Image
	
	//sitk::Show( sitk::Image( croppedImage ) );

	// hard code screen borders
	int botX1 = 98;
	int botX2 = 303;
	int botY1 = 181;
	int botY2 = 380;

	SliceImageType::IndexType botStart;
	botStart[0] = botX1;
	botStart[1] = botY1; 
	SliceImageType::SizeType botSize;
	botSize[0] = botX2 - botX1;
	botSize[1] = botY2 - botY1; 
	SliceImageType::RegionType desiredRegionBot( botStart, botSize );

	ROIType::Pointer ROIBot = ROIType::New();
	ROIBot->SetRegionOfInterest( desiredRegionBot );
	ROIBot->SetInput( croppedImage );
	ROIBot->Update();

	SliceImageType::Pointer botImage = SliceImageType::New();
	botImage = ROIBot->GetOutput();

	//sitk::Show( sitk::Image( botImage ) );
	printf("\nBottom image\n");
	printf("Bottom image size %d %d\n", botImage->GetLargestPossibleRegion().GetSize()[0], botImage->GetLargestPossibleRegion().GetSize()[1]);

	// translate	
	TransformType::Pointer botCropShift = TransformType::New();
	TransformType::OutputVectorType botTranslation;
	botTranslation[0] = 0;//-int( botX1 );
	botTranslation[1] = 0;//-int( botX2 );
	botCropShift->Translate( botTranslation );
	//printf("\nBottom section translation x:%f y:%f\n", translation2[0], translation2[1]);		

	SliceImageType::PointType botOrigin;
	botOrigin[0] = screenX1;//-int( botTranslation[0] );
	botOrigin[1] = screenY1;//-int( botTranslation[1] );

	ResampleFilterType::Pointer resamplerBot = ResampleFilterType::New();
	resamplerBot->SetTransform( botCropShift );	
	resamplerBot->SetInterpolator( interpolator );	
	resamplerBot->SetSize( croppedImageSize );
	resamplerBot->SetOutputOrigin( botOrigin );
	resamplerBot->SetDefaultPixelValue( 0 );
	resamplerBot->SetInput( botImage );
	resamplerBot->Update();	
	
	SliceImageType::Pointer botImageShifted = SliceImageType::New();
	botImageShifted = resamplerBot->GetOutput();

	//sitk::Show( sitk::Image( botImageShifted ) );
	printf("\nBottom image shifted\n");
	printf("Bottom image shifted size %d %d\n", botImageShifted->GetLargestPossibleRegion().GetSize()[0], botImageShifted->GetLargestPossibleRegion().GetSize()[1]);
	
	//----------------------
	// Create Middle Section

	AddFilterType::Pointer addFilter = AddFilterType::New();
	addFilter->SetInput1( topImageShifted );
	addFilter->SetInput2( botImageShifted );
	addFilter->Update();

	SliceImageType::Pointer brightImage = SliceImageType::New();
	brightImage = addFilter->GetOutput();

	sitk::Show( sitk::Image( brightImage ) );
	printf("\nTop and bottom merged\n");
	printf("Bright image cropped out\n");
	

	//-----------
	//Dark region
	//-----------
	
	//-----------------------
	// Duplicate croppedImage

    DuplicatorType::Pointer duplicator = DuplicatorType::New();
    duplicator->SetInputImage(croppedImage);
    duplicator->Update();

    SliceImageType::Pointer clonedImage = SliceImageType::New();
	clonedImage = duplicator->GetOutput();

	printf("\nCloned image\n");
	printf("Cloned image size %d %d\n", clonedImage->GetLargestPossibleRegion().GetSize()[0], clonedImage->GetLargestPossibleRegion().GetSize()[1]);

	//----------
	// Top Image
		
	itk::ImageRegionIterator<SliceImageType> topImageIterator( clonedImage,desiredRegionTop );
 
	while(!topImageIterator.IsAtEnd())
    {
		topImageIterator.Set( 0 ); 
		++topImageIterator;
	}

	//sitk::Show( sitk::Image( clonedImage ) );
	printf("\nDone iterating bottom\n");
	
	//-------------
	// Bottom Image

	itk::ImageRegionIterator<SliceImageType> botImageIterator( clonedImage,desiredRegionBot );
 
	while(!botImageIterator.IsAtEnd())
    {	
		botImageIterator.Set( 0 ); 
		++botImageIterator;
	}

	//sitk::Show( sitk::Image( clonedImage ) );			
	printf("\nDone iterating bottom\n");
	printf("\nDark image\n");
	printf("Dark image size %d %d\n", clonedImage->GetLargestPossibleRegion().GetSize()[0], clonedImage->GetLargestPossibleRegion().GetSize()[1]);
	
	printf("\nDark image cropped out.\n");
	

	//-------------------------------
	// Image processing bright region
	//-------------------------------

	SliceImageType::Pointer brightImageProcessed = SliceImageType::New();

	//-------
	// Median

	if (filterType.compare("median") == 0)
	{
	MedianFilterType::Pointer medianFilter = MedianFilterType::New();  
	medianFilter->SetInput( brightImage ); 
	medianFilter->SetRadius( atoi( argv[2] ) );
	medianFilter->Update();
  
  
	SliceImageType::Pointer medianImage = medianFilter->GetOutput();
	printf("medianImage size %d %d\n", medianImage->GetLargestPossibleRegion().GetSize()[0], medianImage->GetLargestPossibleRegion().GetSize()[1] );
	
	sitk::Show( sitk::Image( medianFilter->GetOutput() ) );	
	brightImageProcessed = medianFilter->GetOutput();

	}

	//----------------------
	// StatisticsImageFilter
	
	typedef itk::StatisticsImageFilter< SliceImageType > StatisticsFilterType;
	StatisticsFilterType::Pointer midStatisticsFilter = StatisticsFilterType::New();
	midStatisticsFilter->SetInput( brightImage );
	midStatisticsFilter->Update();

  
	printf("Max: %f\nMin: %f\nMean: %f\nSD: %f\n\n", 
			midStatisticsFilter->GetMaximum(), midStatisticsFilter->GetMinimum(), 
			midStatisticsFilter->GetMean(), midStatisticsFilter->GetSigma());

	//-----------------------------
	// RecursiveGaussianImageFilter

	if (filterType.compare("gaussian") == 0)
	{
	GaussianFilterType::Pointer gaussianFilterX = GaussianFilterType::New();

	gaussianFilterX->SetSigma( atof(argv[2]) ); 
	gaussianFilterX->SetDirection( 0 ); 
	gaussianFilterX->SetNormalizeAcrossScale( true );

	gaussianFilterX->SetInput( brightImage );
	gaussianFilterX->Update();

	GaussianFilterType::Pointer gaussianFilterY = GaussianFilterType::New();

	gaussianFilterY->SetSigma( atof(argv[2]) ); 
	gaussianFilterY->SetDirection( 1 ); 
	gaussianFilterY->SetNormalizeAcrossScale( true );

	gaussianFilterY->SetInput( gaussianFilterX->GetOutput() );
	gaussianFilterY->Update();

	printf("\nGaussian blur %f\n", atof(argv[2]) );
	sitk::Show( sitk::Image( gaussianFilterY->GetOutput() ) );
	brightImageProcessed = gaussianFilterY->GetOutput();
	}
  

	//-------------------
	// SigmoidImageFilter
  
	if (filterType.compare("sigmoid") == 0)
	{
	SigmoidFilterType::Pointer sigmoid = SigmoidFilterType::New();
	sigmoid->SetOutputMinimum( 0.0 );
	sigmoid->SetOutputMaximum( 1.0 );
	sigmoid->SetBeta( atof( argv[2] ) );
	sigmoid->SetAlpha( atof( argv[3] ) );
  
	sigmoid->SetInput( brightImage );
	sigmoid->Update();
  
	sitk::Show( sitk::Image( sigmoid->GetOutput() ) );
	brightImageProcessed = sigmoid->GetOutput();

	StatisticsFilterType::Pointer winStatisticsFilter = StatisticsFilterType::New();
	winStatisticsFilter->SetInput( sigmoid->GetOutput() );
	winStatisticsFilter->Update();
  
	printf("\nMidSection after windowing \nMax: %f\nMin: %f\nMean: %f\nSD: %f\n\n", 
			winStatisticsFilter->GetMaximum(), winStatisticsFilter->GetMinimum(), 
			winStatisticsFilter->GetMean(), winStatisticsFilter->GetSigma());

	}
  
	//---------------------------
	// IntensityWindowImageFilter

	if (filterType.compare("window") == 0)
	{
	IntensityFilterType::Pointer intensityWindowing = IntensityFilterType::New();

	intensityWindowing->SetWindowMinimum( atof( argv[2] ) );
	intensityWindowing->SetWindowMaximum( atof( argv[3] ) );

	intensityWindowing->SetOutputMinimum( 0 );
	intensityWindowing->SetOutputMaximum( 255 ); // floats but in the range of chars.

	intensityWindowing->SetInput( brightImage );

	sitk::Show( sitk::Image( intensityWindowing->GetOutput() ) );
	brightImageProcessed = intensityWindowing->GetOutput();
	}

	

  return EXIT_SUCCESS;
}






