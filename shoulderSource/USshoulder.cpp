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

#include "itkGaussianOperator.h"
#include "itkNeighborhoodIterator.h"
#include "itkNeighborhoodInnerProduct.h"
#include "itkLaplacianOperator.h"

const double PI = 3.14159265359;
namespace sitk = itk::simple;
typedef itk::Image< float, 2 > SliceImageType;


int main( int argc, char * argv[] )
{

  // declare typedefs here
  	typedef itk::ImageFileReader<SliceImageType> SliceReader;	
	typedef itk::ImageFileWriter< SliceImageType > WriterType;

	//typedef itk::MedianImageFilter< SliceImageType, SliceImageType >  MedianFilterType;
	//typedef itk::RecursiveGaussianImageFilter< SliceImageType, SliceImageType > GaussianFilterType;
	typedef itk::SigmoidImageFilter< SliceImageType, SliceImageType > SigmoidFilterType;
	typedef itk::IntensityWindowingImageFilter< SliceImageType, SliceImageType >  IntensityFilterType;

	typedef itk::LinearInterpolateImageFunction<SliceImageType, double> InterpolatorType;
	typedef itk::ResampleImageFilter<SliceImageType, SliceImageType> ResampleFilterType;
	typedef itk::TranslationTransform< double, 2 >  TransformType;
	typedef itk::AddImageFilter < SliceImageType, SliceImageType > AddFilterType;
	typedef itk::RegionOfInterestImageFilter< SliceImageType, SliceImageType > ROIType;
	typedef itk::ImageDuplicator< SliceImageType > DuplicatorType;
	
	typedef itk::GaussianOperator< SliceImageType::PixelType, 2> GaussianOperatorType;
	typedef itk::ImageRegionIterator<SliceImageType> ImageIteratorType;
	typedef itk::NeighborhoodIterator<SliceImageType> NeighborhoodIteratorType;
	typedef itk::NeighborhoodInnerProduct< SliceImageType> InnerProductType; 

	// one-off functions that keep being reused
	typedef itk::StatisticsImageFilter< SliceImageType > StatisticsFilterType;
	StatisticsFilterType::Pointer statisticsFilter = StatisticsFilterType::New();

	typedef itk::LinearInterpolateImageFunction<SliceImageType, double> InterpolatorType;
	InterpolatorType::Pointer interpolator = InterpolatorType::New();


	//-----------------------
	// Read and process image
	//-----------------------

	std::string filename = "c:\\miia\\yingyinw\\USshoulder\\shoulderSource\\yy_3_26_2014_36.png";
  
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

	//----------------------------
	// hard coded crop imformation
	//----------------------------

	int topX1 = 114;
	int topX2 = 286;
	int topY1 = 0;
	int topY2 = 180;

	SliceImageType::IndexType topStart;
	topStart[0] = topX1;
	topStart[1] = topY1; 
	SliceImageType::SizeType topSize;
	topSize[0] = topX2 - topX1;
	topSize[1] = topY2 - topY1; 
	SliceImageType::RegionType desiredRegionTop( topStart, topSize );

	int botX1 = 97;
	int botX2 = 302;
	int botY1 = 180;
	int botY2 = 380;

	SliceImageType::IndexType botStart;
	botStart[0] = botX1;
	botStart[1] = botY1; 
	SliceImageType::SizeType botSize;
	botSize[0] = botX2 - botX1;
	botSize[1] = botY2 - botY1; 
	SliceImageType::RegionType desiredRegionBot( botStart, botSize );


	//--------------
	// Bright region
	//--------------

	//----------------------
	// Crop Top Bright Image 

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
	botTranslation[0] = 0;
	botTranslation[1] = 0;
	botCropShift->Translate( botTranslation );
	//printf("\nBottom section translation x:%f y:%f\n", translation2[0], translation2[1]);		

	SliceImageType::PointType botOrigin;
	botOrigin[0] = screenX1;
	botOrigin[1] = screenY1;

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

	//sitk::Show( sitk::Image( brightImage ), "brightImage" );
	printf("\nTop and bottom merged\n");
	printf("Bright image cropped out\n");
	

	//-------------------------------
	// Image processing bright region
	//-------------------------------

	SliceImageType::Pointer brightImageProcessed = SliceImageType::New();
		
	statisticsFilter->SetInput( brightImage );
	statisticsFilter->Update();

	printf("\nBright image statistics\n");
	printf("Max: %f\nMin: %f\nMean: %f\nSD: %f\n\n", 
			statisticsFilter->GetMaximum(), statisticsFilter->GetMinimum(), 
			statisticsFilter->GetMean(), statisticsFilter->GetSigma());

	//-------------------
	// SigmoidImageFilter
  
	SigmoidFilterType::Pointer sigmoidTop = SigmoidFilterType::New();
	sigmoidTop->SetOutputMinimum( 0 );
	sigmoidTop->SetOutputMaximum( 255 );
	sigmoidTop->SetBeta( 110 );
	sigmoidTop->SetAlpha( 50 );
  
	sigmoidTop->SetInput( brightImage );
	sigmoidTop->Update();
  
	SliceImageType::Pointer sigmoidImageTop = SliceImageType::New();
	sigmoidImageTop = sigmoidTop->GetOutput();
	//sitk::Show( sitk::Image( sigmoidImageTop ), "sigmoidImageTop");

	statisticsFilter->SetInput( sigmoidImageTop );
	statisticsFilter->Update();
  
	printf("\nBright image statistics after sigmoid\n");
	printf("Max: %f\nMin: %f\nMean: %f\nSD: %f\n\n", 
			statisticsFilter->GetMaximum(), statisticsFilter->GetMinimum(), 
			statisticsFilter->GetMean(), statisticsFilter->GetSigma());


	//-----------
	//Dark region
	//-----------
	
	//-----------------------
	// Duplicate croppedImage

    DuplicatorType::Pointer duplicator = DuplicatorType::New();
    duplicator->SetInputImage(croppedImage);
    duplicator->Update();

    SliceImageType::Pointer darkImage = SliceImageType::New();
	darkImage = duplicator->GetOutput();

	printf("\nCloned image\n");
	printf("Cloned image size %d %d\n", darkImage->GetLargestPossibleRegion().GetSize()[0], darkImage->GetLargestPossibleRegion().GetSize()[1]);

	//----------
	// Top Image
		
	itk::ImageRegionIterator<SliceImageType> topImageIterator( darkImage,desiredRegionTop );
 
	while(!topImageIterator.IsAtEnd())
    {
		topImageIterator.Set( 0 ); 
		++topImageIterator;
	}

	//sitk::Show( sitk::Image( darkImage ) );
	printf("\nDone iterating bottom\n");
	
	//-------------
	// Bottom Image

	itk::ImageRegionIterator<SliceImageType> botImageIterator( darkImage,desiredRegionBot );
 
	while(!botImageIterator.IsAtEnd())
    {	
		botImageIterator.Set( 0 ); 
		++botImageIterator;
	}


	//sitk::Show( sitk::Image( darkImage ), "darkImage" );			
	printf("\nDone iterating bottom\n");
	printf("\nDark image\n");
	printf("Dark image size %d %d\n", darkImage->GetLargestPossibleRegion().GetSize()[0], darkImage->GetLargestPossibleRegion().GetSize()[1]);
	
	printf("\nDark image cropped out.\n");

	//-----------------------------
	// Image processing dark region
	//-----------------------------

	SliceImageType::Pointer darkImageProcessed = SliceImageType::New();

	statisticsFilter->SetInput( darkImage );
	statisticsFilter->Update();

	printf("\nDark image statistics\n");
	printf("Max: %f\nMin: %f\nMean: %f\nSD: %f\n\n", 
			statisticsFilter->GetMaximum(), statisticsFilter->GetMinimum(), 
			statisticsFilter->GetMean(), statisticsFilter->GetSigma());

	//-------------------
	// SigmoidImageFilter
  
	SigmoidFilterType::Pointer sigmoidBot = SigmoidFilterType::New();
	sigmoidBot->SetOutputMinimum( 0 );
	sigmoidBot->SetOutputMaximum( 255 );
	sigmoidBot->SetBeta( 50 );
	sigmoidBot->SetAlpha( 20 );
  
	sigmoidBot->SetInput( darkImage );
	sigmoidBot->Update();
  
	SliceImageType::Pointer sigmoidImageBot = SliceImageType::New();
	sigmoidImageBot = sigmoidBot->GetOutput();
	//sitk::Show( sitk::Image( sigmoidImageBot ), "sigmoidImageBot");

	statisticsFilter->SetInput( sigmoidImageBot );
	statisticsFilter->Update();
  
	printf("\nDark image statistics after sigmoid\n");
	printf("Max: %f\nMin: %f\nMean: %f\nSD: %f\n\n", 
			statisticsFilter->GetMaximum(), statisticsFilter->GetMinimum(), 
			statisticsFilter->GetMean(), statisticsFilter->GetSigma());


	//-------------------
	// Create Final Image
	//-------------------

	AddFilterType::Pointer addFilterFinal = AddFilterType::New();
	addFilterFinal->SetInput1( sigmoidImageTop );
	addFilterFinal->SetInput2( sigmoidImageBot );
	addFilterFinal->Update();

	SliceImageType::Pointer mergedImage = SliceImageType::New();
	mergedImage = addFilterFinal->GetOutput();

	//sitk::Show( sitk::Image( mergedImage ), "mergedImage" );
	printf("Images merged\n");

	//--------------------
	// windowing threshold

	IntensityFilterType::Pointer windowBot = IntensityFilterType::New();
	windowBot->SetOutputMinimum( 0 );
	windowBot->SetOutputMaximum( 255 );
	windowBot->SetWindowMinimum( 0 );
	windowBot->SetWindowMaximum(255);
	windowBot->SetInput( mergedImage );
	windowBot->Update();
	
	SliceImageType::Pointer windowImage = SliceImageType::New();
	windowImage = windowBot->GetOutput();

	//sitk::Show( sitk::Image( windowImage ), "windowImage");
	printf("Images thresholded\n");
	
	//-----------------
	// get rid of lines
	//-----------------

	duplicator->SetInputImage(windowImage);
    duplicator->Update();
    SliceImageType::Pointer windowImageOut = SliceImageType::New();
	windowImageOut = duplicator->GetOutput();
	
	int width = 5;
	int space = 3; 
	double gaussianRadius = 6;

	GaussianOperatorType gaussianOperatorH;
	gaussianOperatorH.SetDirection(0); // Create the operator for the X axis derivative
	gaussianOperatorH.CreateToRadius( gaussianRadius );
	gaussianOperatorH.SetVariance( 3*gaussianRadius );
	GaussianOperatorType gaussianOperatorV;
	gaussianOperatorV.SetDirection(1); // Create the operator for the X axis derivative
	gaussianOperatorV.CreateToRadius( gaussianRadius );
	gaussianOperatorV.SetVariance( 3*gaussianRadius );
	InnerProductType IP;

	//v1 top left
	SliceImageType::IndexType v1Start;
	v1Start[0] = topX1 - width;
	v1Start[1] = topY1; 
	SliceImageType::SizeType v1Size;
	v1Size[0] = 2*width;
	v1Size[1] = (topY2 - topY1) + space;
	SliceImageType::RegionType v1( v1Start, v1Size );	

	ImageIteratorType itv1( windowImageOut, v1 );
	NeighborhoodIteratorType outv1( gaussianOperatorH.GetRadius(), windowImage, v1 );
	
	outv1.GoToBegin();
	for (itv1.GoToBegin(); !itv1.IsAtEnd(); ++itv1, ++outv1)
	{
		itv1.Set( IP( outv1, gaussianOperatorH) );
	}
	//sitk::Show( sitk::Image( windowImageOut ), "windowImageOut");

	//v2 top right
	SliceImageType::IndexType v2Start;
	v2Start[0] = topX2 - width;
	v2Start[1] = topY1; 
	SliceImageType::SizeType v2Size;
	v2Size[0] = 2*width;
	v2Size[1] = (topY2 - topY1) + space;
	SliceImageType::RegionType v2( v2Start, v2Size );	

	ImageIteratorType itv2( windowImageOut, v2 );
	NeighborhoodIteratorType outv2( gaussianOperatorH.GetRadius(), windowImageOut, v2 );
	
	outv2.GoToBegin();
	for (itv2.GoToBegin(); !itv2.IsAtEnd(); ++itv2, ++outv2)
	{
		itv2.Set( IP( outv2, gaussianOperatorH) );
	}
	//sitk::Show( sitk::Image( windowImageOut ), "windowImageOut");

	//v3 bot left
	SliceImageType::IndexType v3Start;
	v3Start[0] = botX1 - width;
	v3Start[1] = botY1 - space; 
	SliceImageType::SizeType v3Size;
	v3Size[0] = 2*width;
	v3Size[1] = (botY2 - botY1) + space;
	SliceImageType::RegionType v3( v3Start, v3Size );	

	ImageIteratorType itv3( windowImageOut, v3 );
	NeighborhoodIteratorType outv3( gaussianOperatorH.GetRadius(), windowImageOut, v3 );
	
	outv3.GoToBegin();
	for (itv3.GoToBegin(); !itv3.IsAtEnd(); ++itv3, ++outv3)
	{
		itv3.Set( IP( outv3, gaussianOperatorH) );
	}
	//sitk::Show( sitk::Image( windowImageOut ), "windowImageOut");

	//v4 bot right
	SliceImageType::IndexType v4Start;
	v4Start[0] = botX2 - width;
	v4Start[1] = botY1 - space; 
	SliceImageType::SizeType v4Size;
	v4Size[0] = 2*width;
	v4Size[1] = (botY2 - botY1) + space;
	SliceImageType::RegionType v4( v4Start, v4Size );	

	ImageIteratorType itv4( windowImageOut, v4 );
	NeighborhoodIteratorType outv4( gaussianOperatorH.GetRadius(), windowImageOut, v4 );
	
	outv4.GoToBegin();
	for (itv4.GoToBegin(); !itv4.IsAtEnd(); ++itv4, ++outv4)
	{
		itv4.Set( IP( outv4, gaussianOperatorH) );
	}
	//sitk::Show( sitk::Image( windowImageOut ), "windowImageOut");

	// h1 top left
	SliceImageType::IndexType h1Start;
	h1Start[0] = botX1 - space;
	h1Start[1] = botY1 - width; 
	SliceImageType::SizeType h1Size;
	h1Size[0] = (topX1 - botX1) + 2*space;
	h1Size[1] = 2*width;
	SliceImageType::RegionType h1( h1Start, h1Size );	

	ImageIteratorType ith1( windowImageOut, h1 );
	NeighborhoodIteratorType outh1( gaussianOperatorV.GetRadius(), windowImageOut, h1 );
	
	outh1.GoToBegin();
	for (ith1.GoToBegin(); !ith1.IsAtEnd(); ++ith1, ++outh1)
	{
		ith1.Set( IP( outh1, gaussianOperatorV) );
	}
	//sitk::Show( sitk::Image( windowImageOut ), "windowImageOut");

	// h1 top left
	SliceImageType::IndexType h2Start;
	h2Start[0] = topX2 - space;
	h2Start[1] = topY2 - width; 
	SliceImageType::SizeType h2Size;
	h2Size[0] = (botX2 - topX2) + 2*space;
	h2Size[1] = 2*width;
	SliceImageType::RegionType h2( h2Start, h2Size );	

	ImageIteratorType ith2( windowImageOut, h2 );
	NeighborhoodIteratorType outh2( gaussianOperatorV.GetRadius(), windowImageOut, h2 );
	
	outh2.GoToBegin();
	for (ith2.GoToBegin(); !ith2.IsAtEnd(); ++ith2, ++outh2)
	{
		ith2.Set( IP( outh2, gaussianOperatorV) );
	}
	//sitk::Show( sitk::Image( windowImageOut ), "windowImageOut");


	sitk::Show( sitk::Image( windowImageOut ), "windowImageOut");

  return EXIT_SUCCESS;
}






